# ==========================================================================
# workflow of rms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_rms <- setClass("job_rms", 
  contains = c("job"),
  prototype = prototype(
    pg = "rms",
    info = c("https://cran.r-project.org/web/packages/rms/refman/rms.html"),
    cite = "",
    method = "",
    tag = "rms",
    analysis = "rms 列线图构建与评估"
    ))

setGeneric("asjob_rms",
  function(x, ...) standardGeneric("asjob_rms"))

setMethod("asjob_rms", signature = c(x = "job_deseq2"),
  function(x, ref, target = "group", format_name = TRUE, levels = rev(.guess_compare_deseq2(x, 1L)))
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    if (!is(ref, "feature")) {
      stop('!is(ref, "feature").')
    }
    ref <- resolve_feature_snapAdd_onExit("x", ref)
    data <- SummarizedExperiment::assay(x$vst)
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    data <- as_tibble(t(data[rownames(data) %in% ref, ]), idcol = "sample")
    metadata <- x$vst@colData
    data <- map(data, "sample", metadata, "sample", target, col = "group")
    if (format_name) {
      message(glue::glue("Colnames: {bind(colnames(data))}"))
      colnames(data) <- formal_name(colnames(data))
      message(glue::glue("Reset as: {bind(colnames(data))}"))
      target <- formal_name(target)
    }
    data[[ target ]] <- factor(data[[ target ]], levels = levels)
    x <- .job_rms(object = data)
    x$metadata <- metadata
    x$target <- target
    return(x)
  })

setMethod("step0", signature = c(x = "job_rms"),
  function(x){
    step_message("Prepare your data with function `job_rms`.")
  })

setMethod("step1", signature = c(x = "job_rms"),
  function(x, target = x$target, seed = x$seed, ...)
  {
    step_message("Logistic ...")
    data <- object(x)[, colnames(object(x)) != "sample" ]
    set.seed(seed)
    x$res_lrm <- new_lrm(data, formula = target)
    x$res_nomo <- new_nomo(x$res_lrm, ...)
    x <- methodAdd(
      x, "以 R 包 `rms` ({packageVersion('rms')}) 依据各独立诊断因子的权重赋值打分，将各因子评分求和得到总评分，作为风险评估模型，进而推断受试者的患病风险。"
    )
    p.nomo <- x$res_nomo$p.nomo_reg
    p.nomo <- set_lab_legend(
      p.nomo,
      glue::glue("{x@sig} nomogram"),
      glue::glue("风险评估列线图|||第一部分为 Points，表示风险分数为某个取值时的单项得分；第二部分为变量，变量后面线段的取值范围表示该变量对结局事件的总贡献值。线段上的刻度表示变量的不同取值；第三部分为 Total Points，表示变量取值的单项得分的总分。")
    )
    x <- methodAdd(x, "以 R 包 `regplot` ({packageVersion('regplot')}) 绘制该评估模型的列线图 (优化的列线图)。")
    x <- snapAdd(x, "列线图{aref(p.nomo)}将复杂的回归方程，转变为了可视化的图形，使预测模型的结果更具有可读性，方便对患者进行评估。如图{aref(p.nomo)}，每一个关键基因对应一个评分，各关键基因评分相加对应总评分，根据总评分预测疾病发病风险。")
    p.rocs <- x$res_lrm$p.rocs
    p.rocs <- set_lab_legend(
      p.rocs,
      glue::glue("{x@sig} ROC evaluation of diagnostic indicators"),
      glue::glue("各诊断指标的ROC|||纵轴为灵敏度，横轴为特异度，虚线为基准线（最低标准），曲线为对应指标的 ROC 曲线。其中 ROC 曲线距离基准线越远，则说明该模型的预测效果越好。ROC 曲线接近左上角，说明模型预测准确率很高。")
    )
    x <- methodAdd(x, "以 R 包 `pROC` ({packageVersion('pROC')}) 绘制该诊断模型的受试者工作特征 (ROC) 曲线，以曲线下面积 (Area Under the Curve，AUC) 评估模型效能。")
    x <- snapAdd(x, "ROC 评估{aref(p.rocs)}各诊断指标：{bind(x$res_lrm$aucs)} (AUC 介于 0.7-1 提示模型具有较好的预测效能)。")
    p.cal <- x$res_lrm$p.cal
    exLegend <- if (packageVersion("rms") >= "8.1.1") {
      "；C.L.为置信区间，如果对角线在 C.L. 内，模型校准良好"
    } else {
      ""
    }
    hl.pvalue <- x$res_lrm$hl.test$p.value
    p.cal <- set_lab_legend(
      p.cal,
      glue::glue("{x@sig} calibration curve"),
      glue::glue("列线图校准曲线|||Ideal 表示理想情况下预测概率与实际概率完全一致的情况，Apparent 表示模型预测概率的表观校准情况，Bias-corrected 表示经过偏差校正后的校准情况{exLegend}。Hosmer-Lemesho 检验 P = {fmt(hl.pvalue)}。")
    )
    x <- methodAdd(x, "以 R 包 `rms` 对列线图绘制校准曲线，评估模型预测风险与实际患病风险的一致性（校准曲线越贴近对角线，表明一致性越高。")
    # x <- snapAdd(x, "对于 rms 校准数据，{clinical_thresholds(x$res_lrm$cal)}")
    x <- snapAdd(x, "校准曲线{aref(p.cal)}斜率接近 1 且 P &gt; 0.05，说明列线图的预测准确性较好。")
    x <- snapAdd(x, "Hosmer-Lemesho 检验 P 为{fmt(hl.pvalue)} (P &gt; 0.05 说明通过 HL 检验，预测值与真实值之间并无非常明显的差异)。")
    x <- plotsAdd(x, p.nomo, p.rocs, p.cal)
    return(x)
  })

setMethod("step2", signature = c(x = "job_rms"),
  function(x, B = 500, ...){
    step_message("rmda.")
    data <- x$res_lrm$data
    markers <- colnames(data)[ colnames(data) != x$target ]
    markers <- c(as.list(markers), list("nomogram"))
    data[[x$target]] <- .check_events_for_factor(data[[x$target]])
    nomogram <- plogis(predict(x$res_lrm$fit, type = "lp"))
    data <- dplyr::mutate(data, nomogram = !!nomogram)
    cli::cli_h1("rmda::decision_curve")
    dcas <- lapply(markers,
      function(marker) {
        rmda::decision_curve(
          as.formula(paste0(x$target, "~", marker)),
          data, study.design = "case", bootstrap = B, ...
        )
      })
    p.dcas <- funPlot(rmda::plot_decision_curve,
      list(x = dcas, curve.names = unlist(markers), confidence.intervals = FALSE)
    )
    p.dcas <- set_lab_legend(
      p.dcas,
      glue::glue("{x@sig} Decision curve analysis"),
      glue::glue("决策曲线分析 (Decision curve analysis)|||DCA 用于评估不同模型在不同阈值下的净收益。横坐标为风险阈值，纵坐标为净获益率，平行于 x 轴的虚线 None 是不对任何人进行干预，抛物线形状的虚线 All 是对所有人进行干预，实线代表各指标的干预效果。")
    )
    x <- methodAdd(x, "以 R 包 `rmda` ({packageVersion('rmda')}) 进行决策曲线分析 (Decision curve analysis，DCA) 并绘制 DCA 曲线。")
    x <- snapAdd(x, "DCA 分析{{aref(p.dcas)}}可知，Nomogram 在阈值范围内的净收益高于 All 和 None 策略，同时高于其他独立诊断因子，显示出其在预测中的潜在价值。")
    x <- plotsAdd(x, p.dcas)
    return(x)
  })

clinical_thresholds <- function(cal, max_deviation = 0.1) {
  cal_data <- data.frame(
    mean_pred = cal[, "predy"],
    mean_obs = cal[, "calibrated.corrected"]
  )
  deviation <- abs(cal_data$mean_pred - cal_data$mean_obs)
  if (all(deviation < max_deviation)) {
    paste0("模型校准良好：所有风险层的预测误差都在", max_deviation*100, "%以内。")
  } else if (mean(deviation) < max_deviation) {
    message("The following risk layers have significant deviations:")
    bad_layers <- which(deviation >= max_deviation)
    for (i in bad_layers) {
      message(" Risk layer", i, ": prediction ", round(cal_data$mean_pred[i]*100, 1), 
        "%, actual ", round(cal_data$mean_obs[i]*100, 1), 
        "%, deviation ", round(deviation[i]*100, 1), "%\n")
    }
    paste0(
      "平均预测误差 (", signif(mean(deviation), 2)*100, "%) 在",
      max_deviation*100, "%以内。"
    )
  } else {
    paste0("模型校准不佳：平均预测误差为", round(mean(deviation)*100, 1), 
      "%, 超过", max_deviation*100, "%的可接受范围。")
  }
}

## logistic
new_lrm <- function(data, formula, rev.level = FALSE, 
  lang = c("en", "cn"), run_boot = FALSE, B = 500, ...)
{
  fun_escape_bug_of_rms <- function() {
    wh <- which(vapply(seq_len(ncol(data)), function(n) is.factor(data[[ n ]]), FUN.VALUE = logical(1)))
    if (any(grpl(colnames(data)[ wh ], "\\s"))) {
      stop("`Due to the bug of `rms`, the names of columns of which is factor, can not contains blank")
    }
  }
  fun_escape_bug_of_rms()
  if (is.character(formula)) {
    message("Detected `formula` input with 'character', use as 'y'.")
    if (length(formula) > 1) {
      stop("length(formula) > 1")
    }
    if (!formula %in% colnames(data)) {
      stop("The 'y' not found in colnames of data.")
    }
    formula <- paste0(formula, " ~ ",
      paste0(
        paste0("`", colnames(data)[ colnames(data) != formula ], "`"),
        collapse = " + "
      )
    )
    message("Guess Formula: ", formula)
    formula <- as.formula(formula)
  }
  isColChar <- apply(data, 2, is.character)
  if (any(isColChar)) {
    message("Convert all character columns as factor.")
    data <- dplyr::mutate(data, dplyr::across(dplyr::where(is.character), as.factor))
  }
  lang <- match.arg(lang)
  y <- as.character(formula[[ 2 ]])
  outcome <- data[[ y ]]
  if (rev.level) {
    data[[ y ]] <- factor(outcome, levels = rev(levels(outcome)))
  }
  message("Check levels: ", paste0(levels <- levels(data[[ y ]]), collapse = ", "))
  set_rms_datadist(data)
  fit <- e(rms::lrm(formula, data = data, x = TRUE, y = TRUE))
  # fit_glm <- stats::glm(formula, data = data, family = stats::binomial())
  # 95\\% CL
  # exp(confint.default(lrm.eff$fit))
  if (run_boot) {
    message("Use boot to calculate average C-index ...")
    ## use boot to calculate average C-index and 95\\% CI
    fun_c_index <- function(data, indices) {
      data <- data[ indices, ]
      fit <- rms::lrm(formula, data = data)
      fit$stats[ c("C", "P") ]
    }
    boot <- boot::boot(data, fun_c_index, R = B)
    boot.ci <- boot::boot.ci(boot, .95, type = "basic")
    fun <- function() {
      means <- apply(boot$t, 2, mean)
      ci <- tail(boot.ci$basic[1, ], n = 2)
      new_lich(list(`Re-sample` = B, `C-index` = means[1],
          `P-value` = means[2], "95\\% CI" = paste0(ci, collapse = " ~ "))
      )
    }
    lich <- fun()
    lich <- .set_lab(lich, "bootstrap others")
    boots <- namel(boot, boot.ci, lich)
  } else {
    boots <- list()
  }
  if (FALSE) {
    old <- getOption("prType")
    options(prType = "html")
    html <- paste0(print(lrm.supp$fit))
    coefs <- get_table.html(html)
    options(prType = old)
  } else {
    coefs <- NULL
  }
  cal <- e(rms::calibrate(fit, method = "boot", B = B))
  if (lang == "cn") {
    xlab <- paste0("预测", levels[2], "概率")
    ylab <- paste0("实际", levels[2], "概率")
  } else {
    xlab <- "Predicted Probability"
    ylab <- "observed Probability"
  }
  # rms::plot.calibrate
  # p <- predict(fit, type = "fitted")
  # y <- fit$y
  # ResourceSelection::hoslem.test(y, p, g = 10)
  fun_plot_calibrate_with_hl <- function(..., p)
  {
    rms:::plot.calibrate.default(...)
    text(
      x = .1, y = .8, labels = glue::glue("Hosmer-Lemeshow P = {fmt(p)}"), cex = 1.2, adj = 0
    )
  }
  hl.test <- e(
    ResourceSelection::hoslem.test(.check_events_for_factor(fit$y), predict(fit, type = "fitted"))
  )
  p.cal <- funPlot(fun_plot_calibrate_with_hl,
    list(
      x = cal, xlim = c(0, 1), ylim = c(0, 1),
      xlab = xlab, ylab = ylab, riskdist = TRUE,
      subtitle = TRUE, p = hl.test$p.value
    )
  )
  if (lang == "cn") {
    message("Try convert English legend as Chinese.")
    labels <- p.cal$children[[15]]$label
    if (identical(labels, c("Apparent", "Bias-corrected", "Ideal"))) {
      p.cal$children[[15]]$label <- c("表观状态", "误差纠正", "理想状态")
    }
  }
  roc <- new_roc(fit$y, predict(fit), lang = lang, ...)
  coefs <- names(fit$coefficients)[ names(fit$coefficients) != "Intercept" ]
  rocs <- lapply(coefs,
    function(coef) {
      pROC::roc(data[[ y ]], data[[ coef ]], ci = TRUE)
    })
  rocs <- setNames(c(rocs, list(roc$roc)), c(coefs, "Nomogram"))
  p.rocs <- plot_roc(rocs)
  aucs <- lapply(rocs, function(roc) signif(roc$auc[[1]], 3))
  # rocs <- 
  .add_internal_job(.job(method = "R package `rms` used for Logistic regression and nomogram visualization"))
  namel(
    fit, coefs, data, levels, cal, p.cal, hl.test, roc, p.rocs, 
    aucs, lang, boots
  )
}

set_rms_datadist <- function(data) {
  .RMS_datadist <- e(rms::datadist(data))
  assign(".RMS_datadist", .RMS_datadist, envir = .GlobalEnv)
  message("A global variable defined herein: `.RMS_datadist`.")
  options(datadist = ".RMS_datadist")
}

new_nomo <- function(lrm, fun_label = lrm$levels[2], lang = lrm$lang,
  lp = FALSE, fun_at = seq(.1, .9, by = .2))
{
  if (!is(lrm, "list")) {
    stop("the `lrm` should be object 'list' return by `new_lrm`")
  }
  lang <- match.arg(lang, c("cn", "en"))
  set_rms_datadist(lrm$data)
  fun_label <- if (lang == "en") {
    paste0("Risk of ", fun_label)
  } else {
    paste0(fun_label, "风险")
  }
  nomo <- e(rms::nomogram(lrm$fit, fun = stats::plogis,
      funlabel = fun_label, lp = lp, fun.at = fun_at))
  if (lang == "en") {
    args <- list(x = nomo,
      lplabel = "Linear Predictor",
      points.label = 'Points', total.points.label = 'Total Points'
    )
  } else {
    args <- list(x = nomo, lplabel = "线性预测",
      points.label = '分数', total.points.label = '总分'
    )
  }
  p.nomo <- wrap(
    funPlot(rms:::plot.nomogram, args), 8, .6 * length(lrm$fit$coefficients) + 2
  )
  p.nomo <- .set_lab(p.nomo, "nomogram plot")
  p.nomo_reg <- wrap(funPlot(regplot::regplot,
    list(
      reg = lrm$fit,
      interval = "confidence",
      observation = TRUE,
      points = TRUE,
      prfail = TRUE,
      title = paste(lrm$levels[2], "Diagnosis"),
      showP = TRUE
    )
  ))
  namel(p.nomo, p.nomo_reg, nomo)
}

new_roc <- function(y, x, ..., plot.thres = NULL, lang = c("en", "cn"), cn.mode = c("zhen", "1-"))
{
  lang <- match.arg(lang)
  roc <- pROC::roc(y, x, ..., ci = TRUE)
  if (TRUE) {
    # try get p-value
    # https://stackoverflow.com/questions/61997453/how-to-get-p-value-after-roc-analysis-with-proc-package
    fun_p <- function() {
      v <- pROC::var(roc)
      b <- roc$auc - .5
      se <- sqrt(v)
      z <- (b / se)
      2 * pt(-abs(z), df = Inf)
    }
    p.value <- fun_p()
  }
  thres <- pROC::coords(roc, "best")
  if (lang == "cn") {
    cn.mode <- match.arg(cn.mode)
    if (cn.mode == "zhen") {
      xlab <- "假阳性率"
      ylab <- "真阳性率"
    } else {
      xlab <- "1-特异性"
      ylab <- "敏感性"
    }
  } else {
    xlab <- "Specificity"
    ylab <- "Sensitivity"
  }
  p.roc <- as_grob(
    expression(pROC::plot.roc(roc, print.thres = plot.thres, print.auc = TRUE,
        print.auc.x = .4, print.auc.y = .05,
        xlab = xlab, ylab = ylab)), environment()
  )
  p.roc <- wrap(p.roc, 7, 7)
  p.roc <- .set_lab(p.roc, "ROC")
  .add_internal_job(.job(method = "R package `pROC` used for building ROC curve"))
  lich <- new_lich(
    list(AUC = as.double(roc$auc),
      "95\\% CI" = roc$ci[-2],
      `P-value` = p.value
    )
  )
  lich <- .set_lab(lich, "ROC others")
  data <- tibble::tibble(Sensitivities = roc$sensitivities, Specificities = roc$specificities)
  namel(p.roc, thres, roc, p.value, lich, data)
}

.check_events_for_factor <- function(events) {
  if (is.factor(events)) {
    if (length(levels(events)) != 2) {
      stop('length(levels(events)) != 2.')
    }
    events <- as.integer(events) - 1L
  }
  events
}

fmt <- function(x) formatC(x, digits = 4, format = "fg")

