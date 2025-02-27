# ==========================================================================
# workflow of survival
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_survival <- setClass("job_survival", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    method = "R package `survival` and `survminer` used for survival analysis",
    tag = "surv",
    analysis = "Survival 生存分析"
    ))

setGeneric("asjob_survival", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_survival"))

setMethod("asjob_survival", signature = c(x = "job_lasso"),
  function(x, use.group = c("mul_cox", "uni_cox"), lambda = c("min", "1se"),
    base_method = c("surv_cutpoint", "median"), fea_coefs = NULL, force = FALSE,
    use_data = c("all", "train", "valid"))
  {
    y <- .job_survival()
    use_data <- match.arg(use_data)
    base_method <- match.arg(base_method)
    if (x@step < 3L && !force) {
      stop("x@step < 3L")
    }
    use.group <- match.arg(use.group)
    message("Check coeffients table.")
    if (missing(fea_coefs)) {
      if (use.group == "mul_cox") {
        lambda <- paste0("lambda.", match.arg(lambda))
        fea_coefs <- x$sig.mul_cox
        fea_coefs <- dplyr::rename(fea_coefs, coef = !!rlang::sym(paste0("coef.", lambda)))
        fea_coefs <- dplyr::filter(fea_coefs, coef != 0)
        y <- snapAdd(y, "选择 {lambda} 时得到的特征集，包含 {nrow(fea_coefs)} 个基因，
          分别为: {bind(fea_coefs$feature)}。")
      }
      if (use.group == "uni_cox") {
        fea_coefs <- x$sig.uni_cox
        y <- snapAdd(y, "选择 Univariate COX 回归得到的显著特征集，包含 {nrow{fea_coefs}} 个基因。")
      }
      y <- snapAdd(y, "以回归系数构建风险评分模型。\n\n$$ Score = \\sum(expr(Gene) \\times coef) $$\n\n")
    } else {
      y <- snapAdd(y, "将风险评分模型运用于数据集 (dataset: {x@sig})。")
    }
    if (is.null(fea_coefs)) {
      stop("`fea_coefs` should be specified.")
    }
    data <- as_tibble(x$get(use_data)$data)
    if (is.null(fea_coefs$coef)) {
      message("Check `fea_coefs`.")
      fea_coefs$coef <- fea_coefs[[ grp(
        colnames(fea_coefs), lambda
        ) ]]
      needLess <- fea_coefs$feature[ fea_coefs$coef == 0 ]
      data <- data[, !colnames(data) %in% needLess]
      fea_coefs <- dplyr::filter(fea_coefs, !feature %in% needLess)
    }
    message("Check genes.")
    use_genes <- fea_coefs$feature
    if (any(!use_genes %in% colnames(data))) {
      not_cover <- use_genes[ !use_genes %in% colnames(data) ]
      message(glue::glue("Genes not cover: {crayon::red(bind(not_cover))}."))
      x$not_cover <- not_cover
      use_genes <- use_genes[ use_genes %in% colnames(data) ]
      fea_coefs <- dplyr::filter(fea_coefs, feature %in% use_genes)
      y <- snapAdd(y, "该数据集 (dataset: {x@sig}) 不包含基因：{bind(not_cover)}，后续计算评分时去除。")
    }
    data <- dplyr::select(data, 1, dplyr::all_of(use_genes))
    if (!identical(colnames(data)[-1], use_genes)) {
      stop("Gene names not match.")
    }
    data <- dplyr::mutate(data,
      risk_score = apply(data[, -1, drop = FALSE], 1,
        function(exprs) {
          sum(fea_coefs$coef * exprs)
        })
    )
    message("Set cutoff point.")
    if (base_method == "median") {
      message("Use Median cutoff point.")
      data <- dplyr::mutate(data,
        score_group = ifelse(risk_score > median(risk_score), "High", "Low")
      )
      cutoff <- median(data$risk_score)
      y <- snapAdd(y, "按中位风险评分，将样本分为 Low 和 High 风险组 (cutoff: {cutoff})
        ({try_snap(data$score_group)})， 随后进行生存分析。")
    } else if (base_method == "surv_cutpoint") {
      cutoff <- survminer::surv_cutpoint(
        data.frame(risk_score = data$risk_score,
          time = x$get(use_data)$time, event = x$get(use_data)$event),
        variables = "risk_score"
        )$cutpoint$cutpoint
      data <- dplyr::mutate(data,
        score_group = ifelse(risk_score > cutoff, "High", "Low")
      )
      y <- snapAdd(y, "按 `survminer::surv_cutpoint` 计算的 cutoff，
        将样本分为 Low 和 High 风险组 (cutoff: {cutoff})
        ({try_snap(data$score_group)})， 随后进行生存分析。")
    }
    message("Be careful, follow_up time should be 'day' as unit.")
    data$days_to_last_follow_up <- x$get(use_data)$time
    data$vital_status <- x$get(use_data)$types
    message("The heatmap plot was used scaled data.")
    p.surv_genes_hp <- plot_genes_heatmap(
      t(scale(x@object))[colnames(x@object) %in% use_genes, ],
      dplyr::select(data, sample = rownames, group = score_group)
    )
    p.surv_genes_hp <- .set_lab(wrap(p.surv_genes_hp, 9, 6), sig(x),
      "risk score related genes heatmap")
    object(y) <- data
    y$p.surv_genes_hp <- p.surv_genes_hp
    y$fun_group <- function(score, ...) ifelse(score > y$cutoff, "High", "Low")
    y$genes_surv <- "risk_score"
    y$fun_status <- function(x) ifelse(x == "Alive", 0, 1) 
    y$time <- "days_to_last_follow_up"
    y$status <- "vital_status"
    y$cutoff <- cutoff
    y$fea_coefs <- fea_coefs
    return(y)
  })

setMethod("asjob_survival", signature = c(x = "job_limma"),
  function(x, genes_surv, use = .guess_symbol(x),
    time = "days_to_last_follow_up", scale = TRUE,
    status = "vital_status", base_method = c("surv_cutpoint", "median"))
  {
    project <- x$project
    if (x@step < 1)
      stop("x@step < 1")
    if (x@step >= 1) {
      object(x) <- x@params$normed_data
    }
    if (is(object(x), "DGEList")) {
      object(x) <- as_EList.DGEList(object(x))
    }
    genes_surv <- resolve_feature_snapAdd_onExit("x", genes_surv)
    i.pos <- gname(object(x)$genes[[ use ]]) %in% genes_surv
    j.pos <- !is.na(object(x)$targets[[ status ]]) & !grpl(object(x)$targets[[ status ]], "not reported", TRUE)
    object(x) <- e(limma::`[.EList`(object(x), i.pos, j.pos))
    object <- object(x)
    counts <- object$E
    if (any(duplicated(names <- object$genes[[ use ]]))) {
      if (TRUE) {
        message("As name duplicated, distinct these..")
        counts <- counts[!duplicated(names), , drop = FALSE]
        names <- names[!duplicated(names)]
      } else {
        message("As name duplicated, mutate with number.")
        whichDup <- duplicated(names)
        names[ whichDup ] <- paste0(
          names[ whichDup ], "_", seq_along(names[ whichDup ])
        )
      }
    }
    rownames(counts) <- names
    data <- t(counts)
    if (scale) {
      data <- scale(data)
    }
    data <- as_tibble(data)
    metadata <- tibble::as_tibble(object$targets)
    if (any(colnames(metadata) == "rownames")) {
      metadata <- dplyr::select(metadata, -rownames)
    }
    data <- dplyr::bind_cols(data, metadata)
    data[[ time ]] <- .mutate_last_follow_up_time(data[[ time ]])
    data <- dplyr::relocate(data, 1, !!rlang::sym(time), !!rlang::sym(status))
    base_method <- match.arg(base_method)
    if (base_method == "surv_cutpoint") {
      message("Use `survminer::surv_cutpoint` to found cutoff.")
      fun_group <- function(x, time, status) {
        data <- data.frame(gene = x, time = time, event = status)
        cutoff <- survminer::surv_cutpoint(data, variables = "gene")$cutpoint$cutpoint
        group <- ifelse(x > cutoff, "High", "Low")
        attr(group, "cutoff") <- cutoff
        group
      }
      snapAdd_onExit("x", "按 `survminer::surv_cutpoint` 计算的 cutoff，将样本分为 Low 和 High 风险组。")
    } else if (base_method == "median") {
      fun_group <- function(x, ...) {
        cutoff <- median(x)
        group <- ifelse(x > cutoff, "High", "Low")
        attr(group, "cutoff") <- cutoff
        group
      }
      snapAdd_onExit("x", "按中位风险评分，将样本分为 Low 和 High 风险组。")
    }
    if (length(meth(x))) {
      methodAdd_onExit("x", meth(x)[[1]])
      methodAdd_onExit("x", "使用标准化过的基因表达数据。")
      snapAdd_onExit("x", "生存数据为{project}，使用标准化过的基因表达数据。")
    } else {
      methodAdd_onExit("x", "生存数据为{project}，使用该基因表达数据。")
      snapAdd_onExit("x", "生存数据为{project}，使用该基因表达数据。")
    }
    x <- .job_survival(object = data)
    snapAdd_onExit("x", "根据元数据信息 (即临床数据) ，去除了生存状态未知的样例。")
    x$time <- time
    x$status <- status
    x$genes_surv <- genes_surv
    x$fun_group <- fun_group
    x$fun_status <- function(x) ifelse(x == "Alive", 0, 1) 
    return(x)
  })

setMethod("regroup", signature = c(x = "job_limma", ref = "job_survival"),
  function(x, ref, feature = names(ref$data_group), ...){
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    if (x@step < 1L || x@step > 1L) {
      stop('x@step < 1L || x@step > 1L.')
    }
    if (length(feature) > 1) {
      stop('length(feature) > 1.')
    }
    message(glue::glue("Use {feature} (survival related expression) to regroup samples."))
    if (is.null(ref$data_group)) {
      stop('is.null(ref$data_group).')
    }
    if (!any(names(ref$data_group) == feature)) {
      stop('!any(names(ref$data_group) == feature).')
    }
    data <- ref$data_group[[ feature ]]
    x <- filter(x, sample %in% !!data$sample, type = "metadata", add_snap = FALSE)
    if (identical(.get_meta(x, "sample"), data$sample)) {
      message("Reset group.")
      x <- modify_job_limma_meta(
        x, group = data$group, fun = dplyr::mutate, modify_object = TRUE
      )
    } else {
      stop('identical(.get_meta(x, "sample"), data$sample) == FALSE')
    }
    x <- set_independence(x)
    x <- set_design(x, ...)
    des <- glue::glue("使用 Survival 分析 ({ref@sig}) 时定义的分组数据 (根据 {feature} 的表达，分为 {try_snap(.get_meta(x, 'group'))})。")
    x@meth[[ "step1" ]] <- des
    x@snap[[ "step1" ]] <- des
    return(x)
  })

as_EList.DGEList <- function(object) {
  new("EList", list(E = object$count,
      targets = object$samples, genes = object$genes))
}

setMethod("step0", signature = c(x = "job_survival"),
  function(x){
    step_message("Prepare your data with function `asjob_survival`.
      In general, job convert from 'job_tcga' to 'job_limma' were
      adapted to this workflow.
      "
    )
  })

setMethod("step1", signature = c(x = "job_survival"),
  function(x, genes_surv = x$genes_surv, fun_group = x$fun_group,
    fun_status = x$fun_status, time = x$time, status = x$status, only_keep_sig = "guess",
    roc_time = c(1, 3, 5))
  {
    step_message("Survival test.")
    data <- dplyr::rename(object(x), time = !!rlang::sym(time),
      o_status = !!rlang::sym(status))
    if (time == "days_to_last_follow_up" || time == "days_to_death") {
      p.density_follow_time <- ggplot(data) +
        geom_density(aes(x = time / 365)) +
        labs(x = "Time (year)")
      p.density_follow_time <- .set_lab(p.density_follow_time, sig(x), "density of follow up time")
      x$p.density_follow_time <- setLegend(
        wrap(p.density_follow_time, 6, 2.5), 
        "为随访时间整体分布趋势 (用于确定Time ROC 时间点) 。"
      )
      data <- dplyr::mutate(data, time = time / 30)
    }
    cli::cli_alert_info("survival::survfit")
    lst <- pbapply::pbsapply(genes_surv, simplify = FALSE,
      function(genes) {
        data <- dplyr::select(
          data, sample = rownames, time, o_status, !!!rlang::syms(genes)
        )
        data_group <- data <- dplyr::mutate(
          data, status = fun_status(o_status),
          group = fun_group(!!!rlang::syms(genes), time, status),
        )
        fit <- survival::survfit.formula(survival::Surv(time, status) ~ group, data = data)
        diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = data)
        title <- paste0(genes, collapse = " & ")
        p.surv <- survminer::ggsurvplot(fit, data = data,
          pval = TRUE, conf.int = TRUE, risk.table = TRUE, title = title,
          ggtheme = theme_bw()
        )
        p.surv$table <- p.surv$table + labs(x = "time (month)")
        p.surv <- wrap(p.surv, 5.5, 5)
        p.surv <- .set_lab(p.surv, sig(x), "survival curve of", title)
        p.surv <- setLegend(p.surv, "为 {title} 生存曲线。")
        # timeROC
        if (length(genes) == 1) {
          require(survival)
          if (any(roc_time == 1) && !any(dplyr::filter(data, time / 12 <= 1)$o_status == "Dead")) {
            message("No event of 'Dead' in year 1, switch to 3, 5, 10.")
            x <<- snapAdd(x, "不存在 1 年以内 Dead 的样本，因此以 3, 5, 10 计算 ROC。")
            roc_time <- c(3, 5, 10)
          }
          cols <- c("blue", "red", "orange")[seq_along(roc_time)]
          roc <- timeROC::timeROC(
            data$time / 12, data$status, data[[ genes ]], 
            cause = 1, times = roc_time, weighting = "cox"
          )
          legends <- mapply(roc_time, cols, seq_along(roc_time),
            FUN = function(time, col, n) {
              plot(roc, time, col = col, add = (n != 1), title = "")
              glue::glue("AUC at {time} year: {round(roc[[ 'AUC' ]][ paste0('t=', time) ], 2)}")
            })
          legend("bottomright", legends, text.col = cols, bty = "n")
          p.roc <- wrap(recordPlot(), 7, 7)
          p.roc <- .set_lab(p.roc, sig(x), "time ROC")
          p.roc <- setLegend(p.roc, "为 {genes} {bind(roc_time)} 年生存分析 ROC 曲线。")
          data <- dplyr::mutate(data,
            var = ifelse(time / 12 <= 1, "year 1",
              ifelse(time / 12 <= 3, "year 3", "year 5")),
            group = as.factor(o_status), value = !!rlang::sym(genes)
          )
          p.boxplot <- try(wrap(.map_boxplot2(data, TRUE, ylab = genes), 6, 6), TRUE)
          if (!inherits(p.boxplot, "try-error")) {
            p.boxplot <- .set_lab(p.boxplot, sig(x), glue::glue("boxplot of {genes}"))
          } else {
            p.boxplot <- NULL
          }
        } else {
          p.roc <- NULL
          roc <- NULL
          p.boxplot <- NULL
        }
        namel(p.surv, fit, p.roc, roc, p.boxplot, diff, data_group)
      })
    x$data_group <- lapply(lst, function(x) x[[ "data_group" ]])
    plots <- sapply(c("p.surv", "p.roc", "p.boxplot"), simplify = FALSE,
      function(name) {
        lapply(lst, function(x) x[[ name ]])
      })
    t.SurvivalPValue <- tibble::tibble(
      name = names(lst), 
      pvalue = vapply(lst, function(x) x$diff$pvalue, double(1)),
      group_low_survival = vapply(lst, function(x) .get_group_mortality.survdiff(x$diff), character(1))
    )
    t.SignificantSurvivalPValue <- dplyr::filter(t.SurvivalPValue, pvalue < .05)
    if (length(genes_surv) > 1) {
      x <- snapAdd(x, "根据 P value &lt; 0.05, 共筛到 {nrow(t.SignificantSurvivalPValue)} 个特征。
        分别为 {bind(t.SignificantSurvivalPValue$name)}。")
    }
    x$models <- lapply(lst, function(x) list(fit = x$fit, roc = x$roc, diff = x$diff))
    if (identical(only_keep_sig, "guess")) {
      only_keep_sig <- length(genes_surv) > 5
    }
    if (only_keep_sig) {
      plots <- lapply(plots,
        function(x) {
          x <- x[ names(x) %in% t.SignificantSurvivalPValue$name ]
          if (!length(x)) {
            NULL
          } else x
        })
      x$models <- x$models[ names(x$models) %in% t.SignificantSurvivalPValue$name ]
    }
    if (length(plots)) {
      plots <- .set_lab(plots, sig(x), c("Survival", "ROC", "Boxplot"), "plots")
      x@plots[[1]] <- plots
    }
    feature(x) <- t.SignificantSurvivalPValue$name
    x <- tablesAdd(x, t.SurvivalPValue, t.SignificantSurvivalPValue)
    x <- methodAdd(x, "以 R 包 `survival` ({packageVersion('survival')}) 生存分析，以 R 包 `survminer` ({packageVersion('survminer')}) 绘制生存曲线。")
    return(x)
  })

collate_dataset_survival <- function(x, name = "survival",
  fun_extract = function(x) x@tables$step1$t.SignificantSurvivalPValue, exclude = NULL, ...)
{
  data <- collate(x, fun_extract, exclude, ...)
  data <- rbind_list(data, .id = "Dataset")
  if (!is.null(name)) {
    if (!any(colnames(data) == "name")) {
      stop('!any(colnames(data) == "name").')
    }
    data <- set_lab_legend(
      data, glue::glue("Genes of {name} in datasets"),
      glue::glue("为基因集 ({less(unique(data$name))}) 在多个数据集中的 {name} 统计。")
    )
    attr(data, "__COLLATE_NAME__") <- name
  }
  return(data)
}


.get_group_mortality.survdiff <- function(object, value = FALSE) {
  if (!is(object, "survdiff")) {
    stop('!is(object, "survdiff").')
  }
  groups <- s(names(object$n), "group=", "")
  mortality <- object$obs / object$n
  if (value) {
    mortality
  } else {
    if (mortality[1] > mortality[2]) groups[1] else groups[2]
  }
}

