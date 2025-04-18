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

.job_survn <- setClass("job_survn", contains = c("job_survival"))

setGeneric("do_survival",
  function(x, ref, ...) standardGeneric("do_survival"))

setMethod("do_survival", signature = c(x = "job_limma", ref = "job_survival"),
  function(x, ref, ...){
    cox <- asjob_lasso(x)
    cox <- step1(cox)
    surv <- asjob_survival(
      cox, fea_coefs = ref$fea_coefs, force = TRUE, ...
    )
    surv <- methodAdd(surv, x@meth$step1)
    surv
  })

setMethod("do_survival", signature = c(x = "list", ref = "job_survival"),
  function(x, ref, ...){
    projects <- vapply(x, function(x) x$project, character(1))
    inst <- x[[1]]
    x <- lapply(x, 
      function(x) {
        if (!is(x, "job_limma")) {
          stop('!is(x, "job_limma").')
        }
        cox <- asjob_lasso(x)
        cox@sig <- x$project
        cox
      })
    data <- x[[1]]
    for (i in seq_along(x[-1])) {
      message(glue::glue("Merge with {x[[i]]$project}"))
      data <- merge(data, x[[i]])
    }
    message("Finished merging, dim: ", bind(dim(object(data))))
    data <- step1(data)
    surv <- asjob_survival(
      data, fea_coefs = ref$fea_coefs, force = TRUE, ...
    )
    surv$objects <- lapply(x, 
      function(obj) {
        obj <- step1(obj)
        obj <- asjob_survival(
          obj, fea_coefs = ref$fea_coefs, force = TRUE, ...
        )
        obj
      })
    names(surv$objects) <- projects
    surv <- methodAdd(surv, inst@meth$step1)
    if (length(x) > 1) {
      surv <- snapAdd(surv, "合并数据集 ({bind(projects)})。{.merge_alias_method(ref, surv, TRUE)}随后，将基因表达数据归一化 (Z-score)。{surv@snap$step0}此外，对于未合并前的各个数据集，以相同的方式生存分析。")
    }
    return(.job_survn(surv))
  })

.merge_alias_method <- function(from, to, merge = FALSE) {
  glue::glue("{if (merge) '对于不同注释来源的基因名，以 `org.Hs.eg.db::org.Hs.eg.db` 获取基因的别名 (ALIAS) ，根据 (ALIAS) 的一致性合并。' else ''}查找预后模型中基因的 ALIAS，在未找到对应基因的情况下，使用该基因的 ALIAS 查找。 (原模型基因：{bind(from$fea_coefs$feature)}；以 ALIAS 匹配后，基因为：{bind(to$fea_coefs$feature)}) 。")
}

setGeneric("asjob_survival", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_survival"))

setMethod("asjob_survival", signature = c(x = "job_lasso"),
  function(x, use.group = c("mul_cox", "uni_cox"), lambda = c("min", "1se"),
    base_method = c("surv_cutpoint", "median"), fea_coefs = NULL, force = FALSE,
    use_data = c("all", "train", "valid"), 
    alias = NULL, use_alias = TRUE, db = org.Hs.eg.db::org.Hs.eg.db)
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
      if (use_alias) {
        merge_alias <- merge(fea_coefs$feature, colnames(data), db = db)
        fea_coefs$feature <- map(fea_coefs$feature, merge_alias, "x")
        colnames(data) <- map(colnames(data), merge_alias, "y")
        use_genes <- fea_coefs$feature
        not_cover <- use_genes[ !use_genes %in% colnames(data) ]
        message(glue::glue("Alias not cover: {crayon::red(bind(not_cover))}."))
      }
      if (!is.null(alias)) {
        if (!is(alias, "list")) {
          stop('!is(alias, "list").')
        }
        fun_reset <- function(db, from, to) {
          ifelse(db == from, to, db)
        }
        message("Try found alias genes.")
        lapply(not_cover, 
          function(name) {
            alias <- alias[[ name ]]
            isThat <- alias %in% colnames(data)
            if (!any(isThat)) {
              message(crayon::red(glue::glue("Alias of {name} neither covered.")))
            } else {
              use <- alias[ which(isThat)[1] ]
              message(glue::glue("Use alias: {use} ({name})"))
              fea_coefs$feature <- fun_reset(fea_coefs$feature, name, use)
              fea_coefs <<- fea_coefs
            }
          })
        use_genes <- fea_coefs$feature
      }
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
    snap <- x@snap$step1
    meth <- x@meth$step1
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
    x <- snapAdd(x, snap)
    x <- methodAdd(x, meth)
    snapAdd_onExit("x", "根据元数据信息 (即临床数据) ，去除了生存状态未知的样例。")
    x$time <- time
    x$status <- status
    x$genes_surv <- genes_surv
    x$fun_group <- fun_group
    x$fun_status <- function(x) ifelse(x == "Alive", 0, 1) 
    return(x)
  })

setMethod("regroup", signature = c(x = "job_limma", ref = "job_survival"),
  function(x, ref, feature = names(ref$data_group)[which], 
    which = 1L, ...)
  {
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
    data <- data[match(.get_meta(x, "sample"), data$sample), ]
    if (!identical(.get_meta(x, "sample"), data$sample)) {
      stop('!identical(.get_meta(x, "sample"), data$sample).')
    }
    message("Reset group.")
    x <- modify_job_limma_meta(
      x, group = data$group, fun = dplyr::mutate, modify_object = TRUE
    )
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

setMethod("step1", signature = c(x = "job_survn"),
  function(x, ...){
    x <- callNextMethod(x, ...)
    message("For individuals ...")
    x$objects <- lapply(x$objects, 
      function(obj) {
        step1(obj, ...)
      })
    return(x)
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
        coxfit <- survival::coxph(
          as.formula(glue::glue("survival::Surv(time, event = status) ~ `{genes}`")), data = data
        )
        ph <- survival::cox.zph(coxfit)
        fit <- survival::survfit.formula(survival::Surv(time, status) ~ group, data = data)
        diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = data)
        title <- paste0(genes, collapse = " & ")
        p.surv <- survminer::ggsurvplot(fit, data = data,
          pval = TRUE, conf.int = TRUE, risk.table = TRUE, title = title,
          risk.table.y.text.col = T,
          risk.table.y.text = FALSE,
          ggtheme = theme_bw()
        )
        p.surv$table <- p.surv$table + labs(x = "time (month)")
        p.surv <- wrap(grid.grabExpr(print(p.surv)), 5.5, 5)
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
          cols <- c("blue", "red", "orange", "green", "black", "pink", "purple")[seq_along(roc_time)]
          expr <- expression({
            roc <- timeROC::timeROC(
              data$time / 12, data$status, data[[ genes ]], 
              cause = 1, times = roc_time, weighting = "cox"
            )
            legends <- mapply(roc_time, cols, seq_along(roc_time),
              FUN = function(time, col, n) {
                plot(roc, time, col = col, add = (n != 1), title = "")
                glue::glue("AUC at {time} year: {round(roc[[ 'AUC' ]][ paste0('t=', time) ], 2)}")
              })
            legend(.5, .2, legends, text.col = cols, bty = "n")
          })
          p.roc <- wrap(as_grob(expr), 7, 7)
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
        namel(
          p.surv, fit, p.roc, roc, p.boxplot, diff, ph, data_group
        )
      })
    x$data_group <- lapply(lst, function(x) x[[ "data_group" ]])
    plots <- sapply(c("p.surv", "p.roc", "p.boxplot"), simplify = FALSE,
      function(name) {
        lapply(lst, function(x) x[[ name ]])
      })
    t.SurvivalPValue <- tibble::tibble(
      name = names(lst), 
      pvalue = vapply(lst, function(x) x$diff$pvalue, double(1)),
      ph = vapply(lst, function(x) x$ph$table[1, "p"], double(1)),
      group_low_survival = vapply(lst, function(x) .get_group_mortality.survdiff(x$diff), character(1))
    )
    t.SignificantSurvivalPValue <- dplyr::filter(
      t.SurvivalPValue, pvalue < .05, ph > .05
    )
    if (length(genes_surv) > 1) {
      x <- snapAdd(x, "根据 P value &lt; 0.05, 共筛到 {nrow(t.SignificantSurvivalPValue)} 个特征。
        分别为 {bind(t.SignificantSurvivalPValue$name)}。")
    } else if (length(genes_surv)) {
      x <- snapAdd(x, "{genes_surv} 生存分析 P value &lt; 0.05, 且满足 ph 假设检验 (P &gt; 0.05)。")
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
    if (is.null(x$fea_coefs)) {
      feature(x) <- t.SignificantSurvivalPValue$name
    } else {
      feature(x) <- dplyr::filter(x$fea_coefs, coef != 0)$feature
    }
    x <- tablesAdd(x, t.SurvivalPValue, t.SignificantSurvivalPValue)
    x <- methodAdd(x, "以 R 包 `survival` ({packageVersion('survival')}) 生存分析，以 R 包 `survminer` ({packageVersion('survminer')}) 绘制生存曲线。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_survival"),
  function(x){
    step_message("Collate results.")
    plots <- x@plots$step1$p.surv
    if (length(plots) <= 1) {
      stop('length(plots) <= 1.')
    }
    p.survs <- smart_wrap(plots, 3)
    p.survs <- set_lab_legend(
      p.survs,
      glue::glue("{x@sig} all significant genes survival curves"),
      glue::glue("为所有基因的生存曲线图 (显著的)。")
    )
    x <- plotsAdd(x, p.survs = p.survs)
    return(x)
  })

setMethod("step2", signature = c(x = "job_survn"),
  function(x){
    step_message("Collate results.")
    fun_collate <- function(name) {
      objects <- c(list(Merged = x), x$objects)
      plots <- lapply(objects, 
        function(obj) {
          obj@plots$step1[[name]][[1]]
        })
      wrap_layout(frame_wrap(plots), c(1, length(objects)), 5)
    }
    p.survs <- fun_collate("p.surv")
    p.survs <- set_lab_legend(
      p.survs,
      glue::glue("{x@sig} all datasets survival plot"),
      glue::glue("为所有数据集的生存分析图。")
    )
    p.rocs <- fun_collate("p.roc")
    p.rocs <- set_lab_legend(
      p.rocs,
      glue::glue("{x@sig} all datasets ROC validation"),
      glue::glue("为所有数据集的 ROC 图。")
    )
    x <- plotsAdd(x, p.survs, p.rocs)
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

