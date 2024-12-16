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

setGeneric("asjob_survival", 
  function(x, ...) standardGeneric("asjob_survival"))

setMethod("asjob_survival", signature = c(x = "job_lasso"),
  function(x, use.group = "uni_cox", base_method = median, fea_coefs = NULL, force = F)
  {
    if (x@step < 3L && !force) {
      stop("x@step < 3L")
    }
    use.group <- match.arg(use.group)
    if (missing(fea_coefs) && use.group == "uni_cox") {
      fea_coefs <- x$sig.uni_cox
      if (is.null(fea_coefs)) {
        stop("`sig.uni_cox` not found.")
      }
    }
    use_genes <- fea_coefs$feature
    data <- as_tibble(x@object)
    if (any(!use_genes %in% colnames(data))) {
      not_cover <- use_genes[ !use_genes %in% colnames(data) ]
      message(glue::glue("Genes not cover: {crayon::red(bind(not_cover))}."))
      x$not_cover <- not_cover
      use_genes <- use_genes[ use_genes %in% colnames(data) ]
      fea_coefs <- dplyr::filter(fea_coefs, feature %in% use_genes)
    }
    data <- dplyr::select(data, 1, dplyr::all_of(use_genes))
    if (!identical(colnames(data)[-1], use_genes)) {
      stop("Gene names not match.")
    }
    data <- dplyr::mutate(data,
      risk_score = apply(data[, -1, drop = F], 1,
        function(exprs) {
          sum(fea_coefs$coef * exprs)
        }),
      score_group = ifelse(risk_score > base_method(risk_score), "High", "Low")
    )
    data$days_to_last_follow_up <- x$time
    data$vital_status <- x$target
    message("The heatmap plot was used scaled data.")
    p.surv_genes_hp <- plot_genes_heatmap(
      t(scale(x@object))[colnames(x@object) %in% use_genes, ],
      dplyr::select(data, sample = rownames, group = score_group)
    )
    p.surv_genes_hp <- .set_lab(wrap(p.surv_genes_hp, 9, 6), sig(x),
      "risk score related genes heatmap")
    x <- .job_survival(object = data)
    x$p.surv_genes_hp <- p.surv_genes_hp
    x$fun_group <- function(score) ifelse(score > x$cutoff, "High", "Low")
    x$genes_surv <- "risk_score"
    x$fun_status <- sur_status
    x$time <- "days_to_last_follow_up"
    x$status <- "vital_status"
    x$cutoff <- base_method(data$risk_score)
    meth(x)$step0 <- glue::glue("将 Univariate COX 回归系数用于风险评分计算，根据中位风险评分 {x$cutoff} 将患者分为低危组和高危组。")
    return(x)
  })

setMethod("asjob_survival", signature = c(x = "job_limma"),
  function(x, genes_surv, use = "gene_name",
    time = "days_to_last_follow_up", status = "vital_status")
  {
    if (x@step < 1)
      stop("x@step < 1")
    if (x@step > 1) {
      object(x) <- x@params$normed_data
    }
    i.pos <- object(x)$genes[[ use ]] %in% genes_surv
    j.pos <- !is.na(object(x)$targets[[ status ]]) & !grpl(object(x)$targets[[ status ]], "not reported", T)
    object(x) <- e(limma::`[.EList`(object(x), i.pos, j.pos))
    object <- object(x)
    counts <- object$E
    rownames(counts) <- object$genes[[ use ]]
    data <- as_tibble(t(counts))
    data <- dplyr::bind_cols(data, select(as_tibble(object$targets), -1))
    if (any(is.na(data[[time]]))) {
      message(glue::glue("NA in {time} found, as max"))
      fun_time_mutate <- function(x) {
        max <- max(x, na.rm = T)
        ifelse(is.na(x), max, x)
      }
      data <- dplyr::mutate(data, dplyr::across(!!rlang::sym(time), fun_time_mutate))
    }
    data <- relocate(data, 1, !!rlang::sym(time), !!rlang::sym(status))
    x <- .job_survival(object = data)
    x$time <- time
    x$status <- status
    x$genes_surv <- genes_surv
    x$fun_group <- sur_group_median
    x$fun_status <- sur_status
    x <- methodAdd(x, "去除了生存状态未知的数据。")
    return(x)
  })

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
    fun_status = x$fun_status, time = x$time, status = x$status, only_keep_sig = T,
    roc_time = c(1, 3, 5))
  {
    step_message("Survival test.")
    data <- dplyr::rename(object(x), time = !!rlang::sym(time),
      o_status = !!rlang::sym(status))
    if (time == "days_to_last_follow_up" || time == "days_to_death") {
      data <- dplyr::mutate(data, time = time / 30)
    }
    cli::cli_alert_info("survival::survfit")
    lst <- pbapply::pbsapply(genes_surv, simplify = F,
      function(genes) {
        data <- dplyr::select(data, time, o_status, !!!rlang::syms(genes))
        data <- dplyr::mutate(data, group = fun_group(!!!rlang::syms(genes)),
          status = fun_status(o_status))
        fit <- survival::survfit(survival::Surv(time, status) ~ group, data = data)
        diff <- survival::survdiff(survival::Surv(time, status) ~ group, data = data)
        title <- paste0(genes, collapse = " & ")
        p.surv <- survminer::ggsurvplot(fit, data = data,
          pval = T, conf.int = T, risk.table = T, title = title,
          ggtheme = theme_bw()
        )
        p.surv$table <- p.surv$table + labs(x = "time (month)")
        p.surv <- wrap(p.surv, 6.5, 5)
        p.surv <- .set_lab(p.surv, sig(x), "survival curve of", title)
        # timeROC
        if (length(genes) == 1) {
          require(survival)
          cols <- c("blue", "red", "orange")[seq_along(roc_time)]
          roc <- timeROC::timeROC(data$time / 12, data$status, data[[ genes ]], cause = 1,
            times = roc_time)
          legends <- mapply(roc_time, cols, seq_along(roc_time),
            FUN = function(time, col, n) {
              plot(roc, time, col = col, add = (n != 1), title = "")
              glue::glue("AUC at {time} year: {round(roc[[ 'AUC' ]][ paste0('t=', time) ], 2)}")
            })
          legend("bottomright", legends, text.col = cols, bty = "n")
          p.roc <- wrap(recordPlot(), 7, 7)
          p.roc <- .set_lab(p.roc, sig(x), "time ROC")
          data <- dplyr::mutate(data,
            var = ifelse(time / 12 <= 1, "year 1",
              ifelse(time / 12 <= 3, "year 3", "year 5")),
            group = as.factor(o_status), value = !!rlang::sym(genes)
          )
          p.boxplot <- try(wrap(.map_boxplot2(data, T, ylab = genes), 6, 6), T)
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
        namel(p.surv, fit, p.roc, roc, p.boxplot, diff)
      })
    plots <- sapply(c("p.surv", "p.roc", "p.boxplot"), simplify = F,
      function(name) {
        lapply(lst, function(x) x[[ name ]])
      })
    plots <- .set_lab(plots, sig(x), c("Survival", "ROC", "Boxplot"), "plots")
    t.SurvivalPValue <- tibble::tibble(
      name = names(lst), pvalue = vapply(lst, function(x) x$diff$pvalue, double(1))
    )
    t.SignificantSurvivalPValue <- dplyr::filter(t.SurvivalPValue, pvalue < .05)
    x$models <- lapply(lst, function(x) list(fit = x$fit, roc = x$roc, diff = x$diff))
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
      x@plots[[1]] <- plots
    }
    x <- tablesAdd(x, t.SurvivalPValue, t.SignificantSurvivalPValue)
    meth(x)$step1 <- glue::glue("以 R 包 `survival` ({packageVersion('survival')}) 生存分析，以 R 包 `survminer` ({packageVersion('survminer')}) 绘制生存曲线。以 R 包 `timeROC` ({packageVersion('timeROC')}) 绘制 1, 3, 5 年生存曲线。")
    return(x)
  })

sur_group_median <- function(x) {
  ifelse(x > median(x), "High", "Low")
}

sur_status <- function(x) {
  ifelse(x == "Alive", 0, 1)
}
