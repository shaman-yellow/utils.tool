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

setMethod("asjob_survival", signature = c(x = "job_limma"),
  function(x, use.filter, use = "gene_name",
    time = "days_to_last_follow_up", status = "vital_status")
  {
    if (x@step < 1)
      stop("x@step < 1")
    if (x@step > 1) {
      object(x) <- x@params$normed_data
    }
    i.pos <- object(x)$genes[[ use ]] %in% use.filter
    j.pos <- !is.na(object(x)$targets[[ status ]])
    object(x) <- e(limma::`[.EList`(object(x), i.pos, j.pos))
    object <- object(x)
    counts <- object$E
    rownames(counts) <- object$genes[[ use ]]
    data <- as_tibble(t(counts))
    data <- dplyr::bind_cols(data, select(as_tibble(object$targets), -1))
    data <- relocate(data, 1, !!rlang::sym(time), !!rlang::sym(status))
    x <- .job_survival(object = data)
    x@params$time <- time
    x@params$status <- status
    x@params$use.filter <- use.filter
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
  function(x, fun_group = sur_group_median, fun_status = sur_status){
    step_message("Survival test.
      "
    )
    data <- rename(object(x), time = !!rlang::sym(x@params$time),
      status = !!rlang::sym(x@params$status),
    )
    sym <- rlang::sym
    lst <- e(sapply(x@params$use.filter, simplify = F,
      function(gene) {
        data <- select(data, time, status, !!sym(gene))
        data <- mutate(data, group = fun_group(!!sym(gene)), status = fun_status(status))
        fit <- survival::survfit(survival::Surv(time, status) ~ group, data = data)
        p.surv <- survminer::ggsurvplot(fit, data = data,
          pval = T, conf.int = T, risk.table = T, title = gene,
          ggtheme = theme_bw()
        )
        p.surv <- wrap(p.surv, 7, 8)
        namel(p.surv, fit)
      }))
    x@plots[[ 1 ]] <- lapply(lst, function(x) x$p.surv)
    x@params$fit <- lapply(lst, function(x) x$fit)
    return(x)
  })

sur_group_median <- function(x) {
  ifelse(x > median(x), "High", "Low")
}

sur_status <- function(x) {
  ifelse(x == "Dead", 0, 1)
}
