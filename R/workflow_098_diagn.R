# ==========================================================================
# workflow of diagn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_diagn <- setClass("job_diagn", 
  contains = c("job"),
  prototype = prototype(
    pg = "diagn",
    info = c(""),
    cite = "",
    method = "",
    tag = "",
    analysis = "RNA 数据集验证"
    ))

job_diagn <- function(x, genes = NULL, ...)
{
  valid_job_list(x, "job_limma", 1L)
  projects <- vapply(x, function(x) x$project, character(1))
  jobs <- lapply(x,
    function(x) {
      asjob_diag(x, genes, ...)
    })
  names(jobs) <- projects
  x <- .job_diagn(object = jobs)
  x <- snapAdd(x, "使用数据集 {bind(projects)} 验证。")
  x$projects <- projects
  x$genes <- genes
  return(x)
}

setMethod("step0", signature = c(x = "job_diagn"),
  function(x){
    step_message("Prepare your data with function `job_diagn`.")
  })

setMethod("step1", signature = c(x = "job_diagn"),
  function(x, pattern_control = "control|ctrl|normal",
    target = "group", levels = "guess", ...)
  {
    step_message("Prepare data.")
    object(x) <- lapply(object(x), 
      function(x) {
        step1(
          x, pattern_control = pattern_control, target = target, 
          levels = levels, ...
        )
      })
    return(x)
  })

setMethod("map", signature = c(x = "job_diagn", ref = "job_diag"),
  function(x, ref, ...){
    object(x) <- lapply(object(x), 
      function(x) {
        map(x, ref, ...)
      })
    x$valid_results <- lapply(object(x), function(x) x$valid_results)
    lambdas <- c("lambda.min", "lambda.1se")
    rocs <- sapply(lambdas, simplify = FALSE, 
      function(lam) {
        roc <- lapply(x$valid_results, 
          function(object) {
            object[[ lam ]]$roc
          })
      })
    p.rocs <- lapply(rocs, plot_roc)
    p.rocs <- .set_lab(p.rocs, sig(x), names(p.rocs), "ROC")
    p.rocs <- setLegend(p.rocs, glue::glue("为 {names(p.rocs)} ROC 曲线。"))
    x$p.hps <- lapply(object(x), function(x) x$p.hp_valid)
    x$p.hps <- .set_lab(x$p.hps, sig(x), names(x$p.hps), "feature heatmap in validation dataset")
    x$p.hps <- setLegend(x$p.hps, glue::glue("为特征基因在验证数据集 ({names(x$p.hps)}) 的表达热图。"))
    x$p.rocs <- p.rocs
    x$rocs <- rocs
    return(x)
  })
