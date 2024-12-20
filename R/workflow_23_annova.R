# ==========================================================================
# workflow of annova
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_annova <- setClass("job_annova", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/annova/wiki"),
    cite = "[@AnnovarFunctiWang2010]",
    method = "`Annova` used for mutation annotation",
    tag = "annova",
    analysis = "Annota 变异注释"
    ))

setGeneric("asjob_annova", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_annova"))

setMethod("asjob_annova", signature = c(x = "job_limma"),
  function(x){
  })

setMethod("step0", signature = c(x = "job_annova"),
  function(x){
    step_message("Prepare your data with function `job_annova`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_annova"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })
