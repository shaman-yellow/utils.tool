# ==========================================================================
# workflow of testThis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_testThis <- setClass("job_testThis", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "testThis",
    info = c(""),
    cite = "",
    method = "",
    tag = "",
    analysis = "测试分析"
    ))

job_testThis <- function(x)
{
  .job_testThis()
}

setMethod("step0", signature = c(x = "job_testThis"),
  function(x){
    step_message("Prepare your data with function `job_testThis`.")
  })

setGeneric("asjob_testThis", 
  function(x, ...) standardGeneric("asjob_testThis"))

setMethod("asjob_testThis", signature = c(x = "feature"),
  function(x){
    feature <- resolve_feature_snapAdd_onExit("x", x)
    snapAdd_onExit("x", "这是之后的内容。")
    x <- .job_testThis()
    return(x)
  })

setMethod("step1", signature = c(x = "job_testThis"),
  function(x, a, b, c, d, ...){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("step2", signature = c(x = "job_testThis"),
  function(x){
    step_message("")
    x$new <- 2:10
    return(x)
  })
