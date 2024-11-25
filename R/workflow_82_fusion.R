# ==========================================================================
# workflow of fusion
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_fusion <- setClass("job_fusion", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "fusion",
    info = c("Tutorial: http://gusevlab.org/projects/fusion/"),
    cite = "[@IntegrativeAppGusev2016]",
    method = "",
    tag = "fusion",
    analysis = "FUSION 用于整合 GWAS 与功能基因集数据"
    ))

job_fusion <- function()
{
  .job_fusion()
}

setMethod("step0", signature = c(x = "job_fusion"),
  function(x){
    step_message("Prepare your data with function `job_fusion`.")
  })

setMethod("step1", signature = c(x = "job_fusion"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
