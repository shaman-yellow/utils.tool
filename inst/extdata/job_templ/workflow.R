# ==========================================================================
# workflow of seurat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat <- setClass("job_seurat", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/seurat/wiki")
    ))

job_seurat <- function()
{
  .job_seurat()
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    callNextMethod()
    step_message("Prepare your data with function `job_seurat`. ",
      ""
    )
  })

setMethod("step1", signature = c(x = "job_seurat"),
  function(x){
    x <- callNextMethod()
    step_message("Quality control (QC).",
    "This do:",
    ""
    )
    return(x)
  })

