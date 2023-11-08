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
    info = c("Tutorial: https://github.com/satijalab/seurat/wiki"),
    cite = ""
    ))

job_seurat <- function()
{
  .job_seurat()
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    step_message("Prepare your data with function `job_seurat`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_seurat"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_seurat"),
  function(x, wd, postfix, run_after_cd){
    x$postfix
    x$run_after_cd
    x$wd
    return(x)
  })
