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
    pg = "seurat",
    info = c("Tutorial: https://github.com/satijalab/seurat/wiki"),
    cite = "",
    method = "",
    tag = "seurat",
    analysis = ""
    ))

job_seurat <- function()
{
  .job_seurat()
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    step_message("Prepare your data with function `job_seurat`.")
  })

setMethod("step1", signature = c(x = "job_seurat"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_seurat"),
  function(x, wd, postfix = NULL, run_after_cd = NULL, tmpdir = NULL, remote = "remote")
  {
    x$wd <- wd
    x$set_remote <- T
    x$remote <- "remote"
    x$map_local <- "seurat_local"
    return(x)
  })
