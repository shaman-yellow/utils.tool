# ==========================================================================
# workflow of seuratSp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seuratSp <- setClass("job_seuratSp", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://satijalab.org/seurat/articles/spatial_vignette.html")
    ))

job_seuratSp <- function()
{
  .job_seuratSp()
}

setMethod("step0", signature = c(x = "job_seuratSp"),
  function(x){
    step_message("Prepare your data with function `job_seuratSp`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_seuratSp"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })
