# ==========================================================================
# workflow of risc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_risc <- setClass("job_risc", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/bioinfoDZ/RISC"),
    cite = "[@RobustIntegratLiuY2021]"
    ))

job_risc <- function(...)
{
  x.lst <- list(...)
  lapply(x.lst,
    function(x) {
      if (!is(x, "job_seurat")) {
        stop("is(x, 'job_seurat') == F")
      }
    })
  x.lst <- e(lapply(x.lst,
      function(x) {
        counts <- object(x)@assays[[ "SCT" ]]
        RISC::readsc()
      }))
  .job_risc()
}

setMethod("step0", signature = c(x = "job_risc"),
  function(x){
    step_message("Prepare your data with function `job_risc`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_risc"),
  function(x){
    step_message("Preprocess before integration.")
    # e(RISC::scFilter())
    # e(RISC::scNormalize())
    e(RISC::scDisperse())

    e(RISC::InPlot())
    return(x)
  })

setMethod("step2", signature = c(x = "job_risc"),
  function(x){
    step_message("Integration and transformed as Seurat object.")
    e(RISC::scMultiIntegrate())
    return(x)
  })
