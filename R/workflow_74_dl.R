# ==========================================================================
# workflow of dl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_dl <- setClass("job_dl", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "dl",
    info = c("https://github.com/JinYSun/D-GCAN"),
    cite = "[@PredictionOfDSunJ2022]",
    method = "The deep learning method `D-GCAN` (python) was used for prediction of drug-likeness"
    ))

job_dl <- function()
{
  activate_base("dl")
  sys <- reticulate::import("sys")
  sys$path <- c(sys$path, pg("dl"))
  predict <- reticulate::import("predict")
  x <- .job_dl()
  x$sys <- sys
  x$predict <- predict
  return(x)
}

setMethod("step0", signature = c(x = "job_dl"),
  function(x){
    step_message("Prepare your data with function `job_dl`.")
  })

setMethod("step1", signature = c(x = "job_dl"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
