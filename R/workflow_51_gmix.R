# ==========================================================================
# workflow of gmix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gmix <- setClass("job_gmix", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/gmix/wiki"),
    cite = "[@TheDisgenetKnPinero2019; @TheGenecardsSStelze2016; @PharmgkbAWorBarbar2018]",
    method = "Databses of `DisGeNet`, `GeneCards`, `PharmGKB` used for collating disease related targets"
    ))

job_gmix <- function()
{
  .job_gmix()
}

setMethod("step0", signature = c(x = "job_gmix"),
  function(x){
    step_message("Prepare your data with function `job_gmix`.")
  })

setMethod("step1", signature = c(x = "job_gmix"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
