# ==========================================================================
# workflow of prr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_prr <- setClass("job_prr",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "prr",
    info = c("https://github.com/paulgeeleher/pRRophetic2; https://osf.io/5xvsg; https://github.com/paulgeeleher/pRRophetic2/tree/master/pRRophetic/vignettes; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4167990/"),
    cite = "[@PrropheticAnGeeleh2014]",
    method = "R Package `pRRophetic` was used for Prediction of Clinical Chemotherapeutic Response"
    ))

job_prr <- function()
{
  .job_prr()
}

setMethod("step0", signature = c(x = "job_prr"),
  function(x){
    step_message("Prepare your data with function `job_prr`.")
  })

setMethod("step1", signature = c(x = "job_prr"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
