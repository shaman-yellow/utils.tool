# ==========================================================================
# workflow of xcms
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_xcms <- setClass("job_xcms", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@XcmsProcessinSmith2006]",
    method = "R package `xcms` used for MASS data preprocessing (Feature Detection)"
    ))

job_xcms <- function()
{
  .job_xcms()
}

setMethod("step0", signature = c(x = "job_xcms"),
  function(x){
    step_message("Prepare your data with function `job_xcms`.")
  })

setMethod("step1", signature = c(x = "job_xcms"),
  function(x){
    step_message("LC-MS run.")
    return(x)
  })

setGeneric("asjob_xcms", 
  function(x, ...) standardGeneric("asjob_xcms"))

setMethod("asjob_xcms", signature = c(x = "job_msconvert"),
  function(x, snthresh = 5, noise = 50000, peakwidth = c(3, 30),
    ppm = 20, minFraction = 0.1)
  {
    mcm <- e(MCnebula2::new_mcmass(x$metadata, snthresh = snthresh,
        noise = noise, peakwidth = peakwidth, ppm = ppm,
        minFraction = minFraction))
    x <- .job_xcms(object = mcm)
    return(x)
  })
