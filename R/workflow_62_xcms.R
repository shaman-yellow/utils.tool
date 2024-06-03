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
    method = "R package `xcms` used for MASS data preprocessing (Feature Detection)",
    tag = "mass"
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
    if (nrow(x$metadata) < 2) {
      # remove the steps which needs Theplural samples
      object(x)@detectFlow %<>% MCnebula2::delete_layers(3)
    }
    object(x) <- MCnebula2::run_lcms(object(x))
    return(x)
  })

setMethod("step2", signature = c(x = "job_xcms"),
  function(x, savepath = timeName(x$ion))
  {
    step_message("Output the mgf file.")
    dir.create(savepath, F)
    x$mgf <- paste0(savepath, "/msms.mgf")
    object(x) <- e(MCnebula2::run_export(object(x), saveMgf = x$mgf))
    x$mgf <- format_mgf.xcms(x$mgf, x$ion)
    x$quant <- e(MCnebula2::features_quantification(object(x)))
    return(x)
  })

setGeneric("asjob_xcms", 
  function(x, ...) standardGeneric("asjob_xcms"))

setMethod("asjob_xcms", signature = c(x = "job_msconvert"),
  function(x, snthresh = 5, noise = 50000, peakwidth = c(3, 30),
    ppm = 20, minFraction = 0.1)
  {
    metadata <- x$metadata
    ion <- x$ion
    mcm <- e(MCnebula2::new_mcmass(metadata, snthresh = snthresh,
        noise = noise, peakwidth = peakwidth, ppm = ppm,
        minFraction = minFraction))
    x <- .job_xcms(object = mcm)
    x$metadata <- metadata
    x$ion <- ion
    return(x)
  })

format_mgf.xcms <- function(mgf, ion, rm.title = T, revise.charge = T) {
  lines <- readLines(mgf)
  if (rm.title) {
    lines <- lines[ !grpl(lines, "TITLE") ]
  }
  if (revise.charge) {
    ion <- match.arg(ion, c("negative", "positive"))
    lines[ grpl(lines, "^CHARGE") ] <- paste0("CHARGE=1", switch(ion, negative = "-", positive = "+"))
  }
  new.mgf <- gs(mgf, "\\.mgf$", "Format.mgf")
  writeLines(lines, new.mgf)
  return(new.mgf)
}
