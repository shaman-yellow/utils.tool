# ==========================================================================
# workflow of msconvert
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_msconvert <- setClass("job_msconvert", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"),
    cite = "[@ACrossPlatforChambe2012]",
    method = "The CLI tool of `msconvert` used for converting between file formats of MASS data (from .raw to .mzML)"
    ))

job_msconvert <- function(wd, files)
{
  x <- .job_msconvert(object = files)
  x$wd <- wd
  return(x)
}

setMethod("step0", signature = c(x = "job_msconvert"),
  function(x){
    step_message("Prepare your data with function `job_msconvert`.")
  })

setMethod("step1", signature = c(x = "job_msconvert"),
  function(x, ion = c("positive", "negative"),
    pattern = NULL,
    dir = paste0(ion, "_mzML"),
    script = paste0(.expath, "/msconvert_pre.sh"))
  {
    step_message("Convert file format.")
    ion <- match.arg(ion)
    pbapply::pblapply(object(x),
      function(file) {
        cdRun("sudo bash ", script,
          " ", ".",
          " ", file,
          " ", dir,
          " ", ion,
          path = x$wd)
      })
    x$mzmls <- list.files(paste0(x$wd, "/", dir), ".mzML$", full.names = T)
    try_metadata <- function(files) {
      if (is.null(pattern)) {
        pattern <- c(group = ".")
      }
      data <- group_strings(files, pattern, "file")
      data <- dplyr::mutate(data,
        sample = get_realname(file),
      )
    }
    x$metadata <- try_metadata(x$mzmls)
    x$ion <- ion
    return(x)
  })

.msconvert_help <- function() {
  cdRun("sudo docker run -it --rm chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert --help")
}
