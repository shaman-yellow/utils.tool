# ==========================================================================
# workflow of epifactor
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_epifactor <- setClass("job_epifactor", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "epifactor",
    info = c("https://epifactors.autosome.org/description"),
    cite = "[@Epifactors2022Maraku2023]",
    method = "Database `EpiFactors` used for screening epigenetic regulators"
    ))

job_epifactor <- function(use = c("protein"), version = "v2.0")
{
  x <- .job_epifactor()
  x$use <- match.arg(use)
  x$version <- version
  return(x)
}

setMethod("step0", signature = c(x = "job_epifactor"),
  function(x){
    step_message("Prepare your data with function `job_epifactor`.")
  })

setMethod("step1", signature = c(x = "job_epifactor"),
  function(x){
    step_message("Obtain epigenetic regulators.")
    file <- switch(x$use, protein = "EpiGenes_main.csv")
    url <- paste0("https://epifactors.autosome.org/public_data/", x$version, "/", file)
    data <- ftibble(RCurl::getURL(url))
    data <- .set_lab(data, sig(x), paste0("all ", x$use), "of epigenetic regulators")
    x@tables[[ 1 ]] <- nl(x$use, list(data))
    return(x)
  })
