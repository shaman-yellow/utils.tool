# ==========================================================================
# workflow of geo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_geo <- setClass("job_geo", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("...")
    ))

job_geo <- function(id)
{
  .job_geo(object = id)
}

setMethod("step0", signature = c(x = "job_geo"),
  function(x){
    step_message("Prepare your data with function `job_geo`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_geo"),
  function(x){
    step_message("Get GEO metadata and information.")
    about <- e(GEOquery::getGEO(object(x)))
    metas <- get_metadata.geo(about)
    prods <- get_prod.geo(metas)
    x@params$about <- about
    x@params$metas <- metas
    x@params$prods <- prods
    return(x)
  })

setMethod("step2", signature = c(x = "job_geo"),
  function(x, filter_regex = NULL){
    step_message("Download geo datasets.")
    e(GEOquery::getGEOSuppFiles(object(x), filter_regex = filter_regex))
    return(x)
  })
