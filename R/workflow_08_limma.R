# ==========================================================================
# workflow of limma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_limma <- setClass("job_limma", 
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

job_limma <- function(DGEList)
{
  if (!is(DGEList, 'DGEList'))
    stop("is(DGEList, 'DGEList') == F")
  .job_limma(object = DGEList)
}

setMethod("step0", signature = c(x = "job_limma"),
  function(x){
    step_message("Prepare your data with function `job_limma`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_limma"),
  function(x, group = x@object$samples$group){
    step_message("Preprocess expression data.")
    x@params$group <- group
    return(x)
  })

