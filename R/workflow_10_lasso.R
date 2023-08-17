# ==========================================================================
# workflow of lasso
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_lasso <- setClass("job_lasso", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/lasso/wiki")
    ))

setGeneric("asjob_lasso", 
  function(x, ...) standardGeneric("asjob_lasso"))

setMethod("asjob_lasso", signature = c(x = "job_limma"),
  function(x, use.filter, use = "gene_name"){
    step_message("The default, 'job_limma' from 'job_tcga' were adapted to
      convertion.
      "
    )
    pos <- object(x)$genes[[use]] %in% use.filter
    object <- e(edgeR::`[.DGEList`(object(x), pos, ))
    .job_lasso(object = object)
  })

setMethod("step0", signature = c(x = "job_lasso"),
  function(x){
    step_message("Prepare your data with function `asjob_lasso`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_lasso"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })

