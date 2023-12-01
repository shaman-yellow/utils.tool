# ==========================================================================
# workflow of rfsrc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_rfsrc <- setClass("job_rfsrc", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "",
    method = "Package randomForestSRC for feature selection"
    ))

job_rfsrc <- function()
{
  .job_rfsrc()
}

setGeneric("asjob_rfsrc", 
  function(x, ...) standardGeneric("asjob_rfsrc"))

setMethod("asjob_rfsrc", signature = c(x = "job_maf"),
  function(x){
    gene.mut <- e(maftools::genesToBarcodes(object(x), object(x)@gene.summary[[1]], justNames = T))
    dat <- dplyr::select(x$clinical, vital_status)
    rownames(dat) <- x$clinical$Tumor_Sample_Barcode
    dat <- dplyr::mutate(dat, vital_status = as.factor(vital_status))
    # rownames, sample
    # colnames, var, genes
    .job_rfsrc(object = object)
  })

setMethod("step0", signature = c(x = "job_rfsrc"),
  function(x){
    step_message("Prepare your data with function `job_rfsrc`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_rfsrc"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
