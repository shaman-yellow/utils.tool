# ==========================================================================
# workflow of biomart2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



.job_biomart <- setClass("job_biomart", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/grimbough/biomaRt/issues/61"),
    cite = "[@MappingIdentifDurinc2009]",
    method = "The `biomart` was used for mapping genes between organism (e.g., mgi_symbol to hgnc_symbol)"
    ))

job_biomart2 <- function()
{
  .job_biomart2()
}

setMethod("step0", signature = c(x = "job_biomart2"),
  function(x){
    step_message("Prepare your data with function `job_biomart2`.")
  })

setMethod("step1", signature = c(x = "job_biomart2"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
