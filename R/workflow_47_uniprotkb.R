# ==========================================================================
# workflow of uniprotkb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_uniprotkb <- setClass("job_uniprotkb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.uniprot.org/help/api_queries"),
    cite = "",
    method = paste0("API of `UniProtKB` (<https://www.uniprot.org/help/api_queries>) ",
      "used for mapping of names or IDs of proteins")
    ))

job_uniprotkb <- function()
{
  .job_uniprotkb()
}

setMethod("step0", signature = c(x = "job_uniprotkb"),
  function(x){
    step_message("Prepare your data with function `job_uniprotkb`.")
  })

setMethod("step1", signature = c(x = "job_uniprotkb"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
