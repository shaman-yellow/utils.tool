# ==========================================================================
# workflow of gsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gsea <- setClass("job_gsea", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://www.gsea-msigdb.org/gsea/downloads.jsp")
    ))

job_gsea <- function()
{
  .job_gsea()
}

setMethod("step0", signature = c(x = "job_gsea"),
  function(x){
    step_message("Prepare your data with function `job_gsea`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_gsea"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })
