# ==========================================================================
# workflow of qiime2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_qiime2 <- setClass("job_qiime2", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("Tutorial: ",
      "https://docs.qiime2.org/2023.7/tutorials/moving-pictures-usage/",
      "https://docs.qiime2.org/2023.7/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq",
      "https://docs.qiime2.org/2023.7/data-resources/",
      collapse = "\n"
    )
    ))

job_qiime2 <- function()
{
  .job_qiime2()
}

setMethod("step0", signature = c(x = "job_qiime2"),
  function(x){
    step_message("Prepare your data with function `job_qiime2`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_qiime2"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })

