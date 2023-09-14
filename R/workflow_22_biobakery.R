# ==========================================================================
# workflow of biobakery
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_biobakery <- setClass("job_biobakery", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/biobakery/biobakery/wiki/biobakery_workflows#2-metagenome-profiling")
    ))

job_biobakery <- function()
{
  .job_biobakery()
}

setMethod("step0", signature = c(x = "job_biobakery"),
  function(x){
    step_message("Prepare your data with function `job_biobakery`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_biobakery"),
  function(x){
    step_message("Quality control (QC).
      "
    )
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_biobakery"),
  function(x, wd, postfix, db = "/data/hlc/biobakery_workflows_databases"){
    x$postfix
    x$run_after_cd <- paste0("export KNEADDATA_DB_HUMAN_GENOME=", db)
    x$wd
    return(x)
  })
