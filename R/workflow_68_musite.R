# ==========================================================================
# workflow of musite
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_musite <- setClass("job_musite", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "musite",
    info = c("https://github.com/duolinwang/MusiteDeep_web"),
    cite = "[@MusitedeepADWang2020]",
    method = "Python tool `MusiteDeep` was used for protein post-translational modification site prediction and visualization"
    ))

job_musite <- function()
{
  .job_musite()
}

setMethod("step0", signature = c(x = "job_musite"),
  function(x){
    step_message("Prepare your data with function `job_musite`.")
  })

setMethod("step1", signature = c(x = "job_musite"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
