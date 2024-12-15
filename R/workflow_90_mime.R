# ==========================================================================
# workflow of mime
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mime <- setClass("job_mime", 
  contains = c("job"),
  prototype = prototype(
    pg = "mime",
    info = c("https://github.com/l-magnificence/Mime"),
    cite = "[@MimeAFlexiblLiuH2024]",
    method = "",
    tag = "mime",
    analysis = "Mime 机器学习构建模型"
    ))

job_mime <- function()
{
  .job_mime()
}

setMethod("step0", signature = c(x = "job_mime"),
  function(x){
    step_message("Prepare your data with function `job_mime`.")
  })

setMethod("step1", signature = c(x = "job_mime"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_mime"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
