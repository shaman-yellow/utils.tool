# ==========================================================================
# workflow of testThis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_testThis <- setClass("job_testThis", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "testThis",
    info = c(""),
    cite = "",
    method = "",
    tag = "",
    analysis = ""
    ))

job_testThis <- function(x)
{
  .job_testThis()
}

setMethod("step0", signature = c(x = "job_testThis"),
  function(x){
    step_message("Prepare your data with function `job_testThis`.")
  })

setGeneric("asjob_testThis", 
  function(x, ...) standardGeneric("asjob_testThis"))

setMethod("asjob_testThis", signature = c(x = "feature"),
  function(x){
    transform(.feature())
  })

setMethod("step1", signature = c(x = "job_testThis"),
  function(x){
    step_message("Quality control (QC).")
    p.ggplotTest <- ggplot(eg, aes(x = x, y = y, color = z)) + geom_point()
    p.ggplotTest <- wrap(p.ggplotTest)
    x <- plotsAdd(x, p.ggplotTest)
    x <- tablesAdd(x, t.eg = eg)
    text <- "__test__"
    transform(.feature())
    message("Run hereing")
    return(x)
  })

setMethod("step2", signature = c(x = "job_testThis"),
  function(x){
    step_message("")
    x$new <- 2:10
    return(x)
  })
