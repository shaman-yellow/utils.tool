# ==========================================================================
# workflow of cardinal
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# http://maldi-msi.org/index.php?option=com_content&view=article&id=189&Itemid=69#ImzMLConverter

.job_cardinal <- setClass("job_cardinal", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "cardinal",
    info = c("https://cardinalmsi.org/"),
    cite = "[@CardinalV3ABemis2023]",
    method = "The R package `Cardinal` used for analyzing mass spectrometry imaging datasets",
    tag = "cardinal",
    analysis = "Cardinal 空间代谢组数据分析"
    ))

job_cardinal <- function()
{
  .job_cardinal()
}

setMethod("step0", signature = c(x = "job_cardinal"),
  function(x){
    step_message("Prepare your data with function `job_cardinal`.")
  })

setMethod("step1", signature = c(x = "job_cardinal"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })


