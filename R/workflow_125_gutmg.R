# ==========================================================================
# workflow of gutmg
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gutmg <- setClass("job_gutmg", 
  contains = c("job"),
  prototype = prototype(
    pg = "gutmg",
    info = c("https://bio-computing.hrbmu.edu.cn/gutmgene/#/Resource"),
    cite = "",
    method = "",
    tag = "gutmg",
    analysis = "gutMGene 数据库预测"
    ))

job_gutmg <- function(dir = ".")
{
  files <- c(
    "Microbial metabolite-Host Gene.csv",
    "Gut Microbe-Host Gene.csv",
    "Gut Microbe-Microbial metabolite.csv"
  )
  .job_gutmg()
}

setMethod("step0", signature = c(x = "job_gutmg"),
  function(x){
    step_message("Prepare your data with function `job_gutmg`.")
  })

setMethod("step1", signature = c(x = "job_gutmg"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_gutmg"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
