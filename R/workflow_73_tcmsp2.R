# ==========================================================================
# workflow of tcmsp2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tcmsp2 <- setClass("job_tcmsp2", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "tcmsp2",
    info = c("Tutorial: https://www.tcmsp-e.com/#/home"),
    cite = "[@TcmspADatabaRuJi2014]",
    method = "Website `TCMSP` <https://tcmsp-e.com/> used for data source"
    ))

job_tcmsp2 <- function(herbs)
{
  .job_tcmsp2(params = list(herbs = herbs))
}

setMethod("step0", signature = c(x = "job_tcmsp2"),
  function(x){
    step_message("Prepare your data with function `job_tcmsp2`.")
  })

setMethod("step1", signature = c(x = "job_tcmsp2"),
  function(x, dir = .prefix("tcmsp2", "db")){
    step_message("Download data.")
    dir.create(dir, F)
    items <- c(
      # "disease" = "Related Diseases",
      "ingredients" = "Ingredients",
      "targets" = "Related Targets"
    )
    lapply(x@params$herbs,
      function(herb) {

      })
    return(x)
  })
