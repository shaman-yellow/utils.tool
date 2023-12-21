# ==========================================================================
# workflow of classyfire
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_classyfire <- setClass("job_classyfire", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@PubchemSubstanKimS2015; @ClassyfireAutDjoumb2016]",
    method = paste0("Database `PubChem` used for querying information (e.g., InChIKey, CID) of chemical compounds; ",
      "Tools of `Classyfire` used for get systematic classification of chemical compounds")
    ))

setGeneric("asjob_classyfire", 
  function(x, ...) standardGeneric("asjob_classyfire"))

setMethod("asjob_classyfire", signature = c(x = "job_herb"),
  function(x){
  })

setMethod("step0", signature = c(x = "job_classyfire"),
  function(x){
    step_message("Prepare your data with function `asjob_classyfire`.")
  })

setMethod("step1", signature = c(x = "job_classyfire"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
