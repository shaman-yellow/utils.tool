# ==========================================================================
# workflow of circr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# .job_circr <- setClass("job_circr", 
#   contains = c("job"),
#   representation = representation(
#     object = "ANY",
#     params = "list",
#     plots = "list",
#     tables = "list",
#     others = "ANY"),
#   prototype = prototype(
#     info = c(""),
#     cite = "[@CircrAComputDori2022]",
#     method = "Standalone python script `Circr` used for identifying miRNA-circRNA Associations"
#     ))

# job_circr <- function()
# {
#   .job_circr()
# }

# setMethod("step0", signature = c(x = "job_circr"),
#   function(x){
#     step_message("Prepare your data with function `job_circr`.
#       "
#     )
#   })

# setMethod("step1", signature = c(x = "job_circr"),
#   function(x){
#     step_message("Quality control (QC).")
#     return(x)
  # })
