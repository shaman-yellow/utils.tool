# ==========================================================================
# workflow of abismal
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# .job_abismal <- setClass("job_abismal", 
#   contains = c("job"),
#   representation = representation(
#     object = "ANY",
#     params = "list",
#     plots = "list",
#     tables = "list",
#     others = "ANY"),
#   prototype = prototype(
#     pg = "abismal",
#     info = c("https://github.com/smithlabcode/abismal"),
#     cite = "[@FastAndMemoryDeSen2021]",
#     method = "The tool of `abismal` used to map bisulfite-converted reads for methylation analysis"
#     ))

# job_abismal <- function()
# {
#   .job_abismal()
# }

# setMethod("step0", signature = c(x = "job_abismal"),
#   function(x){
#     step_message("Prepare your data with function `job_abismal`.")
#   })

# setMethod("step1", signature = c(x = "job_abismal"),
#   function(x){
#     step_message("Quality control (QC).")
#     return(x)
#   })

# setMethod("set_remote", signature = c(x = "job_abismal"),
#   function(x, wd, postfix = NULL, run_after_cd = NULL, tmpdir = NULL, remote = "remote")
#   {
#     x$wd <- wd
#     x$set_remote <- TRUE
#     x$remote <- "remote"
#     x$map_local <- "abismal_local"
#     return(x)
  # })
