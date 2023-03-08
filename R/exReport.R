# ==========================================================================
# extra tool for report system
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# gather_sections <- function(prefix = "s", envir = parent.frame(),
#   sort = T, get = T)
# {
#   objs <- ls(envir = envir)
#   sections <- objs[ grepl(paste0("^", prefix, "[0-9]"), objs) ]
#   if (sort) {
#     num <- stringr::str_extract(
#       sections, paste0("(?<=", prefix, ")[0-9.]{0,}[0-9]")
#     )
#     sections <- sections[order(as.numeric(num))]
#   }
#   if (get) {
#     sections <- sapply(sections, get, envir = envir, simplify = F)
#   }
#   sections
# }

#' @export .workflow_name
.workflow_name <-
  substitute(
    c("Abstract" = 1, "Introduction" = 1, "Set-up" = 1,
      "Integrate data and Create Nebulae" = 1,
      "Initialize analysis" = 2,
      "Filter candidates" = 2,
      "Filter chemical classes" = 2,
      "Create Nebulae" = 2,
      "Visualize Nebulae" = 2,
      "Nebulae for Downstream analysis" = 1,
      "Statistic analysis" = 2,
      "Set tracer in Child-Nebulae" = 2,
      "Quantification in Child-Nebulae" = 2,
      "Annotate Nebulae" = 2,
      "Query compounds" = 2,
      "Pathway enrichment" = 2,
      "Session infomation" = 1
      ))

