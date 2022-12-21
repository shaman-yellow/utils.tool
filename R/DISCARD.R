
## to deal with the bug in `grid.grep`
# valide_vp <- function(ch, chs, sort = F){
  # check <- grid.grep(ch, viewports = T, strict = T, no.match = F)
  # if (is(check, "logical")) {
  #   return(F)
  # }
  # if (sort) {
  #   chs <- sort_vpPaths(chs)
  # }
  # nums <- vapply(chs, FUN.VALUE = 0,
  #                function(ch) {length(grepRaw("::", ch, all = T))})
  # parent.n <- min(nums)
  # match <- grepRaw("::", ch, all = T)
  # end <- match[ parent.n + 1 ] - 1
  # if (!is.na(end)) {
  #   ch.parent <- stringr::str_sub(ch, 1L, end)
  # } else {
  #   ch.parent <- ch
  # }
  # chs <- chs[which( nums > parent.n )]
  # if (length(chs) == 0)
  #   return(F)
  # if (any(!grepl(ch.parent, chs))) {
  #   return(F)
  # }
  # return(ch)
# }

# files <- c("mcn_serum6501.rdata", "serum.tar.gz")
#   urls <- paste0("https://raw.githubusercontent.com/Cao-lab-zcmu/utils_tool/master/",
#                 "inst/extdata/", files)
#   lapply(1:length(files),
#          function(n) {
#            url <- urls[n]
#            data <- RCurl::getURLContent(url)
#            file <- paste0(tmp, "/", files[n])
#            target <- file(file, "wb")
#            writeBin(data, target)
#            close(target)
#          })

# check_pkg <-
#   function(
#            packages = c("data.table", "dplyr", "pbapply", "RCurl", "XML")
#            ){
#     lapply(packages,
#          function(pkg){
#            if (!requireNamespace(pkg, quietly = T))
#              install.packages(pkg)
#          })
#     message("Job Done")
# }

## ---------------------------------------------------------------------- 
# add_spanner <- 
#   function(
#            t,
#            names,
#            group,
#            col_spanner
#            ){
#     envir <- environment()
#     mapply(base_add_spanner, group, col_spanner,
#            MoreArgs = list(envir = envir))
#     return(t)
#   }
# base_add_spanner <- 
#   function(
#            the_group,
#            spanner,
#            envir,
#            names = get("names", envir = envir),
#            t = get("t", envir = envir)
#            ){
#     columns <- names %>%
#       .[grepl(the_group, .)]
#     t <- t %>%
#       tab_spanner(label = spanner,
#                   columns = columns)
#     assign("t", t, envir = envir)
#   }
## ---------------------------------------------------------------------- 

