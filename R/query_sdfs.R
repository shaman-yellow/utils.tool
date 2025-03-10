# ==========================================================================
# query sdfs for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_sdfs
#'
#' @title Query SDF files of compounds via CID
#'
#' @description Bulk search for compound sdfs via pubchem API.
#' <https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/SDF>
#' @family queries
#'
#' @name query_sdfs
NULL
#> NULL

#' @export query_sdfs
#' @aliases query_sdfs
#' @description \code{query_sdfs}: ...
#' @rdname query_sdfs
query_sdfs <- function(cid, dir = "SDFs", curl_cl = NULL, 
  group_number = 50, filename = "all_compounds.sdf", ...)
{
  dir.create(dir, FALSE)
  group <- grouping_vec2list(cid, group_number = group_number, TRUE)
  files <- pbapply::pbvapply(group, pubchem_get_sdf,
    dir = dir, ..., cl = curl_cl, FUN.VALUE = character(1))
  lines <- lapply(files, readLines)
  lines <- lapply(lines, function(ls) c(ls, ""))
  lines <- unlist(lines)
  writeLines(lines, file <- file.path(dir, filename))
  file.remove(files)
  file
}

#' @export pubchem_get_sdf
#' @aliases pubchem_get_sdf
#' @description \code{pubchem_get_sdf}: ...
#' @rdname query_sdfs
pubchem_get_sdf <- function(cid, dir, ...)
{
  savename <- attr(cid, "name")
  file <- paste0(dir, "/", savename, ".sdf")
  cid <- paste(cid, collapse = ",")
  url_start <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  url_end <- "/SDF"
  url <- paste0(url_start, cid, url_end)
  check <- 0
  while(check == 0 | inherits(check, "try-error")){
    check <- try(text <- RCurl::getURL(url), silent = TRUE)
  }
  while(grepl("Status:	503", text)){
    text <- RCurl::getURL(url)
  }
  writeLines(text, file)
  return(file)
}
