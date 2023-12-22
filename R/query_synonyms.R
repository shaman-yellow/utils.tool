# ==========================================================================
# query synonyms for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_synonyms
#'
#' @title Query synonyms of compounds via CID
#'
#' @description Bulk search for compound synonyms via pubchem API.
#' (https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/.../synonyms/XML)
#' @family queries
#'
#' @name query_synonyms
NULL
#> NULL

#' @export query_synonyms
#' @aliases query_synonyms
#' @description \code{query_synonyms}: ...
#' @rdname query_synonyms
query_synonyms <- function(cid, dir, rdata.name = "synonyms.rdata", curl_cl = NULL,
  gather_as_rdata = T, group_number = 50, ...)
{
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  rdata <- paste0(dir, "/", rdata.name)
  cid_set <- extract_rdata_list(rdata)
  if (!is.null(cid_set)) {
    cid_set <- data.table::rbindlist(cid_set)
    if (nrow(cid_set) > 0)
      extra <- list(cid_set)
    else
      extra <- NULL
    if("cid" %in% colnames(cid_set)){
      cid_set <- dplyr::distinct(cid_set, cid, syno)
    }
    cid <- cid[!cid %in% cid_set$cid]
    if(length(cid) == 0) 
      return(paste0(dir, "/", rdata.name))
  } else {
    extra <- NULL
  }
  group <- grouping_vec2list(cid, group_number = group_number)
  pbapply::pblapply(group, pubchem_get_synonyms,
    dir = dir, ..., cl = curl_cl)
  if (gather_as_rdata) {
    cat("## gather data\n")
    packing_as_rdata_list(dir, pattern = "^G[0-9]{1,}$",
      dedup = F,
      rdata = rdata.name,
      extra = extra)
  }
  return(paste0(dir, "/", rdata.name))
}

#' @export pubchem_get_synonyms
#' @aliases pubchem_get_synonyms
#' @description \code{pubchem_get_synonyms}: ...
#' @rdname query_synonyms
pubchem_get_synonyms <- function(cid, dir, ...)
{
  savename <- attr(cid, "name")
  file <- paste0(dir, "/", savename)
  cid <- paste(cid, collapse = ",")
  url_start <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  url_end <- "/synonyms/XML"
  url <- paste0(url_start, cid, url_end)
  check <- 0
  while(check == 0 | inherits(check, "try-error")){
    check <- try(text <- RCurl::getURL(url), silent = T)
  }
  while(grepl("Status:	503", text)){
    text <- RCurl::getURL(url)
  }
  text <- XML::xmlToList(text)
  text <- text[names(text) == "Information"]
  text <-
    lapply(text,
      function(list){
        syno <- list[names(list) == "Synonym"]
        syno <-
          lapply(syno,
            function(char){
              if(is.null(char)){
                return(NA)
              }else{
                return(char)
              }
            })
        data.table::data.table(cid = list$CID, syno = unlist(syno))
      })
  text <- data.table::rbindlist(text, fill = T)
  write_tsv(text, filename = file)
}

#' @export extract_rdata_list
#' @aliases extract_rdata_list
#' @description \code{extract_rdata_list}: extract results from .rdata
#' @rdname query_synonyms
extract_rdata_list <- function(rdata, names = NULL){
  path <- get_path(rdata)
  if (!dir.exists(path)) {
    message("The parent directory of `rdata` not exists, so make it herein.")
    dir.create(path)
  }
  if(!file.exists(rdata))
    return()
  load(rdata)
  if (!is.null(names)){
    list <- list[names(list) %in% names]
  }
  return(list)
}

#' @export packing_as_rdata_list
#' @aliases packing_as_rdata_list
#' @description \code{packing_as_rdata_list}: gather table as .rdata
#' @rdname query_synonyms
packing_as_rdata_list <- function(path, pattern, rdata,
  extra = NULL, rm_files = T, dedup = T)
{
  file_set <- list.files(path, pattern = pattern)
  if (!length(file_set))
    return()
  list <- pbapply::pblapply(paste0(path, "/", file_set),
    function(file) {
      res <- try(ftibble(file), T)
      if (inherits(res, "try-error")) {
        data.frame()
      } else {
        res
      }
    })
  names(list) <- file_set
  list <- c(extra, list)
  if (dedup){
    df <- data.table::data.table(name = names(list), n = 1:length(list))
    df <- dplyr::distinct(df, name, .keep_all = T)
    list <- list[df$n]
  }
  if (rm_files){
    lapply(paste0(path, "/", file_set), file.remove)
  }
  save(list, file = paste0(path, "/", rdata))
}

