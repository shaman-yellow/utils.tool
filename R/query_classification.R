# ==========================================================================
# query classification for compounds using classyfire API
# run after query_inchikey
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_classification
#'
#' @title Query classification of compounds via packgage of 'classifyerR'
#'
#' @description The used function is:
#' \code{classyfireR::get_classification(inchikey)}
#' @family queries
#'
#' @name query_classification
NULL
#> NULL

#' @export query_classification
#' @aliases query_classification
#' @description \code{query_classification}: ...
#' @rdname query_classification
query_classification <- function(query, dir, 
  inchikey.rdata = paste0(dir, "/inchikey.rdata"), rdata.name = "classification.rdata",
  classyfire_cl = NULL, gather_as_rdata = TRUE, ...)
{
  rdata <- paste0(dir, "/", rdata.name)
  db_classes <- extract_rdata_list(rdata)
  fun_format.dbInchi <- function(sets) {
    sets <- lapply(sets,
      function(df){
        if("InChIKey" %in% colnames(df))
          return(df)
      })
    sets <- data.table::rbindlist(sets)
    sets <- dplyr::mutate(sets, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]{1,}"))
    sets
  }
  if (all(grpl(query, "^[0-9]+$"))) {
    ## CID to inchikey2d
    inchikey_set <- extract_rdata_list(inchikey.rdata, as.character(query))
    inchikey_set <- fun_format.dbInchi(inchikey_set)
    query <- inchikey_set$inchikey2d[ match(as.integer(query), as.integer(inchikey_set$CID)) ]
  }
  if (!is.null(db_classes))
    query <- query[!query %in% names(db_classes)]
  if (!length(query))
    return(paste0(dir, "/", rdata.name))
  if (!exists("inchikey_set")) {
    inchikey_set <- extract_rdata_list(inchikey.rdata, query)
    sets <- fun_format.dbInchi(inchikey_set)
  } else {
    sets <- dplyr::filter(inchikey_set, inchikey2d %in% !!query)
  }
  l <- classyfire_get_classification(sets, dir, classyfire_cl = classyfire_cl, ...)
  if (is.logical(l))
    return(paste0(dir, "/", rdata.name))
  if (gather_as_rdata) {
    cat("## gather data\n")
    packing_as_rdata_list(dir, pattern = "^[A-Z]{14}$",
      rdata = rdata.name, extra = db_classes)
  }
  return(paste0(dir, "/", rdata.name))
}

#' @export classyfire_get_classification
#' @aliases classyfire_get_classification
#' @description \code{classyfire_get_classification}: ...
#' @rdname query_classification
classyfire_get_classification <- function( sets, dir, classyfire_cl = NULL,
  log_file = paste0(dir, "/classyfire.log"), ...)
{
  if (file.exists(log_file)){
    log_df <- data.table::fread(log_file)
    sets <- dplyr::filter(sets, !InChIKey %in% log_df$log)
    if(nrow(sets) == 0)
      return(FALSE)
  }
  sets <- split(data.frame(sets), ~ inchikey2d)
  log <- pbapply::pblapply(names(sets), cl = classyfire_cl,
    function(inchikey2d) {
      set <- sets[[ inchikey2d ]]
      unlist(lapply(set[["InChIKey"]], .get_classification,
          file = paste0(dir, "/", inchikey2d)),
        use.names = FALSE)
    })
  log <- unlist(log, use.names = FALSE)
  log <- data.frame(log = log)
  if (exists("log_df"))
    log <- dplyr::bind_rows(log_df, log)
  write_tsv(log, file = log_file)
}

.get_classification <- function(inchikey, file)
{
  if(!file.exists(file)){
    ch <- classyfireR::get_classification(inchikey)
  }else{
    return()
  }
  if(is.null(ch)){
    return(inchikey)
  }else{
    ch <- classyfireR::classification(ch)
    write_tsv(ch, file)
  }
}
