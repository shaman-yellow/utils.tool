# ==========================================================================
# query inchikey for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_inchikey
#'
#' @title Query InChIkey of compounds via 'InChIkey 2D'
#'
#' @description 
#' The API:
#' url_start = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/", type, "/")
#' url_end = paste0("/property/", paste(get, collapse = ","), "/CSV")
#' url = paste0(url_start, "/", inchikey2d, "/", url_end)
#' 
#' @family queries
#'
#' @name query_inchikey
NULL
#> NULL

#' @export query_inchikey
#' @aliases query_inchikey
#' @description \code{query_inchikey}: ...
#' @rdname query_inchikey
query_inchikey <- 
  function(
           inchikey2d,
           dir,
           rdata.name = "inchikey.rdata",
           curl_cl = NULL,
           gather_as_rdata = T,
           ...
           ){
    rdata <- paste0(dir, "/", rdata.name)
    inchikey_set <- extract_rdata_list(rdata)
    if (!is.null(inchikey_set))
      inchikey2d <- inchikey2d[!inchikey2d %in% names(inchikey_set)]
    if(length(inchikey2d) == 0)
      return(paste0(dir, "/", rdata.name))
    pbapply::pblapply(inchikey2d, pubchem_get_inchikey,
                      dir = dir, cl = curl_cl, ...)
    if (gather_as_rdata) {
      cat("## gather data\n")
      packing_as_rdata_list(dir, pattern = "^[A-Z]{14}$",
                            rdata = rdata.name, extra = inchikey_set)
    }
    return(paste0(dir, "/", rdata.name))
  }

#' @export pubchem_get_inchikey
#' @aliases pubchem_get_inchikey
#' @description \code{pubchem_get_inchikey}: ...
#' @rdname query_inchikey
pubchem_get_inchikey <- 
  function(
           inchikey2d,
           dir,
           type = "inchikey",
           get = "InChIkey",
           ...
           ){
    file <- paste0(dir, "/", inchikey2d)
    if(file.exists(file)){
      csv <- read_tsv(file)
      if("CID" %in% colnames(csv))
        return()
    }
    url_start = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/", type, "/")
    url_end = paste0("/property/", paste(get, collapse = ","), "/CSV")
    url = paste0(url_start, "/", inchikey2d, "/", url_end)
    check <- 0
    while(check == 0 | inherits(check, "try-error")){
      check <- try(csv <- RCurl::getURL(url), silent = T)
    }
    if(grepl("Status: 404", csv)){
      write_tsv(csv, file = file)
      return()
    }
    while(grepl("Status:	503", csv)){
      csv <- RCurl::getURL(url)
    }
    csv <- data.table::fread(text = csv)
    write_tsv(csv, file = file)
  }
