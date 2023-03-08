# ==========================================================================
# query other property for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases query_iupac
#'
#' @title Query IUPAC name of compounds via 'InChIkey 2D'
#'
#' @description Similar to [query_inchikey()], but get 'IUPACName'.
#' @family queries
#'
#' @name query_iupac
NULL
#> NULL

#' @export query_iupac
#' @aliases query_iupac
#' @description \code{query_iupac}: ...
#' @rdname query_iupac
query_iupac <-
  function(inchikey2d,
           dir,
           rdata.name = "iupac.rdata",
           curl_cl = NULL,
           gather_as_rdata = T,
           ...
           ) {
    query_inchikey(inchikey2d, dir, rdata.name, curl_cl, gather_as_rdata,
                   get = "IUPACName")
}


