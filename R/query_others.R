# ==========================================================================
# query other property for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
query_iupuc <-
  function(inchikey2d,
           dir,
           rdata.name = "iupuc.rdata",
           curl_cl = NULL,
           gather_as_rdata = T,
           ...
           ) {
    query_inchikey(inchikey2d, dir, rdata.name, curl_cl, gather_as_rdata,
                   get = "IUPACName")
}


