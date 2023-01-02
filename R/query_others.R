# ==========================================================================
# query other property for compounds using pubchem API
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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


