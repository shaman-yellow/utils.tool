# ==========================================================================
# Following a certain algorithm to get a unique value from the candidate items.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pick_synonym <- 
  function(inchikey2d, inchikey.rdata, synonym.rdata, iupac.rdata) {
    syno <- extract_rdata_list(synonym.rdata)
    syno <- data.frame(data.table::rbindlist(syno))
    if (!is.null(inchikey2d)) {
      inchikey <- extract_rdata_list(inchikey.rdata, inchikey2d)
      meta <- sapply(inchikey, simplify = F, function(x) as.character(x$CID))
      syno$cid <- as.character(syno$cid)
      syno <- group_switch(syno, meta, by = "cid")
    } else {
      syno <- split(syno, ~ cid)
    }
    return(syno)
  }
