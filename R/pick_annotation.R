# ==========================================================================
# Following a preset algorithm to get a unique value from the candidate items.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pick_class <- 
  function(inchikey2d, class.rdata, filter = .filter_pick.class, 
           fun = PickClass){
    class <- extract_rdata_list(class.rdata, inchikey2d)
    if (!is.null(filter)) {
      class <- data.frame(data.table::rbindlist(class, idcol = T))
      class <- dplyr::filter(class, !!!filter)
      class <- split(class, ~ .id)
    }
    class <- sapply(inchikey2d, simplify = F,
                    function(key2d) {
                      class[[ key2d ]]$Classification
                    })
    if (!is.null(fun)) {
      class <- lapply(class, fun)
    }
    unlist(class)
  }

.filter_pick.class <- 
  list(quote(!Level %in% dplyr::all_of(c("kingdom", "level 7", "level 8", "level 9"))),
       quote(!grepl("[0-9]|Organ", Classification))
  )

PickClass <- 
  function(class){
    if (is.null(class)) NA
    else tail(class, n = 1)
  }

pick_synonym <- 
  function(inchikey2d = NULL, inchikey.rdata = NULL,
           synonym.rdata, iupac.rdata = NULL,
           filter = .filter_pick.general, fun = PickGeneral) {
    syno <- extract_rdata_list(synonym.rdata)
    syno <- data.frame(data.table::rbindlist(syno))
    if (!is.null(filter)) {
      syno <- dplyr::filter(syno, !!!filter)
    }
    if (!is.null(inchikey2d)) {
      inchikey <- extract_rdata_list(inchikey.rdata, inchikey2d)
      meta <- sapply(inchikey, simplify = F, function(x) as.character(x$CID))
      syno$cid <- as.character(syno$cid)
      syno <- group_switch(syno, meta, by = "cid")
      syno <- sapply(inchikey2d, simplify = F,
                     function(key2d) {
                       if (is.null(syno[[ key2d ]]))
                         return()
                       else
                         syno[[ key2d ]]$syno
                     })
    } else {
      syno <- lapply(split(syno, ~ cid), function(set) set$syno)
    }
    if (!is.null(iupac.rdata)) {
      iupac <- extract_rdata_list(iupac.rdata, inchikey2d)
      if (is.null(inchikey2d)) {
        iupac <- data.table::rbindlist(iupac)
        iupac <- lapply(split(iupac, ~ CID), function(set) set$IUPACName)
      } else {
        iupac <- lapply(iupac, function(set) set$IUPACName)
      }
      syno <- sapply(names(syno), simplify = F,
                     function(name) {
                       c(syno[[ name ]], iupac[[ name ]])
                     })
    }
    if (!is.null(fun)) {
      syno <- lapply(syno, fun)
    }
    unlist(syno)
  }

.filter_pick.general <-
  list(quote(!is.na(syno)),
       quote(!grepl('[0-9]{3}', syno)),
       quote(!grepl('^[A-Z-]{1,5}$', syno)),
       quote(!grepl('^[A-Z0-9]{1,}$', syno)),
       quote(!grepl('(?<=-)[A-Z0-9]{5,}$', syno, perl = T)),
       quote(!grepl('^[0-9-]*$', syno))
  )

PickGeneral <- function(syno,
                        ps = c("^[a-zA-Z]*$", "^[a-zA-Z-]*$",
                               "^[a-zA-Z0-9-]*$", "^[^:]*$")
                        ){
  if (is.null(syno)) return(NA)
  unlist(lapply(ps, function(p) syno[grepl(p, syno)]))[1]
}
