# ==========================================================================
# query classification for compounds using classyfire API
# run after query_inchikey
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
query_classification <- 
  function(
           inchikey2d,
           dir,
           inchikey.rdata = paste0(dir, "/inchikey.rdata"),
           rdata.name = "classification.rdata",
           classyfire_cl = NULL,
           gather_as_rdata = T,
           ...
           ){
    rdata <- paste0(dir, "/", rdata.name)
    classes <- extract_rdata_list(rdata)
    if (!is.null(classes))
      inchikey2d <- inchikey2d[!inchikey2d %in% names(classes)]
    if(length(inchikey2d) == 0)
      return(paste0(dir, "/", rdata.name))
    inchikey_set <- extract_rdata_list(inchikey.rdata, inchikey2d)
    if (is.null(inchikey_set))
      stop("is.null(inchikey_set) == T. File `inchikey.rdata` may not exists.")
    sets <- lapply(inchikey_set, function(df){
                     if("InChIKey" %in% colnames(df))
                       return(df)
           })
    sets <- data.table::rbindlist(sets)
    sets <- dplyr::mutate(sets, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]{1,}"))
    l <- classyfire_get_classification(sets, dir, classyfire_cl = classyfire_cl, ...)
    if (is.logical(l))
      return(paste0(dir, "/", rdata.name))
    if (gather_as_rdata) {
      cat("## gather data\n")
      packing_as_rdata_list(dir, pattern = "^[A-Z]{14}$",
                            rdata = rdata.name, extra = classes)
    }
    return(paste0(dir, "/", rdata.name))
  }

classyfire_get_classification <-
  function(
           sets,
           dir,
           classyfire_cl = NULL,
           log_file = paste0(dir, "/classyfire.log"),
           ...
           ){
    if (file.exists(log_file)){
      log_df <- data.table::fread(log_file)
      sets <- dplyr::filter(sets, !InChIKey %in% log_df$log)
      if(nrow(sets) == 0)
        return(F)
    }
    sets <- split(data.frame(sets), ~ inchikey2d)
    log <- pbapply::pblapply(names(sets), cl = classyfire_cl,
                             function(inchikey2d) {
                               set <- sets[[ inchikey2d ]]
                               unlist(lapply(set[["InChIKey"]], .get_classification,
                                             file = paste0(dir, "/", inchikey2d)),
                                      use.names = F)
                             })
    log <- unlist(log, use.names = F)
    log <- data.frame(log = log)
    if (exists("log_df"))
      log <- dplyr::bind_rows(log_df, log)
    write_tsv(log, file = log_file)
  }

.get_classification <-
  function(inchikey, file){
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
