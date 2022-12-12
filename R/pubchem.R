check_pkg <-
  function(
           packages = c("data.table", "dplyr", "pbapply", "RCurl", "XML")
           ){
    lapply(packages,
         function(pkg){
           if (!requireNamespace(pkg, quietly = T))
             install.packages(pkg)
         })
    message("Job Done")
}
pubchem_get_synonyms <- 
  function(
           cid,
           dir,
           curl_cl = NULL,
           gather_as_rdata = T,
           group_number = 50,
           ...
           ){
    rdata <- paste0(dir, "/", "cid.rdata")
    ## extract as list
    cid_set <- extract_rdata_list(rdata)
    ## as data.table
    cid_set <- data.table::rbindlist(cid_set)
    ## !duplicated
    if("cid" %in% colnames(cid_set)){
      cid_set <- cid_set %>% 
        dplyr::distinct(cid, syno)
    }
    ## -------------------------------------
    ## exclude existing
    cid <- cid[!cid %in% cid_set$cid]
    ## ------------------------------------- 
    if(length(cid) == 0)
      return()
    ## ------------------------------------- 
    group <- grouping_vec2list(cid, group_number = group_number)
    group <<- group
    ## ------------------------------------- 
    pbapply::pblapply(group, base_pubchem_get_synonyms,
                      dir = dir, cl = curl_cl, ...)
    ## ------------------------------------- 
    cat("## gather synonyms\n")
    packing_as_rdata_list(dir, pattern = "^G[0-9]{1,}$",
                          dedup = F,
                          rdata = "cid.rdata",
                          ## data.table as list
                          extra = list(cid_set))
  }
base_pubchem_get_synonyms <- 
  function(
           cid,
           dir,
           ...
           ){
    savename <- attr(cid, "name")
    file <- paste0(dir, "/", savename)
    ## gather cid and sep by ,
    cid <- paste(cid, collapse = ",")
    ## use cid
    url_start <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    ## get as XML
    url_end <- "/synonyms/XML"
    ## paste as url
    url <- paste0(url_start, cid, url_end)
    ## ------------------------------------- 
    check <- 0
    while(check == 0 | class(check)[1] == "try-error"){
      check <- try(text <- RCurl::getURL(url), silent = T)
    }
    ## ------------------------------------- 
    while(grepl("Status:	503", text)){
      text <- RCurl::getURL(url)
    }
    ## "PUGREST.BadRequest"
    ## ------------------------------------- 
    ## convert to list
    text <- XML::xmlToList(text)
    ## only 'information'
    text <- text[names(text) == "Information"]
    ## in list to separate data
    text <- lapply(text, function(list){
                     syno <- list[names(list) == "Synonym"]
                     syno <- lapply(syno,
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
    ## ------------------------------------- 
    ## save data
    write_tsv(text, filename = file)
  }
## ------------------------------------- 
grouping_vec2list <- 
  function(
           vector,
           group_number,
           byrow = F
           ){
    if(length(vector) < group_number){
      attr(vector, "name") <- "G1"
      return(list(vector))
    }
    ## if grouped, the rest number
    rest <- length(vector) %% group_number
    ## assign group
    group <- matrix(vector[1:(length(vector) - rest)],
                    ncol = group_number,
                    byrow = byrow) %>% 
      ## use apply to multiple list
      apply(1, c, simplify = F) %>% 
      ## gather the rest vector
      c(., list(tail(vector, n = rest))) %>% 
      ## add group name
      mapply(FUN = function(vec, name){
               attr(vec, "name") <- name
               return(vec)
           }, ., paste0("G", 1:length(.)),
           SIMPLIFY = F)
    ## ------------------------------------- 
    if(rest == 0)
      group <- group[1:(length(group) - 1)]
    return(group)
  }
## ---------------------------------------------------------------------- 
extract_rdata_list <- 
  function(
           rdata,
           names = NA
           ){
    if(!file.exists(rdata))
      return()
    load(rdata)
    if(!is.na(names[1])){
      list <- list[names(list) %in% names]
    }
    return(list)
  }
packing_as_rdata_list <- 
  function(
           path,
           pattern,
           rdata,
           extra = NULL,
           rm_files = T,
           dedup = T
           ){
    file_set <- list.files(path, pattern = pattern)
    if(length(file_set) == 0)
      return()
    ## read as list
    list <- pbapply::pblapply(paste0(path, "/", file_set), read_tsv)
    names(list) <- file_set
    ## merge
    list <- c(extra, list)
    ## according to name, unique
    if(dedup){
      df <- data.table::data.table(name = names(list), n = 1:length(list))
      df <- dplyr::distinct(df, name, .keep_all = T)
      list <- list[df$n]
    }
    ## rm origin file sets
    if(rm_files){
      lapply(paste0(path, "/", file_set), file.remove)
    }
    ## save as rdata
    save(list, file = paste0(path, "/", rdata))
  }
write_tsv <-
  function(x, filename, col.names = T, row.names = F){
    write.table(x, file = filename, sep = "\t",
                col.names = col.names, row.names = row.names, quote = F)
  }
read_tsv <- function(path){
  file <- data.table::fread(input = path, sep = "\t",
                            header = T, quote = "", check.names = F)
  return(file)
}
