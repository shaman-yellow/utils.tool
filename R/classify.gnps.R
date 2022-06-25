classify.gnps <- 
  function(
           inchikey,
           ...
           ){
    ## url path
    url_start <- "https://gnps-classyfire.ucsd.edu/entities/"
    url_end <- ".json"
    ## gather
    url <- paste0(url_start, inchikey, url_end)
    ## ------------------------------------- 
    check <- 0
    while(check == 0 | class(check)[1] == "try-error"){
      check <- try(json <- RCurl::getURL(url), silent = T)
    }
    if(grepl("Key not found", json)){
      return(NA)
    }
    ## read json
    list <- rjson::fromJSON(json_str = json, ...)
    return(list)
  }
