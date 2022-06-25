find_latest_script <- 
  function(
           path = "~/outline",
           pattern = ".*R$"
           ){
    set <- list.files(path = path,
                     pattern = pattern,
                     full.names = T, recursive = T)
    df <- file.info(set) %>%
      dplyr::mutate(file = rownames(.)) %>% 
      dplyr::relocate(file) %>%
      dplyr::arrange(desc(mtime)) %>% 
      dplyr::as_tibble()
    return(df)
  }
ssource <- 
  function(
           path = "~/outline",
           pattern = ".*R$"
           ){
    script <- find_latest_script(path, pattern)$file[1]
    cat("[INFO]: Latest script is:", script, "\n")
    source(script)
  }
rrender <- 
  function(
           path = "~/outline",
           pattern = "md$",
           format = "pdf_document"
           ){
    script <- find_latest_script(path, pattern)$file[1]
    cat("[INFO]: Latest script is:", script, "\n")
    rmarkdown::render(script, format)
  }
msource <- 
  function(
           path = "~/outline",
           pattern = ".*R$",
           block = "symbol",
           symbol = "^## ========== Run block ==========",
           extra = NA,
           source = T
           ){
    script <- find_latest_script(path, pattern)$file[1]
    ## ------------------------------------- 
    df <- data.table::fread(file = script, header = F, sep = NULL) %>% 
      dplyr::rename(code = 1) %>% 
      dplyr::mutate(valid = ifelse(!grepl(paste0("^## -{1,}"), code), T, F),
                    valid = ifelse(grepl(symbol, code), T, valid),
                    row = rownames(.),
                    line = ifelse(valid, row, ","))
    ## ------------------------------------- 
    cluster <- paste(df$line, collapse = ",") %>% 
      gsub(",{2,}", "-", .) %>% 
      strsplit(split = "-") %>% 
      unlist() %>% 
      strsplit(split = ",") %>% 
      lapply(as.numeric)
    ## ------------------------------------- 
    run.block <- grep(symbol, df$code)
    if(length(class) == 0){
      run.block <- length(df$code)
    }
    if(!is.na(extra)){
      run.block <- c(run.block, extra)
    }
    ## ------------------------------------- 
    script <- lapply(cluster,
                     function(vec){
                       if(T %in% (run.block %in% vec))
                         return(vec)
                     }) %>% 
      unlist() %>% 
      df$code[.] %>% 
      paste(collapse = "\n") %>% 
      parse(text = .)
    ## ------------------------------------- 
    if(source){
      source(exprs = script)
    }else{
      return(cluster)
    }
    ## ------------------------------------- 
    # source(exprs = parse(text = script))
  }
