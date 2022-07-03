set_export.no <- 
  function(
           df,
           col = "name"
           ){
    tmp <- data.table::data.table(name = unique(df[[col]])) %>% 
      dplyr::mutate(., No = 1:nrow(.)) %>% 
      dplyr::select(No, name)
    df <- merge(df, tmp, by.x = col, by.y = "name", all.x = T, sort = F) %>% 
      dplyr::relocate(No) %>% 
      dplyr::as_tibble()
    return(df)
  }
