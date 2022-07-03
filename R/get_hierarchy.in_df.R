get_hierarchy.in_df <- 
  function(
           df,
           col = "Classification"
           ){
    deep <- mutate_get_parent_class(df[[col]], class_cutof = 1) %>% 
      lapply(function(vec){length(vec) + 1}) %>% 
      unlist()
    df <- dplyr::mutate(df, hierarchy = unlist(lapply(eval(parse(text = col)),
                                           function(class){
                                             deep[[class]]
                                           })))
    return(df)
  }
