get_metadata_from_names <- 
  function(
           names,
           meta.group
           ){
    metadata <- meta.group %>% 
      lapply(function(vec){
               str <- .meta_find_and_sort(names, vec)
           })
    metadata <- mapply(metadata, names(metadata), SIMPLIFY = F,
                       FUN = function(vec, name){
                         df <- data.table::data.table(group = name, sample = vec)
                         return(df)
                       })
    metadata <- data.table::rbindlist(metadata)
    return(metadata)
  }
