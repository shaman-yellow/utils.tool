order_list <- 
  function(
           list
           ){
    lt <- list()
    length(lt) <- length(list)
    names(lt) <- sort(names(list))
    for(i in names(lt)){
      lt[[i]] <- list[[i]]
    }
    return(lt)
  }
