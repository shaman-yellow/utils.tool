gett <- 
  function(
           obj
           ){
    system(paste("echo", obj, "| xsel -b -i"))
  }
