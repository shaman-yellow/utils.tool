try_do <- 
  function(
           text,
           envir
           ){
    check <- 0
    n <- 0
    while(check == 0){
      n <- n + 1
      check <- try(res <- eval(parse(text = text), envir = envir), silent = T)
      if(class(check)[1] == "try-error"){
        check <- 0
      }else{
        check <- 1
      }
      cat("##", "Try...", n, "\n")
    }
    return(res)
  }

