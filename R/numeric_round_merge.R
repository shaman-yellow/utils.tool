numeric_round_merge <- 
  function(
           main = NULL,
           sub = NULL,
           list = NULL,
           main_col = "mass",
           sub_col = "mz",
           num.tol = 0.002,
           round_factor = 1,
           noise = F
           ){
    if(is.list(list)){
      ## as a list to store data, so that lead to easily apply 'pbapply::pblapply'
      ## for multi-threads running
      main <- list[[1]]
      sub <- list[[2]]
    }
    .Expr <- paste0("as.numeric(", main_col, ")")
    ## to reduce computation, round numeric for limitation
    mutate_main <- mutate(main, .id.m = round(eval(parse(text = .Expr)), round_factor))
    ## ------------------------------------- 
    .Expr <- paste0("as.numeric(", sub_col, ")")
    sub.x <- mutate(sub, .id.m = round(eval(parse(text = .Expr)), round_factor))
    ## escape from missing
    sub.y <- mutate(sub.x, .id.m = .id.m + (1 * 10^-(round_factor)))
    sub <- bind_rows(sub.x, sub.y)
    ## ------------------------------------- 
    ## expand merge
    df <- merge(mutate_main, sub, by = ".id.m", all.x = T, allow.cartesian = T)
    ## ------------------------------------- 
    ## expression
    .Expr <- paste0("abs(as.numeric(", main_col, ") - as.numeric(", sub_col, "))")
    ## eval do expression
    df <- mutate(df, .diff = eval(parse(text = .Expr)))
    ## filter
    df <- filter(df, .diff <= num.tol)
    ## remove the assist col
    df <- select(df, -.id.m, -.diff)
    ## ------------------------------------- 
    if(noise){
      ## as relative intensity
      main <- mutate(main, re.inte = as.numeric(inte) / max(as.numeric(inte)))
      ## eval expression
      .Expr <- paste0("!", main_col, " %in% df$", main_col)
      ## get the noise peak (invalid for sirius)
      df <- filter(main, eval(parse(text = .Expr)))
    }
    ## ------------------------------------- 
    return(df)
  }
