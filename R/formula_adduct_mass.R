formula_adduct_mass <- 
  function(
           formula = NA,
           compound_weight = NA,
           get_formula_weight = F,
           iontype = "neg",
           db_adduct = "[M+H]+,[M+K]+,[M+Na]+,[M+H-H2O]+,[M+H-H4O2]+,[M+NH4]+,[M-H]-,[M+Cl]-,[M-H2O-H]-,[M+Br]-,[M+FA-H]-,[M+ACN-H]-"
           ){
    ## ---------------------------------------------------------------------- 
    welement <- c(H = 1.007825,
                  C = 12.0,
                  N = 14.003074,
                  O = 15.994915,
                  F = 18.998403,
                  P = 30.973762,
                  S = 31.972071,
                  Cl = 34.968853,
                  Br = 78.918336,
                  Na = 22.989770,
                  K = 38.963708)
    ## ---------------------------------------------------------------------- 
    db_adduct <- db_adduct %>%
      strsplit(split = ",") %>% 
      unlist()
    ## ------------------------------------- 
    if(iontype == "neg"){
      db_adduct <- db_adduct %>% 
        .[grepl("(?<=\\])-$", ., perl = T)]
    }else{
      db_adduct <- db_adduct %>% 
        .[grepl("(?<=\\])\\+$", ., perl = T)]
    }
    ## ------------------------------------- 
    db_adduct <- db_adduct %>% 
      gsub("FA", "CO2H2", .) %>% 
      gsub("ACN", "C2H3N", .)
    cat(paste0("[INFO] use adduct: ", paste(db_adduct, collapse = " | "), " \n\n"))
    ## ------------------------------------- 
    ## calculate adduct mass
    adduct_mass <- lapply(db_adduct, get_adduct_mass) %>% 
      unlist()
    ## ---------------------------------------------------------------------- 
    if(is.na(compound_weight)){
      compound_weight <- lapply(formula, element_extract) %>% 
        lapply(element_calculate, welement = welement) %>% 
        unlist()
    }
    if(get_formula_weight){
      return(compound_weight)
    }
    list <- lapply(compound_weight, function(x, plus){x + plus}, adduct_mass) %>% 
      lapply(function(mass, adduct){data.table::data.table(adduct = adduct, mass = mass)},
             adduct = db_adduct)
    names(list) <- formula
    return(list)
  }
get_adduct_mass <- 
  function(
           adduct
           ){
    ## ------------------------------------- 
    welement <- c(H = 1.007825,
                  C = 12.0,
                  N = 14.003074,
                  O = 15.994915,
                  F = 18.998403,
                  P = 30.973762,
                  S = 31.972071,
                  Cl = 34.968853,
                  Br = 78.918336,
                  K = 38.963708,
                  Na = 22.989770)
    ## ------------------------------------- 
    com <- stringr::str_extract_all(adduct, "(?<=[\\-\\+])[A-Z&a-z&0-9]{1,}(?=[\\-\\+\\]])")
    com <- unlist(com)
    ufunc <- stringr::str_extract_all(adduct, "(?<=[0-9|a-z|A-Z])\\+|-(?=[0-9|a-z|A-Z])")
    ufunc <- unlist(ufunc)
    mass <- lapply(com, element_extract)
    mass <- lapply(mass, element_calculate, welement = welement)
    mass <- unlist(mass)
    df <- data.table::data.table(ufunc = ufunc, mass = mass)
    df <- dplyr::mutate(df, mass = ifelse(ufunc == "+", mass, mass * (-1)))
    sum <- sum(df$mass)
    return(sum)
  }
element_calculate <- 
  function(
           df,
           welement
           ){
    weight <- mapply(function(element, number, welement){
                       welement[[element]] * number},
                     df$element, df$number,
                     MoreArgs = list(welement = welement))
    weight <- sum(weight)
    return(weight)
  }
element_extract <- 
  function(
           formula
           ){
    element <- unlist(stringr::str_extract_all(formula, "[A-Z]{1}[a-z]{0,1}"))
    number <- unlist(lapply(paste0("(?<=", element, ")[0-9]{0,}[0-9]{0,}[0-9]{0,}"), mutate_extract,
                     db = formula))
    df <- data.table::data.table(element = element, number = number)
    ## if is NA, set as 1
    df <- dplyr::mutate(df, number = ifelse(number == "", 1, as.numeric(number)))
    if(T %in% duplicated(df$element))
      error
    return(df)
  }
mutate_extract <- 
  function(
           pattern,
           db
           ){
    ch <- stringr::str_extract(db, pattern)
    return(ch)
  }
get_adduct_df <- 
  function(
           adduct
           ){
    com <- stringr::str_extract_all(adduct, "(?<=[\\-\\+])[A-Z&a-z&0-9]{1,}(?=[\\-\\+\\]])")
    com <- unlist(com)
    com <- lapply(com, element_extract)
    ufunc <- stringr::str_extract_all(adduct, "(?<=[0-9|a-z|A-Z])\\+|-(?=[0-9|a-z|A-Z])")
    ufunc <- unlist(ufunc)
    names(com) <- ufunc
    com <- data.table::rbindlist(com, idcol = T)
    return(com)
  }
formula_reshape_with_adduct <- 
  function(
           formula,
           adduct,
           order = F
           ){
    adduct <- get_adduct_df(adduct)
    if(nrow(adduct) == 0)
      return(formula)
    df <- element_extract(formula)
    meta <- environment()
    df <- merge(df, adduct, by = "element", all = T)
    df <- dplyr::mutate(df, number = ifelse(is.na(.id), number.x,
                                            ifelse(.id == "-", number.x - number.y,
                                                   ifelse(is.na(number.x), number.y, number.x + number.y))))
    df <- df[, c("element", "number")]
    df <- summarise_at(group_by(df, element), "number", sum)
    if(order){
      levels <- c("C", "H", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S")
      levels <- levels[levels %in% df$element]
      df <- arrange(df, factor(element, levels = levels))
    }
    df$number <- as.character(df$number)
    ch <- apply(df, 1, paste0)
    ch <- paste(ch, collapse = "")
    return(ch)
    # apply(adduct, 1, base_formula_reshape_with_adduct, envir = meta)
    # return(df)
  }
# base_formula_reshape_with_adduct <- 
  # function(
  #          vector,
  #          envir
  #          ){
  #   df <- get("df", envir = envir)
  #   addn <- as.numeric(vector[3])
  #   if(vector[2] %in% df$element){
  #     num <- df[which(df$element == vector[2]), "number"]
  #     df[which(df$element == vector[2]), "number"] <- ifelse(vector[1] == "+", num + addn, num - addn)
  #   }else{
  #     df <- add_row(df, element = vector[2], number = addn)
  #   }
  #   assign("df", df, envir = envir)
  # }
