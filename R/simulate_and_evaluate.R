# ==========================================================================
# simulate .mgf from .msp for evaluation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## main function
msp_to_mgf <-
  function(
    name,
    id_prefix,
    path = "~/Downloads/msp/MoNA/",
    write_meta_data = paste0(path, "/", name, ".meta.tsv"),
    fun = "mutate_deal_with_msp_record",
    pre_modify = T,
    ...
    ){
    if(pre_modify == T){
      system(paste0("sed -i 's/\r//g' ", path, "/", name))
    }
    msp <- read_msp(paste0(path, "/", name))
    cache <- new.env()
    store <- new.env()
    assign("id", 0, envir = cache)
    mgf <- paste0(path, "/", name, ".mgf")
    assign("envir_meta", environment(), envir = parent.env(environment()))
    cat("", file = mgf)
    ms_fun <- match.fun(fun)
    pbapply::pblapply(msp[[1]], ms_fun, 
      id_prefix = id_prefix,
      cache = cache,
      store = store,
      ...)
    set <- ls(envir = store)
    meta_data <- lapply(set, get_envir_df,
      envir = store)
    meta_data <- data.table::rbindlist(meta_data, fill = T)
    if(is.null(write_meta_data) == F){
      write_tsv(meta_data, write_meta_data)
    }
    return(meta_data)
  }
read_msp <-
  function(
    filepath
    ){
    msp <- data.table::fread(filepath, sep = NULL, header = F)
  }
get_envir_df <-
  function(
    var,
    envir
    ){
    df <- get(var, envir = envir)
    return(df)
  }

## deal with each line
mutate_deal_with_msp_record <- 
  function(
    ...
    ){
    args <- list(...,
      mass_sep = " ",
      input = c(
        name = "Name",
        mass = "PrecursorMZ", 
        adduct = "Precursor_type",
        formula = "Formula",
        rt = "NA"),
      other = c(
        "Name", "Synon", "DB#", "InChIKey",
        "Precursor_type", "Spectrum_type", "PrecursorMZ",
        "Instrument_type", "Instrument", "Ion_mode",
        "Collision_energy", "Formula",
        "MW", "ExactMass", "Comments")
    )
    do.call(deal_with_msp_record, args)
  }
deal_with_msp_record <-
  function(
    string,
    id_prefix,
    cache,
    store,
    mass_level = 2,
    set_rt = NA,
    mass_sep = "\t",
    id = get("id", envir = cache),
    input = c(name = "NAME",
      mass = "PRECURSORMZ",
      adduct = "PRECURSORTYPE",
      formula = "FORMULA",
      rt = "RETENTIONTIME"),
    other = c("NAME", "PRECURSORMZ", "PRECURSORTYPE",
      "FORMULA", "Ontology", "INCHIKEY", "SMILES",
      "RETENTIONTIME", "CCS", "IONMODE",
      "INSTRUMENTTYPE","INSTRUMENT",
      "COLLISIONENERGY", "Comment", "Num Peaks"),
    output = c(begin = "BEGIN IONS",
      id = "FEATURE_ID=",
      mass = "PEPMASS=",
      charge = "CHARGE=",
      rt = "RTINSECONDS=",
      level = "MSLEVEL=",
      end = "END IONS"),
    add_scans = F
    ){
    ## get name and value
    name = get_name(string)
    name = ifelse(is.na(name) == T, "", name)
    if(grepl("^[A-Z]", name) == T){
      value = get_value(string)
    }
    cat = 0
    if(name == input[["name"]]){
      catapp(output[["begin"]], "\n")
      ## id update
      id = id + 1
      assign("id", id, envir = cache)
      assign("ion", 1, envir = cache)
      ## output
      cat = 1
      p = output[["id"]]
      s = paste0(id_prefix, id)
      ## new var in envir: store
      info <- data.table::data.table(.id = s, name = value)
      assign(paste0(id), info, envir = store)
    }else if(name == input[["mass"]]){
      cat = 1
      p = output[["mass"]]
      s = value
    }else if(name == input[["adduct"]]){
      cat = 1
      p = output[["charge"]]
      s = ifelse(grepl("]-|]+", value) == F, "0",
        ifelse(grepl("]-", value), "1-", "1+"))
      id <- get("id", envir = cache)
      info = get(paste0(id), envir = store)
      info[["charge"]] = s
      assign(paste0(id), info, envir = store)
      assign("adduct", value, envir = cache)
    }else if(name == input[["formula"]]){
      assign("formula", value, envir = cache)
    }else if(name == input[["rt"]]){
      cat = 1
      p = output[["rt"]]
      if(is.na(set_rt)){
        s = value
      }else{
        s = set_rt
      }
    }else if(name == "Num Peaks"){
      cat = 0
      id <- get("id", envir = cache)
      info = get(paste0(id), envir = store)
      if(mass_level == "all"){
        catapp(output[["level"]], "1\n")
        catapp(info[["PRECURSORMZ"]], "100\n", sep = " ")
        ## here, use rcdk to simulate calculate the isotope pattern
        adduct <- get("adduct", envir = cache)
        if(grepl("FA|ACN", adduct)){
          adduct <- gsub("FA", "CO2H2", adduct)
          adduct <- gsub("ACN", "C2H3N", adduct)
        }
        if(adduct != "[M+H-99]+"){
          formula <- get("formula", envir = cache)
          ## according to adduct to revise formula
          formula <- formula_reshape_with_adduct(formula, adduct)
          ## rcdk function
          array <- rcdk::get.isotopes.pattern(rcdk::get.formula(formula))
          apply(array, 1, cat_isotope)
        }
        catapp(output[["end"]], "\n")
        catapp("\n")
        ## begin mass level 2
        catapp(output[["begin"]], "\n")
        catapp(output[["id"]], info[[".id"]], "\n")
        catapp(output[["mass"]], info[["PRECURSORMZ"]], "\n")
        catapp(output[["charge"]], info[["charge"]], "\n")
      }
      if(!is.na(set_rt))
        info[["RETENTIONTIME"]] = set_rt
      catapp(output[["rt"]], info[["RETENTIONTIME"]], "\n")
      catapp(output[["level"]], "2\n")
      catapp("MERGED_STATS=1 / 1 (0 removed due to low quality, 0 removed due to low cosine)\n")
    }else if(grepl("^[0-9]", string)){
      cat = 2
      p = get_name(string, sep = mass_sep)
      s = get_value(string, sep = mass_sep)
    }else if(string == ""){
      ion <- get("ion", envir = cache)
      if(ion == 0){
        return()
      }
      assign("ion", 0, envir = cache)
      cat = 1
      p = output[["end"]]
      s = "\n"
    }
    if(cat == 1){
      catapp(p, s, "\n")
      if(add_scans == T){
        if(p == output[["mass"]]){
          id <- get("id", envir = cache)
          catapp("SCANS=", id, "\n")
        }
      }
    }else if(cat == 2){
      catapp(p, s, "\n", sep = " ")
    }
    ## data store
    if(name %in% other == T){
      id <- get("id", envir = cache)
      info = get(paste0(id), envir = store)
      info[[name]] = value
      assign(paste0(id), info, envir = store)
    }
    return()
    ## output
  }
catapp <-
  function(
    ...,
    sep = "",
    mgf = get("mgf", envir = get("envir_meta"))
    ){
    cat(paste(..., sep = sep), file = mgf, append = T)
  }
cat_isotope <- 
  function(
    vector
    ){
    catapp(vector[1], vector[2] * 100, "\n", sep = " ")
  }
get_value <-
  function(
    string,
    sep = ": "
    ){
    string <- unlist(strsplit(string, split = sep))
    return(string[2])
  }
get_name <-
  function(
    string,
    sep = ": "
    ){
    string <- unlist(strsplit(string, split = sep))
    return(string[1])
  }

# ==========================================================================
# Calculate ion mass of formula with adduct
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

formula_adduct_mass <- 
  function(
    formula = NA,
    compound_weight = NA,
    get_formula_weight = F,
    iontype = "neg",
    db_adduct = "[M+H]+,[M+K]+,[M+Na]+,[M+H-H2O]+,[M+H-H4O2]+,[M+NH4]+,[M-H]-,[M+Cl]-,[M-H2O-H]-,
    [M+Br]-,[M+FA-H]-,[M+ACN-H]-"
    ){
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
    db_adduct <- unlist(strsplit(db_adduct, split = ","))
    if(iontype == "neg"){
      db_adduct <- db_adduct[grepl("(?<=\\])-$", db_adduct, perl = T)]
    }else{
      db_adduct <- db_adduct[grepl("(?<=\\])\\+$", db_adduct, perl = T)]
    }
    db_adduct <- gsub("FA", "CO2H2", db_adduct)
    db_adduct <- gsub("ACN", "C2H3N", db_adduct)
    cat(paste0("[INFO] use adduct: ", paste(db_adduct, collapse = " | "), " \n\n"))
    ## calculate adduct mass
    adduct_mass <- unlist(lapply(db_adduct, get_adduct_mass))
    if(is.na(compound_weight)){
      compound_weight <- lapply(formula, element_extract)
      compound_weight <- unlist(lapply(compound_weight, element_calculate, welement = welement))
    }
    if(get_formula_weight){
      return(compound_weight)
    }
    list <- lapply(compound_weight, function(x, plus){x + plus}, adduct_mass)
    list <- lapply(list,
      function(mass, adduct){
        data.table::data.table(adduct = adduct, mass = mass)
      }, adduct = db_adduct)
    names(list) <- formula
    return(list)
  }
get_adduct_mass <- 
  function(
    adduct
    ){
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
    number <- unlist(lapply(paste0("(?<=", element, ")[0-9]{0,}[0-9]{0,}[0-9]{0,}"),
        mutate_extract, db = formula))
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
    df <- dplyr::summarise_at(dplyr::group_by(df, element), "number", sum)
    if(order){
      levels <- c("C", "H", "Cl", "F", "I", "K", "N", "Na", "O", "P", "S")
      levels <- levels[levels %in% df$element]
      df <- dplyr::arrange(df, factor(element, levels = levels))
    }
    df$number <- as.character(df$number)
    ch <- apply(df, 1, paste0)
    ch <- paste(ch, collapse = "")
    return(ch)
    # apply(adduct, 1, base_formula_reshape_with_adduct, envir = meta)
    # return(df)
  }

# ==========================================================================
# add noise to mgf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

collate_as_noise_pool <- 
  function(
    origin_list,
    valid_list
    ){
    ## filter origin_list
    args <- list(list = origin_list, discard_level1 = T, only_peak_info = T, mass_shift = F)
    ## get mz and intensity
    cat("## Catch main peak information\n")
    origin_list <- do.call(spectrum_add_noise, args)
    ## discard the NULL data
    cat("## Discard empty dataset\n")
    order_list <- origin_list[vapply(origin_list, is.data.frame, logical(1), USE.NAMES = F)]
    ## order the origin_list and valid_list according to .id
    ## first, filter the origin_list, only the .id in valid_list is reserved.
    origin_list <- origin_list[names(origin_list) %in% names(valid_list)]
    ## keep identical
    valid_list <- valid_list[names(valid_list) %in% names(origin_list)]
    ## order
    cat("## Order the lists...\n")
    origin_list <- order_list(origin_list)
    valid_list <- order_list(valid_list)
    origin_list <- lapply(origin_list,
      function(df) {
        df[[ "mass" ]] <- as.double(df[[ "mass" ]])
        df[[ "inte" ]] <- as.double(df[[ "inte" ]])
        df
      })
    cat("## Merge to get noise list\n")
    noise_list <- pbapply::pblapply(names(origin_list),
      function(name) {
        df <- tol_mergeEx(
          origin_list[[ name ]], valid_list[[ name ]],
          main_col = "mass", sub_col = "mz"
        )
        df$rel.int. <- df$inte / max(origin_list[[ name ]]$inte)
        df
      })
    noise_df <- data.table::rbindlist(noise_list, fill = T)
    dplyr::mutate(noise_df, mass = as.numeric(mass), inte = as.numeric(inte))
  }

filter_mgf <- 
  function(
           filter_id = prapare_inst_data(.MCn.structure_set)$.id,
           file
           ){
    mgf <- read_msp(file)
    start <- which(mgf$V1 == "BEGIN IONS")
    end <- which(mgf$V1 == "")
    id <- mgf[grepl("FEATURE_ID", mgf$V1), ]
    id <- stringr::str_extract(id, "(?<==).*$")
    list <- pbapply::pbmapply(
      function(start, end, mgf) {
        dplyr::slice(mgf, start:end)
      }, start, end, MoreArgs = list(mgf = mgf), SIMPLIFY = F
    )
    names(list) <- id
    if(is.null(filter_id) == F){
      list <- list[names(list) %in% filter_id]
    }
    return(list)
  }

## mutate function of tol_merge, get the non-merged data
tol_mergeEx <- 
  function(main,
           sub,
           main_col = "mz",
           sub_col = "mz",
           tol = 0.002,
           bin_size = 1
           ){
    if (main_col == sub_col) {
      new_name <- paste0(sub_col, ".sub")
      colnames(sub)[colnames(sub) == sub_col] <- new_name
      sub_col <- new_name
    }
    main$...seq <- 1:nrow(main)
    backup <- main
    ## to reduce computation, round numeric for limitation
    ## main
    main$...id <- round(main[[ main_col ]], bin_size)
    ## sub
    sub.x <- sub.y <- sub
    sub.x$...id <- round(sub.x[[ sub_col ]], bin_size)
    sub.y$...id <- sub.x$...id + ( 1 * 10^-bin_size )
    sub <- rbind(sub.x, sub.y)
    ## expand merge
    df <- merge(main, sub, by = "...id", all.x = T, allow.cartesian = T)
    df$...diff <- abs(df[[ main_col ]] - df[[ sub_col ]])
    df <- dplyr::filter(df, ...diff <= !!tol)
    ## get the non-merged
    df <- backup[!backup$...seq %in% df$...seq, ]
    df$...seq <- NULL
    df
  }

mass_shift <- 
  function(
           df,
           merge = T,
           sep = " ",
           int.sigma = 1,
           re.ppm = 1e-6,
           global.sigma = 10/3 * re.ppm,
           indivi.sigma = 10/3 * re.ppm,
           sub.factor = 0.03,
           .noise_pool = noise_pool,
           alpha = 0.2,
           ...
           ){
    df <- dplyr::mutate(df, mass = as.numeric(mass), inte = as.numeric(inte))
    ## intensity variation
    var <- rnorm(nrow(df), 1, int.sigma)
    df <- dplyr::mutate(df, inte = inte * var)
    ## subtract according to max intensity
    df <- dplyr::mutate(df, inte = round(inte - max(inte) * sub.factor, 0))
    ## if intensity less than 0, discard
    df <- dplyr::filter(df, inte > 0)
    ## almost one peak, discard the data
    if(nrow(df) <= 1)
      return()
    ## global shift
    var <- rnorm(1, 0, global.sigma)
    df <- dplyr::mutate(df, mass = mass + mass * var)
    ## individual shift
    var <- rnorm(nrow(df), 0, indivi.sigma)
    df <- dplyr::mutate(df, mass = round(mass + mass * var, 4))
    ## add noise peak
    ## random drawn noise peak from noise pool
    noise <- .noise_pool[sample(1:nrow(.noise_pool), round(alpha * nrow(df))), ]
    ## re-size intensity
    noise <- dplyr::mutate(noise, inte = max(df$inte) * rel.int.)
    ## bind into df
    df <- bind_rows(df, dplyr::select(noise, mass, inte))
    if(merge){
      df <- dplyr::mutate(df, V1 = paste0(mass, sep, inte))
      df <- dplyr::select(df, V1)
    }
    return(df)
  }

# ==========================================================================
# spectra add noise
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## for .mgf data
spectrum_add_noise <- 
  function(list, cl = NULL, filter_empty = T, ...)
  {
    list <- pbapply::pblapply(list,
      function(df, discard_level1 = F, mass_process_level_1 = F,
        mass_process_level_2 = T, ...)
      {
        mass_level <- df$V1[grepl("MSLEVEL", df$V1)]
        ## process level 1
        if(mass_level == "MSLEVEL=1"){
          if(discard_level1)
            return()
          if(mass_process_level_1)
            df <- mass_process_level_1(df, ...)
          return(df)
          ## process level 2
        }else{
          if(mass_process_level_2)
            df <- mass_process_level_2(df, ...)
          return(df)
        }
      }, cl = cl, ...)
    ## filter the empty spectrum
    if(filter_empty){
      list <- list[!vapply(list, is.null, logical(1), USE.NAMES = F)]
    }
    return(list)
  }

discard_level1 <- 
  function(list) {
    spectrum_add_noise(
      list = list, mass_process_level_2 = F, discard_level1 = T
    )
  }

mass_process_level_1 <- 
  function(df, ...){}

mass_process_level_2 <- 
  function(df, mass_shift = T, ...){
    list <- separate_peak_info(df, ...)
    if(!mass_shift)
      return(list)
    list[[2]] <- mass_shift(list[[2]], merge = T, ...)
    if(length(list) == 2)
      return()
    df <- rbindlist(list)
    return(df)
  }

separate_peak_info <- 
  function(df, sep = " ", only_peak_info = F, ...)
  {
    peak_row <- grep("^[0-9]", df$V1)
    peak_info <- dplyr::slice(df, peak_row)
    peak_info <- tidyr::separate(peak_info, col = "V1", into = c("mass", "inte"), sep = sep)
    if(only_peak_info == T)
      return(peak_info)
    list <- list(dplyr::slice(df, 1:(min(peak_row) - 1)),
      peak_info,
      dplyr::slice(df, (max(peak_row) + 1):nrow(df))
    )
    return(list)
  }

# ==========================================================================
# other
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mgf_add_anno.gnps <- 
  function(df){
    slice_line <- list("1:3", "4:6", "7:nrow(df)")
    list <- lapply(slice_line, function(lines){
      list <- dplyr::slice(df, eval(parse(text = lines)))
      return(list)
    })
    ## scans
    scans <- stringr::str_extract(list[[1]][2, 1], "[0-9]{1,}$")
    scans <- c(V1 = paste0("SCANS=", scans))
    ## merge
    merge <- c(
      V1 = paste0("MERGED_STATS=1 / 1 ",
        "(0 removed due to low quality, 0 removed due to low cosine)"
      )
    )
    dplyr::bind_rows(list[[1]], scans, list[[2]], merge, list[[3]])
  }

simulate_gnps_quant <- 
  function(
           meta,
           save_path,
           file = paste0(save_path, "/", "quant.csv"),
           rt = 1000,
           area = 10000,
           id = ".id",
           mz = "PRECURSORMZ",
           simu_id = "row ID",
           simu_mz = "row m/z",
           simu_rt = "row retention time",
           simu_quant = "sample.mzML Peak area",
           return_df = F
           ){
    meta <- dplyr::select(meta, all_of(c(id, mz)))
    meta <- dplyr::mutate(meta, rt = rt, sample = area)
    colnames(meta) <- 
      mapply_rename_col(
        .find_and_sort_strings(colnames(meta), c(id, mz, "rt", "sample")),
        c(simu_id, simu_mz, simu_rt, simu_quant),
        colnames(meta)
      )
    if(return_df)
      return(meta)
    write.table(meta, file = file, sep = ",", row.names = F, col.names = T, quote = F)
  }

# ==========================================================================
# evaluate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

stat_classify <- 
  function(id_set, resClass, id2Key, mcn, ref)
  {
    set <- get_parent_classes(resClass, mcn)
    cat("stat:", resClass, "\n")
    list <- pbapply::pblapply(
      id_set,
      function(id) {
        class <- ref[[ id2Key[[ id ]] ]]
        if (is.null(class)) {
          stat <- data.frame(id = id, evaluate = NA)
          return(stat)
        }
        if (resClass %in% class[["Classification"]]){
          stat <- data.frame(id = id, evaluate = "true")
        } else {
          if (class[3,]$Classification %in% set){
            ## at least the cluster is T in "class" level
            evaluate <- "latent"
          } else {
            evaluate <- "false"
          }
          stat <- data.frame(id = id, evaluate = evaluate)
        }
        return(stat)
      }
    )
    df <- data.table::rbindlist(list)
    return(df)
  }

table_app <- 
  function(df, col = "evaluate", prop = T){
    if (!is.data.frame(df))
      return()
    stat <- table(df[[col]])
    if (prop) {
      sum <- sum(stat)
      stat <- prop.table(stat)
    }
    stat <- dplyr::bind_rows(stat)
    stat$sum <- sum
    return(stat)
  }

stat_identification <- 
  function(id_lst, id_key2d, ref)
  {
    id_stat <- lapply(id_lst, merge, y = id_key2d, by = ".features_id", all.x = T)
    id_stat <- lapply(id_stat, merge, y = ref, by = ".features_id", all.x = T)
    id_stat <- lapply(id_stat, dplyr::mutate,
      evaluate = ifelse(inchikey2d.x == inchikey2d.y, "true", "false"))
    lst <- lapply(id_stat, table_app)
    df <- data.table::rbindlist(lst, idcol = T, fill = T)
    dplyr::rename(df, class.name = .id)
  }

# ==========================================================================
# visualize the results of evaluation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

visualize_stat <- 
  function(df, mcn,
           palette = ggsci::pal_npg()(9),
           group_palette = ggsci::pal_rickandmorty()(12),
           ylab = "ratio",
           xlab = "classes",
           fill_lab = "type",
           col_max = 500,
           weight = c(pl = 1, pm = 1.2, pr = 1)
           )
  {
    parent_class <- unlist(
      lapply(get_parent_classes(df$class.name, mcn),
        function(vec) if (length(vec) == 0) NA else tail(vec, n = 1)),
      use.names = F
    )
    df.back <- df
    df <- dplyr::mutate(
      df, parent_class = ifelse(is.na(parent_class), class.name, !!parent_class)
    )
    annotation <- df
    df <- tidyr::gather(df, "type", "value", -class.name, -parent_class, -sum)
    df <- dplyr::mutate(
      df, class.name = stringr::str_wrap(class.name, width = 25),
      parent_class = stringr::str_wrap(parent_class, width = 25),
      type = as.character(type),
      type = Hmisc::capitalize(type)
    )
    pm <- ggplot(data = df, aes(x = class.name, y = value, fill = type)) +
      geom_col(width = 0.7, position = "stack") +
      scale_fill_manual(values = palette) +
      labs(y = Hmisc::capitalize(ylab),
        x = Hmisc::capitalize(xlab),
        fill = Hmisc::capitalize(fill_lab)
        ) +
      coord_flip()
    pm.theme <- pm.theme2 <- theme(
      axis.text.y = element_blank(),
      text = element_text(family = .font, face = "bold"),
      plot.title = element_text(hjust = 0.3)
    )
    pm.theme$legend.position <- "none"
    pm.theme2$legend.position <- "bottom"
    pm.legend <- MCnebula2:::.get_legend(pm + pm.theme2)
    pm <- pm + pm.theme
    pr <- ggplot(df.back) +
      geom_col(aes(x = class.name,
          y = ifelse(sum >= col_max, col_max, sum)),
        width = 0.7, fill = "#709AE1FF", alpha = 0.7) +
      ylim(0, col_max) +
      labs(x = "", y = "Classified number") +
      coord_flip() +
      theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = .font, face = "bold")) +
      geom_blank()
    pl <- ggplot(annotation) +
      geom_col(aes(x = stringr::str_wrap(class.name, width = 25), y = 1,
          fill = stringr::str_wrap(parent_class, width = 25)),
        width = .7) +
      labs(fill = "Super Classes", x = "", y = "Classes") +
      theme_minimal() +
      scale_fill_manual(
        values = colorRampPalette(group_palette)(length(unique(annotation$parent_class)))
        ) +
      theme(
        text = element_text(face = "bold", family = .font),
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_blank(),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "left",
        panel.grid = element_blank()) +
      coord_flip()
    lst <- lapply(namel(pl, pm, pr), as_grob)
    lst <- frame_col(weight, lst)
    frame_row(c(lst = 10, pm.legend = .5), namel(lst, pm.legend))
  }

visualize_statComplex <- 
  function(
    df_list,
    mcn,
    weight = c(pl = 1, pm = 1, pr = 1),
    ylab = "ratio",
    xlab = "classes",
    fill_lab = "type",
    palette = ggsci::pal_npg()(9),
    mutate_palette = c("true" = palette[3], "latent" = palette[2], "false" = palette[1],
      "medium_noise" = "#FED439FF", "high_noise" = "#8A4198FF"),
    extra_palette = c("sum" = "#95CC5EFF",
      "medium_noise" = "#FED439FF", "high_noise" = "#8A4198FF"),
    group_palette = ggsci::pal_rickandmorty()(12),
    y_cut_left = c(50, 450),
    y_cut_right = NULL,
    y_cut_left_breaks = c(50, seq(100, 450, by = 100)),
    y_cut_right_breaks = NULL
  )
  {
    ## get parent class
    df_list.back <- df_list
    df_list <- lapply(df_list,
      function(df){
        parent_class <- unlist(
          lapply(get_parent_classes(df$class.name, mcn),
            function(vec) if (length(vec) == 0) NA else tail(vec, n = 1)),
          use.names = F
        )
        df <- dplyr::mutate(
          df, parent_class = ifelse(is.na(parent_class), class.name, !!parent_class),
          st.true = 0, en.true = true,
          st.latent = en.true, en.latent = st.latent + latent,
          st.false = en.latent, en.false = st.false + false
        )
      })
    ## group draw
    annotation <- df_list[["origin"]]
    pl <- ggplot(annotation) +
      geom_col(aes(x = stringr::str_wrap(class.name, width = 25), y = 1,
          fill = stringr::str_wrap(parent_class, width = 25)),
        width = .7) +
      labs(fill = "Super classes", x = "", y = "Classes") +
      theme_minimal() +
      scale_fill_manual(
        values = colorRampPalette(group_palette)(length(unique(annotation$parent_class)))
        ) +
      theme(
        text = element_text(face = "bold", family = .font),
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_blank(),
        legend.key.height = unit(1.5, "cm"),
        legend.position = "left",
        panel.grid = element_blank()) +
      coord_flip()
    ## initial stat
    mutate_origin <- tidyr::gather(
      df_list[[ "origin" ]], "type", "value", true, false, latent
    )
    mutate_origin <- dplyr::mutate(
      mutate_origin,
      y = as.numeric(apply(mutate_origin, 1, function(v) v[[paste0("st.", v[["type"]])]] )),
      yend = as.numeric(apply(mutate_origin, 1, function(v) v[[paste0("en.", v[["type"]])]] ))
    )
    ## medium noise dirft
    mergeMutate <-
      function(v1, v2, df_list){
        df <- merge(df_list[[v1]], df_list[[v2]], by = "class.name", all.x = T)
        df <- dplyr::mutate(df, flow1 = "true", flow2 = "latent")
        df <- tidyr::gather(df, "type", "value", flow1, flow2)
        ## calculate segment from y to yend
        dplyr::mutate(
          df, y = as.numeric(apply(df, 1, function(v) v[[paste0("en.", v[["value"]], ".x")]] )),
          yend = as.numeric(apply(df, 1, function(v) v[[paste0("en.", v[["value"]], ".y")]] )),
          exclude = ifelse(is.na(yend), T, F), y = ifelse(is.na(yend), 0, y),
          yend = ifelse(is.na(yend), 1, yend)
        )
      }
    medium_noise_df <- mergeMutate("origin", "medium_noise", df_list)
    medium_noise_df <- dplyr::filter(
      medium_noise_df, y != yend, !exclude, class.name %in% mutate_origin$class.name
    )
    ## high noise drift
    high_noise_df <- mergeMutate("medium_noise", "high_noise", df_list)
    high_noise_df <- dplyr::filter(
      high_noise_df, y != yend, !exclude, class.name %in% mutate_origin$class.name
    )
    pm <- ggplot() +
      ## origin
      geom_segment(
        data = mutate_origin,
        aes(x = class.name, xend = class.name, y = y, yend = yend, color = type),
        size = 7) +
      ## medium noise drift
      geom_segment(
        data = medium_noise_df,
        aes(x = class.name, xend = class.name, y = y, yend = yend, color = "medium_noise"),
        size = 7) +
      ## high noise drift
      geom_segment(
        data = high_noise_df,
        aes(x = class.name, xend = class.name, y = y, yend = yend, color = "high_noise"),
        size = 7) +
      ## the point indicate the start of noise drift
      geom_segment(
        data = medium_noise_df,
        aes(x = class.name, xend = class.name, y = ifelse(yend > y, y - 0.001, y + 0.001),
          yend = y, color = "medium_noise"),
        arrow = arrow(length = unit(10, "pt")), size = 0.5, lineend = "round") +
      ## the point indicate the start of high noise drift
      geom_segment(
        data = high_noise_df,
        aes(x = class.name, xend = class.name, y = ifelse(yend > y, y - 0.001, y + 0.001),
          yend = y, color = "high_noise"),
        arrow = arrow(length = unit(10, "pt")), size = 0.5, lineend = "round") +
      scale_color_manual(
        values = mutate_palette,
        labels = c(sum = "Sum",
          true = "True",
          false = "False",
          latent = "Latent",
          medium_noise = "Medium noise",
          high_noise = "High noise")) +
      labs(
        y = Hmisc::capitalize(ylab),
        x = Hmisc::capitalize(xlab),
        color = Hmisc::capitalize(fill_lab)) +
      coord_flip()
    pm.theme <- theme(
      legend.position = "bottom",
      axis.text.y = element_blank(),
      text = element_text(family = .font, face = "bold"),
      plot.title = element_text(hjust = 0.3)
    )
    pm.legend <- MCnebula2:::.get_legend(pm + pm.theme)
    pm.theme$legend.position <- "none"
    pm <- pm + pm.theme
    extra_list <- lapply(df_list.back, dplyr::select, class.name, sum)
    extra.noise_df <- merge(
      extra_list[[ "medium_noise" ]], extra_list[[ "high_noise" ]],
      by = "class.name", all.x = T
    )
    extra.noise_df <- merge(
      extra.noise_df, extra_list[[ "origin" ]], by = "class.name", all.y = T
    )
    pr <- ggplot() +
      ## origin sum
      geom_segment(
        data = extra_list[["origin"]],
        aes(x = class.name, xend = class.name, y = 0, yend = sum, color = "sum"),
        size = 7) +
      ## medium_noise drift
      geom_segment(
        data = dplyr::mutate(extra.noise_df, sum.x = ifelse(is.na(sum.x), 0, sum.x)),
        aes(x = class.name, xend = class.name, y = sum, yend = sum.x, color = "medium_noise"),
        size = 7) +
      ## high_noise drift
      geom_segment(
        data = dplyr::mutate(
          dplyr::filter(extra.noise_df, is.na(sum.x) == F),
          sum.x = ifelse(is.na(sum.y), 0, sum.x),
          sum.y = ifelse(is.na(sum.y), sum.x, sum.y)),
        aes(x = class.name, xend = class.name, y = sum.x, yend = sum.y, color = "high_noise"),
        size = 7) +
      scale_color_manual(
        values = extra_palette,
        labels = c(sum = "Sum", medium_noise = "Medium noise", high_noise = "High noise")) +
      geom_blank()
    pr.lab <- pr.lab1 <- pr.lab2 <- labs(x = NULL, y = NULL, color = "Type")
    pr.lab1$y <- "Classified number"
    pr.lab2$y <- "..."
    pr.theme <- theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      text = element_text(family = .font, face = "bold")
    )
    pr.legend <- MCnebula2:::.get_legend(pr + pr.theme + pr.lab)
    pr.theme$legend.position <- "none"
    pr <- pr + pr.theme
    pr1 <- pr + pr.lab1 +
      coord_flip(ylim = y_cut_left) +
      geom_hline(yintercept = c(50), linetype = "dashed", size = 0.7,
        color = "grey") +
      scale_y_continuous(breaks = y_cut_left_breaks)
    pr1 <- as_grob(pr1)
    if (is.null(y_cut_right)) {
      pr2 <- nullGrob()
      pr.w <- c(pr1 = 2, pr.legend = 1)
    } else {
      pr2 <- pr + pr.lab2 +
        coord_flip(ylim = y_cut_right) +
        scale_y_continuous(breaks = y_cut_right_breaks)
      pr2 <- as_grob(pr2)
      pr.w <- c(pr1 = 2, pr2 = 1, pr.legend = 1)
    }
    pr <- frame_col(pr.w, namel(pr1, pr2, pr.legend))
    ## gather all
    pl <- as_grob(pl)
    pm <- as_grob(pm)
    lst <- frame_col(weight, namel(pl, pm, pr))
    frame_row(c(lst = 10, pm.legend = .5), namel(lst, pm.legend))
  }

visualize_comparison <- 
  function(
    list1,
    list2,
    ylim_min = 50,
    from = c("MCnebula2", "GNPS"),
    ylab = "Classified number",
    xlab = "Classes",
    color_lab = "Methods",
    palette = ggsci::pal_npg()(9)
    )
  {
    ## select in common classification
    common.class <- dplyr::select(merge(list1[[1]], list2[[1]], by = "class.name"), 1)
    ## filter classification
    list <- list(list1, list2)
    list <- lapply(list, function(lst){
      lst <- lapply(lst, dplyr::select, class.name, sum)
      lst <- lapply(lst, merge, y = common.class, by = "class.name", all.y = T)
      lst <- lapply(lst, dplyr::mutate, sum = ifelse(is.na(sum), 0, sum))
      df <- data.table::rbindlist(lst, idcol = T)
      dplyr::rename(df, group = .id)
    })
    list <- mapply(list, from, SIMPLIFY = F,
      FUN = function(df, VALUE) {
        dplyr::mutate(df, from = !!VALUE)
      })
    ## for segment
    df <- merge(list[[1]], list[[2]], by = c("group", "class.name"))
    group_levels = c(origin = "Origin", medium_noise = "Medium noise",
      high_noise = "High noise")
    df$group <- vapply(df$group, function(x) group_levels[[ x ]], character(1))
    ## for point
    df2 <- dplyr::bind_rows(list[[1]], list[[2]])
    df2$group <- vapply(df2$group, function(x) group_levels[[ x ]], character(1))
    ## plot figure
    p <- ggplot() +
      geom_segment(
        data = df,
        aes(x = class.name, xend = class.name, y = sum.x, yend = sum.y),
        color = "black") +
      geom_point(
        data = df2,
        aes(x = class.name, y = sum, color = from),
        size = 3, position = "identity") +
      scale_color_manual(values = palette) +
      scale_y_continuous(breaks = c(ylim_min, 300, 600, 900, 1200)) +
      labs(
        y = Hmisc::capitalize(ylab),
        x = Hmisc::capitalize(xlab),
        color = Hmisc::capitalize(color_lab)) +
      coord_flip(ylim = c(ylim_min, max(df2$sum) * 1.1)) +
      facet_wrap(~ factor(group, levels = group_levels)) +
      theme(
        legend.position = "bottom",
        text = element_text(family = .font, face = "bold"),
        plot.title = element_text(hjust = 0.3)) +
      geom_blank()
    grob <- as_grob(p)
    attr(grob, "data") <- list
    return(grob)
  }

visualize_summary <- 
  function(summary){
    data <- data.frame(do.call(rbind, summary))
    data <- dplyr::mutate_if(data, is.list, unlist)
    data <- dplyr::mutate(
      data, type = form(rownames(data)),
      type = factor(type, levels = type)
    )
    ## draw plot of classified number
    theme1 <- theme2 <- theme(text = element_text(family = .font))
    theme2$legend.position <- "none"
    p.num <- ggplot(dplyr::mutate(data, type = factor(type, levels = rev(type)))) +
      geom_col(aes(x = type, y = sum, fill = type), width = .5, fill = "#709AE1FF") +
      geom_text(aes(x = type, y = sum, label = round(sum)),
        family = .font, hjust = -1) +
      labs(x = "", y = "Sum") + 
      coord_flip(ylim = c(0, 250)) +
      theme_classic() + theme2
    theme3 <- theme1
    theme3$legend.position <- "bottom"
    theme3$panel.spacing <- u(0, line)
    data.ratioFal <- dplyr::mutate(data, `non-false` = 100 - false)
    data.ratioFal <- tidyr::gather(data.ratioFal, evaluate, value, -sum, -type)
    data.ratioSt <- dplyr::mutate(
      data.ratioFal, change = (max(sum) - sum ) / max(sum) * 100,
      evaluate = ifelse(evaluate == "false", "lost", "non-lost"),
      change = ifelse(evaluate == "lost", change, 100 - change),
      change = round(change, 1)
    )
    ## draw plot of stability
    data.ratioSt <- dplyr::filter(data.ratioSt, type != "Origin")
    p.ratioSt <- ggplot(data.ratioSt) +
      geom_col(aes(x = 0, y = change, fill = form(evaluate))) +
      geom_text(data = dplyr::filter(data.ratioSt, evaluate == "lost"),
        aes(x = -2, y = 0, label = paste0(change, "%")),
        hjust = .5, family = .font) +
      coord_polar(theta = "y", direction = -1) +
      xlim(-2, .5) +
      ylim(0, 100) +
      labs(fill = "") +
      scale_fill_manual(values = c("#D9D9D9", "#91D1C2")) +
      facet_wrap(~ type, nrow = 1) +
      theme_void()
    p.ratioSt <- sep_legend(p.ratioSt, theme3)
    lostRate <- mean(dplyr::filter(data.ratioSt, evaluate == "lost")$change) / 100
    data.ratioRelFal <- dplyr::mutate(
      data.ratioFal, rel.value = ifelse(evaluate == "false",
        100 - (100 - value) * (1 - lostRate),
        value * (1 - lostRate)
      )
    )
    p.ratioRelFal <- ggplot(data.ratioRelFal) +
      geom_col(aes(x = 0, y = rel.value, fill = form(evaluate))) +
      geom_text(data = dplyr::filter(data.ratioRelFal, evaluate == "false"),
        aes(x = -2, y = 0, label = paste0(round(rel.value, 1), "%")),
        hjust = .5, family = .font) +
      coord_polar(theta = "y", direction = -1) +
      xlim(-2, .5) +
      ylim(0, 100) +
      labs(fill = "") +
      scale_fill_manual(values = c("#E64B35FF", "#FED439FF")) +
      facet_wrap(~ type, nrow = 1) +
      theme_void()
    p.ratioRelFal <- sep_legend(p.ratioRelFal, theme3)
    namel(p.num, p.ratioSt, p.ratioRelFal)
  }

visualize_idRes <- 
  function(
    list,
    palette = ggsci::pal_simpsons()(9),
    ylab = "Ratio",
    xlab = "Classes",
    color_lab = "type"
  )
  {
    list.name <- names(list)
    list <- lapply(list, tidyr::gather,
      "type", "value", -class.name
    )
    list <- lapply(list, dplyr::mutate,
      class.name = stringr::str_wrap(class.name, width = 25),
      type = as.character(type),
      type = Hmisc::capitalize(type))
    df <- data.table::rbindlist(list, idcol = T)
    df <- dplyr::filter(df, type == "True")
    line_df <- tidyr::spread(df, .id, value)
    p <- ggplot(data = df, aes(x = class.name, y = value, color = .id)) +
      geom_segment(data = line_df,
        aes(x = class.name, xend = class.name,
          y = eval(parse(text = paste0("`", list.name[1], "`"))),
          yend = eval(parse(text = paste0("`", list.name[2], "`")))),
        color = "black") +
      geom_point(size = 3, position = "identity") +
      scale_color_manual(values = palette) +
      labs(y = Hmisc::capitalize(ylab),
        x = Hmisc::capitalize(xlab),
        color = Hmisc::capitalize(color_lab)) +
      coord_flip() +
      theme(legend.position = "bottom",
        text = element_text(family = .font, face = "bold"),
        plot.title = element_text(hjust = 0.3)) + 
      geom_blank()
    as_grob(p)
  }

gtext90 <- function(label, fill, rot = 90) {
  rect_mcnebula <- roundrectGrob(gp = gpar(col = "transparent", fill = fill))
  gtext <- gtext(label, list(cex = 1.2, col = "white"), rot = rot)
  ggather(rect_mcnebula, gtext)
}
