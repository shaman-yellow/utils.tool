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
