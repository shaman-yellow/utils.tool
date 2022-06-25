qi_get_format <- 
  function(
           file,
           metadata = F
           ){
    df <- data.table::fread(file) %>% 
      dplyr::as_tibble()
    ## ------------------------------------- 
    meta.ori <- dplyr::slice(df, 1:3)
    ## group
    meta.group <- meta.ori[2, ] %>% 
      unlist() %>% 
      .[which(. != "")] %>% 
      .[1:(length(.) / 2 + 1)] %>% 
      data.table::data.table(group = ., col = names(.)) %>% 
      ## from col to col_end is the sample of group
      dplyr::mutate(col_end = c(col[2:length(col)], NA),
                    col = stringr::str_extract(col, "[0-9]{1,}"),
                    col_end = stringr::str_extract(col_end, "[0-9]{1,}"),
                    col = as.numeric(col),
                    col_end = as.numeric(col_end) - 1) %>% 
      dplyr::slice(1:(nrow(.) - 1))
    ## ------------------------------------- 
    if(metadata){
      meta.sample <- mapply(function(from, to){
                              name <- unlist(meta.ori[3, ], use.names = F)
                              data.table::data.table(sample = name[from:to])
                    }, meta.group$col, meta.group$col_end,
                    SIMPLIFY = F)
      names(meta.sample) <- meta.group$group
      meta.sample <- data.table::rbindlist(meta.sample, idcol = T) %>% 
        dplyr::rename(group = .id)
      return(dplyr::as_tibble(meta.sample))
    }
    ## ---------------------------------------------------------------------- 
    df <- fread(file, skip = 2) %>% 
      dplyr::select(grep("^Compound$|m/z|Retention time", colnames(.)),
                    meta.group$col[1]:meta.group$col_end[nrow(meta.group)])
    return(dplyr::as_tibble(df))
  }
## ------------------------------------- 
qi_as_metabo_inte.table <- 
  function(
           df,
           metadata,
           select
           ){
    select.sam <- dplyr::filter(metadata, group %in% all_of(select))
    ## ---------------------------------------------------------------------- 
    anno <- metadata$group
    names(anno) <- metadata$sample
    ## ------------------------------------- 
    mz_rt <- dplyr::select(df, 2:3) %>% 
      dplyr::rename(mz = 1, rt = 2) %>% 
      dplyr::mutate(Sample = paste0(mz, "__", rt)) %>% 
      dplyr::select(Sample)
    ## ------------------------------------- 
    df.data <- df[, 4:ncol(df)] %>% 
      dplyr::summarise_all(as.character)
    ## ---------------------------------------------------------------------- 
    ## bind id
    qi_format.id <- dplyr::bind_rows(c(Sample = "Lable"), mz_rt)
    ## bind data
    qi_format.df <- dplyr::bind_rows(anno, df.data) %>% 
      dplyr::select(colnames(.)[colnames(.) %in% select.sam$sample])
    ## ------------------------------------- 
    qi_format <- dplyr::bind_cols(qi_format.id, qi_format.df)
    return(qi_format)
  }
