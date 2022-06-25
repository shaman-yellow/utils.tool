metabo_collate <- 
  function(
           path = "~/Desktop"
           ){
    ## ----------------------------------------------------------------------
    ## read file
    compound <- list.files(path, pattern = "compound_all.{0,5}.csv$", full.names = T) %>% 
      data.table::fread()
    pathway <- list.files(path, pattern = "pathway_enrichment.{0,5}.csv$", full.names = T) %>% 
      data.table::fread()
    ## ------------------------------------- 
    pathway <- metabo_collate_pathway(pathway)
    ## ------------------------------------- 
    compound <- metabo_collate_compound(compound)
    ## ----------------------------------------------------------------------
    ## gather pathway and compound
    list <- lapply(pathway, merge, y = compound,
                   by.x = "compound", by.y = "Empirical.Compound", all.x = T) %>% 
      lapply(dplyr::as_tibble)
    return(list)
  }
metabo_collate_pathway <- 
  function(
           pathway
           ){
    db <- dplyr::rename(pathway, pathway = V1) %>% 
      by_group_as_list("pathway") %>% 
      lapply(add_row_via_separate_col, only_split = F)
    return(db)
  }
add_row_via_separate_col <- 
  function(
           df_row,
           col = "EC.Hits",
           only_split = F
           ){
    vector <- df_row[[col]] %>% 
      unlist(use.names = F) %>% 
      strsplit(split = ";") %>% 
      unlist(use.names = F)
    ## ------------------ 
    if(only_split == T)
      return(vector)
    ## ------------------ 
    df <- vector %>% 
      data.table::data.table(compound = .)
    ## ------------------
    df_row <- df_row %>%
      .[, which(colnames(.) != col)] %>%
      .[rep(1, nrow(df)), ] %>% 
      dplyr::bind_cols(df) %>%
      dplyr::select(pathway, Hits.sig, Gamma, compound) %>%
      dplyr::as_tibble()
    return(df_row)
  }
## ----------------------------------------------------------------------
metabo_collate_compound <- 
  function(
           compound
           ){
    ## BiGG database download
    if(file.exists("bigg_compound.tsv") == F){
      system("curl http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt > bigg_compound.tsv")
    }
    bigg <- read_tsv("bigg_compound.tsv", fill = T) %>% 
      dplyr::select(universal_bigg_id, name)
    ## ---------------------------------------------------------------------- 
    part_bigg <- compound[["Matched.Compound"]] %>% 
      data.table::data.table(bigg = .) %>% 
      merge(bigg, by.x = "bigg", by.y = "universal_bigg_id") %>% 
      distinct(bigg, .keep_all = T)
    ## ---------------------------------------------------------------------- 
    ## main header: universal_bigg_id name
    ## API: KEGGREST::keggGet(vector) ## return a list
    part_kegg <- compound[["Matched.Compound"]] %>% 
      .[grepl("^C[0-9]{1,100}$", .)] %>% 
      unique() %>%
      data.table::data.table(kegg = .) %>% 
      dplyr::mutate(name = batch_kegg_get(kegg))
    ## ---------------------------------------------------------------------- 
    compound <- compound %>% 
      ## merge bigg compound name
      merge(part_bigg, by.x = "Matched.Compound", by.y = "bigg", all.x = T, sort = F) %>% 
      ## merge kegg compound name
      merge(part_kegg, by.x = "Matched.Compound", by.y = "kegg", all.x = T, sort = F) %>% 
      dplyr::distinct(Query.Mass, Retention.Time, Matched.Compound, .keep_all = T) %>% 
      dplyr::mutate(name = ifelse(is.na(name.y), name.x, name.y)) %>% 
      dplyr::select(-name.x, -name.y) %>% 
      dplyr::as_tibble()
    ## ---------------------------------------------------------------------- 
    return(compound)
  }
## ------------------------------------- 
batch_kegg_get <- 
  function(
           kegg
           ){
    cat("## kegg compound query\n")
    db <- data.table::data.table(kegg = kegg, seq = 1:length(kegg)) %>% 
      dplyr::mutate(index = (seq - seq %% 10) / 10) %>%
      by_group_as_list("index") %>%
      lapply(select, kegg) %>% 
      lapply(unlist) %>% 
      lapply(unname) %>% 
      pbapply::pblapply(.meta_kegg_check_kegg_get) %>%
      unlist() %>% 
      unname()
  }
.meta_kegg_check_kegg_get <- 
  function(
           id_set
           ){
    db <- KEGGREST::keggGet(id_set) %>% 
      lapply(.meta_kegg_check_get_name) %>% 
      unlist() %>% 
      data.table::data.table(kegg = names(.), compound = unname(.))
    df <- data.table::data.table(kegg = id_set) %>% 
      merge(db, by = "kegg", all.x = T, sort = F)
    return(df$compound)
  }
.meta_kegg_check_get_name <- 
  function(
           list,
           id = "ENTRY",
           name = "NAME"
           ){
    id <- list[[id]][1]
    compound <- list[[name]][1]
    if(is.null(compound))
      compound <- NA
    names(compound) <- id
    return(compound)
  }
## ------------------------------------- 
.meta_sort_uniq_df <- 
  function(
           df,
           col,
           pattern_set
           ){
    levels <- df[[col]] %>% 
      .meta_find_and_sort(., pattern_set) %>% 
      unique()
    df[[col]] <- factor(df[[col]], levels = levels)
    df <- df[order(df[[col]]), ] %>% 
      .[!duplicated(.[[col]]), ]
    return(df)
  }
