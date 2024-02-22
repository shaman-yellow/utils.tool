# ==========================================================================
# workflow of metabo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_metabo <- setClass("job_metabo", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://www.metaboanalyst.ca"),
    cite = "[@Metaboanalyst4Chong2018]",
    method = "`MetaboAnalyst` used for metabolomic data analysis"
    ))

job_metabo <- function()
{
  .job_metabo()
}

setMethod("step0", signature = c(x = "job_metabo"),
  function(x){
    step_message("Prepare your data with function `job_metabo`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_metabo"),
  function(x, cpds, libname = "kegg_pathway"){
    step_message("Enrichment analysis.")
    require(MetaboAnalystR)
    mSet <- e(MetaboAnalystR::InitDataObjects("conc", "msetora", FALSE))
    mSet <- e(MetaboAnalystR::Setup.MapData(mSet, cpds))
    mSet <- e(MetaboAnalystR::CrossReferencing(mSet, "name"))
    mSet <- e(MetaboAnalystR::CreateMappingResultTable(mSet))
    mSet <- e(MetaboAnalystR::SetMetabolomeFilter(mSet, F))
    mSet <- e(MetaboAnalystR::SetCurrentMsetLib(mSet, libname, 2))
    if (T) {
      mSet <- e(MetaboAnalystR::CalculateHyperScore(mSet))
      mtfig <- .mtb_file("metabolites_ORA_dot_", libname)
      mSet <- e(MetaboAnalystR::PlotEnrichDotPlot(mSet, "ora", mtfig$pn, "pdf", width = NA))
      ## plots
      x@plots[[ 1 ]] <- c(mtfig$file_fig)
    } else {
      mtfig <- list(name = "Name")
    }
    ## tables
    mapped <- tibble::as_tibble(data.frame(mSet$dataSet$map.table))
    mapped <- dplyr::filter(mapped, KEGG != "NA")
    if (!is.null(mSet$analSet$ora.mat)) {
      table <- as_tibble(mSet$analSet$ora.mat)
    } else {
      table <- c()
    }
    x@tables[[ 1 ]] <- c(n.l(mtfig$name, list(table)), namel(mapped))
    ## params
    hits <- n.l(mtfig$name, mSet$analSet$ora.hits)
    x$hits <- hits
    x$mSet <- mSet
    return(x)
  })

setMethod("step2", signature = c(x = "job_metabo"),
  function(x, query){
    step_message("Get KEGG infomation to Genes.")
    mapped <- filter(x@tables$step1$mapped, Query %in% dplyr::all_of(query))
    x$db_cpd <- e(KEGGREST::keggGet(mapped$KEGG), text = "compounds")
    x$db_module <- e(lapply(x$db_cpd,
      function(db) {
        KEGGREST::keggGet(names(db$MODULE))
      }), text = "modules")
    extract <- stringr::str_extract_all
    db_genes.kegg <- e(lapply(x$db_module,
        function(db) {
          lapply(db,
            function(sub) {
              qr <- extract(names(sub$ORTHOLOGY), "K[0-9]+")
              qr <- unique(unlist(qr))
              if (F) {
                res <- lapply(grouping_vec2list(qr, 10, T),
                  function(x) KEGGREST::keggGet(x))
                res <- unlist(res, recursive = F)
                lapply(res, function(x) {
                  .plist(data = x$symbol)
                  })
              } else {
                qr
              }
            })
        }), text = "genes")
    x$db_genes.kegg <- unlist(db_genes.kegg, use.names = F)
    # extract <- stringr::str_extract
    # pblapply <- pbapply::pblapply
    # db_genes <- e(pblapply(x$db_genes.kegg,
    #   function(obj) {
    #     extract(obj@data)
    # }))
    return(x)
  })

.plist <- setClass("plist", 
  contains = c(),
  representation = representation(data = "ANY"),
  prototype = NULL)

.mtb_file <- function(..., postfix = "dpi72.pdf") {
  name <- paste0(unlist(list(...)), collapse = "")
  pn <- paste0(get_savedir("figs"), "/", name, "_")
  file <- paste0(pn, postfix)
  file_fig <- file_fig(name, file)
  namel(pn, file, name, file_fig)
}

meta_metabo_pathway <- function(
  export = NA, mz_rt = NA, p_col = NA, extra_entity = NA,
  only_return = F,
  ## from name involves the character rename into the columns
  key = c("mz", "q_value", "log2.fc", "rt"),
  ## note that both key and as_col must be ordered
  as_col = c("m.z", "p.value", "t.score", "r.t"),
  ion_mode = "negative", ppm = 10, p_cutoff = 0.05, db_pathway = "hsa_mfn")
{
  ## ---------------------------------------------------------------------- 
  if(is.data.frame(extra_entity)){
    df_mz_rt <- extra_entity
  }else{
    df_mz_rt <- mz_rt %>% 
      dplyr::filter(id %in% export$id) %>% 
      merge(export[, c("id", p_col)], all.x = T, by = "id") %>% 
      dplyr::as_tibble()
  }
  ## ---------------------------------------------------------------------- 
  ## note that this step is setting up to auto find and rename
  colnames(df_mz_rt) <- colnames(df_mz_rt) %>% 
    ## find and sort as order of `key`
    .meta_find_and_sort(., key) %>% 
    ## rename columns of df_mz_rt
    mapply_rename_col(., as_col, colnames(df_mz_rt))
  ## ------------------------------------- 
  ## then, select and relocate (sort) the columns of df
  df <- dplyr::select(df_mz_rt, all_of(as_col)) %>% 
    ## convert rt from min to secounds
    dplyr::mutate(r.t = r.t * 60)
  ## ---------------------------------------------------------------------- 
  ## save to file
  write_tsv(df, file = "tmp.txt")
  ## ------------------------------------- 
  ## get the submit file
  if(only_return == T)
    return(list(id = df_mz_rt, submit = df))
  ## ------------------------------------- 
  cat("## submit to MetaboAnalyst\n")
  print(dplyr::as_tibble(df))
  ## submit to MetaboAnalyst
  mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
  mSet <- SetPeakFormat(mSet, "mprt")
  mSet <- UpdateInstrumentParameters(mSet, ppm, ion_mode, "yes", 0.02);
  mSet <- Read.PeakListData(mSet, "tmp.txt");
  mSet <- SetRTincluded(mSet, "seconds")
  mSet <- SanityCheckMummichogData(mSet)
  mSet <- SetPeakEnrichMethod(mSet, "mum", "v2")
  mSet <- SetMummichogPval(mSet, p_cutoff)
  mSet <- PerformPSEA(mSet, db_pathway, "current", 3 , 100)
  return(mSet)
}


