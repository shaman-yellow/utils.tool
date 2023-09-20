# ==========================================================================
# workflow of vina
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_vina <- setClass("job_vina", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorials: https://autodock-vina.readthedocs.io/en/latest/docking_basic.html")
    ))

job_vina <- function(cids, hgnc_symbols)
{
  .job_vina(object = namel(cids, hgnc_symbols))
}

setGeneric("asjob_vina", 
  function(x, ...) standardGeneric("asjob_vina"))

setMethod("asjob_vina", signature = c(x = "job_stringdb"),
  function(x, cids, job_herb = NULL, compounds = NULL, hubs = 10){
    hgnc_symbols <- head(x@tables$step1$hub_genes$hgnc_symbol, n = hubs)
    if (!is.null(job_herb)) {
      compounds_targets <- dplyr::filter(job_herb@tables$step2$compounds_targets,
        Target.name %in% dplyr::all_of(hgnc_symbols))
      compounds <- dplyr::filter(object(job_herb)$component,
        Ingredient_id %in% compounds_targets$Ingredient_id)
      cids <- compounds$PubChem_id
      cids <- cids[ !is.na(cids) ]
      from_job_herb <- T
    } else {
      from_job_herb <- NULL
    }
    x <- job_vina(cids, hgnc_symbols)
    x$compounds <- compounds
    x$from_job_herb <- from_job_herb
    return(x)
  })

setMethod("step0", signature = c(x = "job_vina"),
  function(x){
    step_message("Prepare your data with function `job_vina`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_vina"),
  function(x, bdb = "../BindingDB_All.tsv", each_target = 3){
    step_message("Prepare Docking Combination.
      "
    )
    if (is.null(x$mart)) {
      mart <- new_biomart()
    } else {
      mart <- x$mart
    }
    x$targets_annotation <- filter_biomart(mart, c("hgnc_symbol", "pdb"), "hgnc_symbol",
      object(x)$hgnc_symbols, distinct = F)
    x$targets_annotation <- dplyr::filter(x$targets_annotation, pdb != "")
    if (!is.null(bdb)) {
      bdb <- ld_cutRead(bdb, c("PubChem CID", "PDB ID(s) of Target Chain"))
      bdb <- dplyr::filter(bdb, `PubChem CID` %in% object(x)$cids)
      colnames(bdb) <- c("PubChem_id", "pdb_id")
      x@params$bdb_compounds_targets <- bdb
      bdb <- nl(bdb$PubChem_id, strsplit(bdb$pdb_id, ","))
      bdb <- lst_clear0(bdb)
      bdb <- lapply(bdb, function(v) v[ v %in% x$targets_annotation$pdb ])
      bdb <- lst_clear0(bdb)
      x$dock_layout <- bdb
    } else {
      pdbs <- split(x$targets_annotation$pdb, x$targets_annotation$hgnc_symbol)
      pdbs <- unlist(lapply(pdbs, head, n = each_target), use.names = F)
      x$dock_layout <- sapply(object(x)$cids, function(cid) pdbs, simplify = F)
      names(x$dock_layout) <- object(x)$cids
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_vina"),
  function(x, cl = 10){
    step_message("Download sdf files and convert as pdbqt.")
    sdfFile <- query_sdfs(unique(names(x$dock_layout)), curl_cl = cl)
    res.pdbqt <- sdf_as_pdbqts(sdfFile) 
    x$res.ligand <- nl(res.pdbqt$pdbqt.cid, res.pdbqt$pdbqt)
    return(x)
  })

setMethod("step3", signature = c(x = "job_vina"),
  function(x, cl = 10, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS"){
    step_message("Dowload pdb files of targets.")
    pdb.files <- get_pdb(unique(unlist(x$dock_layout, use.names = F)), cl = cl)
    pdb.files <- filter_pdbs(pdb.files, pattern)
    x$res.receptor <- prepare_receptor(pdb.files)
    return(x)
  })

setMethod("step4", signature = c(x = "job_vina"),
  function(x, time = 3600, savedir = "vina_space", log = "/tmp/res.log", save.object = "vn3.rds")
  {
    step_message("Run vina ...")
    runs <- tibble::tibble(
      Ligand = rep(names(x$dock_layout), lengths(x$dock_layout)),
      Receptor = tolower(unlist(x$dock_layout, use.names = F))
    )
    x$show_layout <- runs
    runs <- apply(runs, 1, unname, simplify = F)
    x$runs <- runs
    x$savedir <- savedir
    n <- 0
    if (file.exists(log)) {
      file.remove("/tmp/res.log")
    }
    saveRDS(x, save.object)
    pbapply::pblapply(x$runs,
      function(v) {
        lig <- x$res.ligand[[ v[1] ]]
        recep <- x$res.receptor[[ tolower(v[2]) ]]
        n <<- n + 1
        vina_limit(lig, recep, time, dir = savedir)
      }
    )
    return(x)
  })

setMethod("step5", signature = c(x = "job_vina"),
  function(x, compounds, by.y, facet = "Ingredient_name"){
    step_message("Summary and visualization for results.")
    x$summary_vina <- summary_vina(x$savedir)
    res_dock <- dplyr::mutate(x$summary_vina, PubChem_id = as.integer(PubChem_id))
    res_dock <- tbmerge(res_dock, dplyr::mutate(x$targets_annotation, pdb = tolower(pdb)),
      by.x = "PDB_ID", by.y = "pdb", all.x = T)
    if (!is.null(x$from_job_herb)) {
      res_dock <- tbmerge(res_dock, x$compounds, by = "PubChem_id", all.x = T)
      facet <- "Ingredient_name"
    } else {
      res_dock <- tbmerge(res_dock, compounds,
        by.x = "PubChem_id", by.y = by.y, all.x = T)
    }
    res_dock <- dplyr::arrange(res_dock, Affinity)
    data <- dplyr::distinct(res_dock, PubChem_id, hgnc_symbol, .keep_all = T)
    p.res_vina <- ggplot(data) + 
      geom_col(aes(x = hgnc_symbol, y = Affinity, fill = Affinity), width = .7) +
      labs(x = "", y = "Affinity (kcal/mol)") +
      coord_flip() +
      facet_wrap(as.formula(paste0("~ Hmisc::capitalize(paste0(", facet, "))")),
        ncol = 1, scales = "free_y") +
      theme()
    x@tables[[ 5 ]] <- namel(res_dock, unique_tops = data)
    x@plots[[ 5 ]] <- namel(p.res_vina)
    return(x)
  })

setMethod("step6", signature = c(x = "job_vina"),
  function(x, time = 2){
    step_message("Use pymol for all visualization.")
    pbapply::pbapply(x@tables$step5$res_dock, 1,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], timeLimit = time)
      }
    )
    return(x)
  })

setMethod("step7", signature = c(x = "job_vina"),
  function(x, time = 120, backup = "./figs"){
    step_message("Custom visualization.")
    pbapply::pbapply(x@tables$step5$unique_tops, 1,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], timeLimit = time, backup = backup)
      }
    )
    return(x)
  })
