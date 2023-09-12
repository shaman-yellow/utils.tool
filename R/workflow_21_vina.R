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
  function(x, cids, job_herb = NULL, hubs = 10){
    hgnc_symbols <- head(x@tables$step1$hub_genes$hgnc_symbol, n = hubs)
    if (!is.null(job_herb)) {
      compounds_targets <- dplyr::filter(job_herb@tables$step2$compounds_targets,
        Target.name %in% dplyr::all_of(hgnc_symbols))
      compounds <- dplyr::filter(object(job_herb)$component,
        Ingredient_id %in% compounds_targets$Ingredient_id)
      cids <- compounds$PubChem_id
      cids <- cids[ !is.na(cids) ]
    }
    x <- job_vina(cids, hgnc_symbols)
    x@params$compounds <- compounds
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
  function(x){
    step_message("")
    runs <- tibble::tibble(
      Ligand = rep(names(x$dock_layout), lengths(x$dock_layout)),
      Receptor = tolower(unlist(x$dock_layout, use.names = F))
    )
    x$show_layout <- runs
    runs <- apply(runs, 1, unname, simplify = F)
    return(x)
  })
