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
    info = c("Tutorials: https://autodock-vina.readthedocs.io/en/latest/docking_basic.html"),
    cite = "[@AutodockVina1Eberha2021]",
    method = "`AutoDock vina` used for molecular docking"
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
  function(x, bdb = "../BindingDB_All.tsv", each_target = 3, pdbs = NULL){
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
    if (!is.null(pdbs)) {
      pdbs <- list(hgnc_symbol = names(pdbs), pdb = unname(pdbs))
      if (!is.character(x$targets_annotation$pdb)) {
        x$targets_annotation <- dplyr::mutate(x$targets_annotation, pdb = as.character(pdb))
      }
      x$targets_annotation <- tibble::add_row(x$targets_annotation, !!!pdbs)
      x$targets_annotation <- dplyr::distinct(x$targets_annotation, hgnc_symbol, .keep_all = T)
    }
    if (!is.null(bdb)) {
      if (file.exists(bdb)) {
        bdb <- ld_cutRead(bdb, c("PubChem CID", "PDB ID(s) of Target Chain"))
        bdb <- dplyr::filter(bdb, `PubChem CID` %in% object(x)$cids)
        colnames(bdb) <- c("PubChem_id", "pdb_id")
        x@params$bdb_compounds_targets <- bdb
        bdb <- nl(bdb$PubChem_id, strsplit(bdb$pdb_id, ","))
        bdb <- lst_clear0(bdb)
        bdb <- lapply(bdb, function(v) v[ v %in% x$targets_annotation$pdb ])
        bdb <- lst_clear0(bdb)
        x$dock_layout <- bdb
      }
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
    step_message("Download sdf files and convert as pdbqt for ligands.")
    sdfFile <- query_sdfs(unique(names(x$dock_layout)), curl_cl = cl)
    res.pdbqt <- sdf_as_pdbqts(sdfFile) 
    x$res.ligand <- nl(res.pdbqt$pdbqt.cid, res.pdbqt$pdbqt)
    return(x)
  })

setMethod("step3", signature = c(x = "job_vina"),
  function(x, cl = 10, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS",
    extra_pdb.files = NULL, extra_layouts = NULL, extra_symbols = NULL)
  {
    step_message("Dowload pdb files for Receptors.")
    ids <- unique(unlist(x$dock_layout, use.names = F))
    if (length(ids)) {
      pdb.files <- get_pdb(ids, cl = cl)
    } else {
      pdb.files <- NULL
    }
    if (!is.null(extra_layouts)) {
      if (is(extra_layouts, "character")) {
        extra_layouts <- as.list(extra_layouts)
      }
      x$dock_layout <- c(x$dock_layout, extra_layouts)
    }
    if (!is.null(extra_pdb.files)) {
      pdb.files <- c(pdb.files, extra_pdb.files)
      if (is.null(extra_symbols)) {
        extra_symbols <- nl(toupper(names(extra_pdb.files)), names(extra_pdb.files), F)
      }
      if (!nrow(x@params$targets_annotation)) {
        x@params$targets_annotation <- dplyr::mutate(x@params$targets_annotation, pdb = character(0))
      }
      x@params$targets_annotation %<>%
        tibble::add_row(hgnc_symbol = names(extra_symbols), pdb = unname(extra_symbols))
    }
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
      if (missing(compounds)) {
        compounds <- data.frame(PubChem_id = as.integer(object(x)$cids))
        if (is.null(names(object(x)$cids))) {
          compounds$Ingredient_name <- paste0("CID:", object(x)$cids)
        } else {
          compounds$Ingredient_name <- names(object(x)$cids)
        }
        by.y <- "PubChem_id"
        facet <- "Ingredient_name"
      }
      res_dock <- tbmerge(res_dock, compounds,
        by.x = "PubChem_id", by.y = by.y, all.x = T)
    }
    res_dock <- dplyr::arrange(res_dock, Affinity)
    data <- dplyr::distinct(res_dock, PubChem_id, hgnc_symbol, .keep_all = T)
    p.res_vina <- ggplot(data) + 
      geom_col(aes(x = reorder(hgnc_symbol, Affinity, decreasing = T), y = Affinity, fill = Affinity), width = .7) +
      labs(x = "", y = "Affinity (kcal/mol)") +
      coord_flip() +
      ylim(zoRange(c(-1, data$Affinity), 1.4)) +
      facet_wrap(as.formula(paste0("~ Hmisc::capitalize(paste0(", facet, "))")),
        ncol = 1, scales = "free_y") +
      theme()
    p.res_vina <- wrap(p.res_vina)
    p.res_vina <- .set_lab(p.res_vina, sig(x), "Overall combining Affinity")
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

vina_limit <- function(lig, recep, timeLimit = 120, dir = "vina_space", ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit, dir = dir), T)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoring = "ad4", stout = "/tmp/res.log", timeLimit = 60)
{
  if (!file.exists(dir)) {
    dir.create(dir, F)
  } 
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
  if (!file.exists(paste0(wd, "/", subdir, "_out.pdbqt"))) {
    dir.create(wd, F)
    file.copy(c(recep, lig), wd)
    .message_info("Generating affinity maps", subdir)
    .cdRun <- function(...) cdRun(..., path = wd)
    files <- get_filename(c(lig, recep))
    .cdRun("prepare_gpf.py -l ", files[1], " -r ", files[2], " -y")
    .cdRun("autogrid4 -p ", reals[2], ".gpf ", " -l ", reals[2], ".glg")
    cat("\n$$$$\n", date(), "\n", subdir, "\n\n", file = stout, append = T)
    try(.cdRun("timeout ", timeLimit, 
        " vina  --ligand ", files[1],
        " --maps ", reals[2],
        " --scoring ", scoring,
        " --exhaustiveness ", exhaustiveness,
        " --out ", subdir, "_out.pdbqt",
        " >> ", stout), T)
  }
}

vinaShow <- function(Combn, recep, subdir = Combn, dir = "vina_space",
  timeLimit = 3, backup = NULL)
{
  if (!file.exists(path <- paste0(dir, "/", subdir))) {
    stop("file.exists(path <- paste0(dir, \"/\", subdir))")
  } 
  wd <- paste0(dir, "/", subdir)
  out <- paste0(Combn, "_out.pdbqt")
  recep <- paste0(recep, ".pdbqt")
  res <- paste0(Combn, ".png")
  .cdRun <- function(...) cdRun(..., path = wd)
  try(.cdRun("timeout ", timeLimit, 
      " pymol ", out, " ", recep,
      " -g ", res), T)
  if (is.character(backup)) {
    dir.create(backup, F)
    file.copy(paste0(wd, "/", res), backup, T)
  }
}

summary_vina <- function(space = "vina_space", pattern = "_out\\.pdbqt$")
{
  files <- list.files(space, pattern, recursive = T, full.names = T)
  res_dock <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (length(lines) >= 1) {
        name <- gsub("_out\\.pdbqt", "", get_filename(file))
        top <- stringr::str_extract(lines[2], "[\\-0-9.]{1,}")
        top <- as.double(top)
        names(top) <- name
        top
      }
    })
  res_dock <- unlist(res_dock)
  res_dock <- tibble::tibble(
    Combn = names(res_dock), Affinity = unname(res_dock)
  )
  res_dock <- dplyr::mutate(
    res_dock, PubChem_id = stringr::str_extract(Combn, "^[^_]{1,}"),
    PDB_ID = stringr::str_extract(Combn, "[^_]{1,}$"),
    dir = paste0(space, "/", Combn),
    file = paste0(dir, "/", Combn, "_out.pdbqt")
  )
  dplyr::select(res_dock, PubChem_id, PDB_ID, Affinity, dir, file, Combn)
}

smiles_as_sdfs.obabel <- function(smiles) {
  lst.sdf <- pbapply::pbsapply(smiles, simplify = F,
    function(smi) {
      ChemmineOB::convertFormat("SMI", "SDF", source = test)
    })
  lst.sdf
}

tbmerge <- function(x, y, ...) {
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)
  tibble::as_tibble(data.table::merge.data.table(x, y, ...))
}

ld_cols <- function(file, sep = "\t") {
  line <- readLines(file, n = 1)
  strsplit(line, sep)[[ 1 ]]
}

lst_clear0 <- function(lst, len = 0) {
  lst[ vapply(lst, function(v) if (length(v) > len) T else F, logical(1)) ]
}

ld_cutRead <- function(file, cols, abnum = T, sep = "\t", tmp = "/tmp/ldtmp.txt") {
  if (is.character(cols)) {
    names <- ld_cols(file, sep)
    pos <- which( names %in% cols )
    if (abnum) {
      pos <- head(pos, n = length(cols))
    }
  } else {
    pos <- cols
  }
  cdRun("cut -f ", paste0(pos, collapse = ","), " ", file, " > ", tmp)
  ftibble(tmp)
}

sdf_as_pdbqts <- function(sdf_file, mkdir.pdbqt = "pdbqt", check = F) {
  dir.create(mkdir.pdbqt, F)
  check_sdf_validity <- function(file) {
    lst <- sep_list(readLines(file), "^\\${4,}$")
    lst <- lst[ - length(lst) ]
    osum <- length(lst)
    lst <- lapply(lst,
      function(line) {
        pos <- grep("^\\s*-OEChem|M\\s*END", line)
        if (pos[2] - pos[1] > 5)
          line
      })
    lst <- lst[ !vapply(lst, is.null, logical(1)) ]
    nsum <- length(lst)
    list(data = lst, osum = osum, nsum = nsum, dec = osum - nsum)
  }
  if (check) {
    lst <- check_sdf_validity(sdf_file)
    lines <- unlist(lst$data)
    lst <- lst[ -1 ]
    writeLines(lines, sdf_file <- gsub("\\.sdf", "_modified.sdf", sdf_file))
  } else {
    lst <- list(file = sdf_file)
  }
  system(paste0("mk_prepare_ligand.py -i ", sdf_file, " --multimol_outdir ", mkdir.pdbqt))
  lst$file <- sdf_file
  lst$pdbqt <- list.files(mkdir.pdbqt, "\\.pdbqt$", full.names = T)
  lst$pdbqt.num <- length(lst$pdbqt)
  lst$pdbqt.cid <- stringr::str_extract(lst$pdbqt, "(?<=/|^)[0-9]{1,}")
  return(lst)
}

select_files_by_grep <- function(files, pattern){
  res <- lapply(files,
    function(file) {
      line <- readLines(file)
      if (any(grepl(pattern, line, T)))
        return(file)
      else NULL
    })
  res <- res[ !vapply(res, is.null, logical(1)) ]
  res
}

filter_pdbs <- function(files, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS") {
  select_files_by_grep(files, pattern)
}

prepare_receptor <- function(files, mkdir.pdbqt = "protein_pdbqt") {
  dir.create(mkdir.pdbqt, F)
  file <- lapply(files,
    function(file) {
      if (!is.null(file)) {
        newfile <- gsub("\\.pdb$", ".pdbqt", get_filename(file))
        newfile <- paste0(mkdir.pdbqt, "/", newfile)
        system(paste0("prepare_receptor -r ", file, " -o ", newfile, " -A hydrogens"))
        return(newfile)
      }
    })
  file <- file[ !vapply(file, is.null, logical(1)) ]
  file[ vapply(file, file.exists, logical(1)) ]
}

