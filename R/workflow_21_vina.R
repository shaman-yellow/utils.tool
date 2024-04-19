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
    cite = "[@AutodockVina1Eberha2021; @AutogridfrImpZhang2019; @AutodockCrankpZhang2019; @AutositeAnAuRavind2016; @AutodockfrAdvRavind2015]",
    method = "The CLI tools of `AutoDock vina` and `ADFR` software used for auto molecular docking"
    ))

job_vina <- function(cids, hgnc_symbols, .layout = NULL)
{
  if (missing(cids) || missing(hgnc_symbols)) {
    message("Get from `.layout` first 3 columns: cpd names, hgnc_symbols, cids")
    cids <- nl(.layout[[ 1 ]], .layout[[ 3 ]])
    hgnc_symbols <- .layout[[ 2 ]]
  }
  x <- .job_vina(object = namel(cids, hgnc_symbols))
  x$.layout <- .layout
  x
}

setGeneric("asjob_vina", 
  function(x, ...) standardGeneric("asjob_vina"))

setMethod("asjob_vina", signature = c(x = "job_stringdb"),
  function(x, cids, job_herb = NULL, compounds = NULL, hubs = 10)
  {
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
  function(x, bdb = F, each_target = 1, pdbs = NULL, bdb_file = .prefix("BindingDB_All_202401.tsv", "db"))
  {
    step_message("Prepare Docking Combination.")
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
    if (any(isThat <- !object(x)$hgnc_symbols %in% x$targets_annotation$hgnc_symbol)) {
      message("PDB not found:\n\t", paste0(unique(object(x)$hgnc_symbols[ isThat ]), collapse = ", "))
      if (!usethis::ui_yeah("Continue?")) {
        stop("Consider other ways to found PDB files for docking.")
      }
    } else {
      message("Got PDB for all `hgnc_symbol`.")
    }
    if (bdb) {
      if (file.exists(bdb_file)) {
        message("Database `BindingDB` used for pre-filtering of docking candidates.")
        bdb <- ld_cutRead(bdb_file, c("PubChem CID", "PDB ID(s) of Target Chain"))
        bdb <- dplyr::filter(bdb, `PubChem CID` %in% object(x)$cids)
        colnames(bdb) <- c("PubChem_id", "pdb_id")
        x@params$bdb_compounds_targets <- bdb
        bdb <- nl(bdb$PubChem_id, strsplit(bdb$pdb_id, ","))
        bdb <- lst_clear0(bdb)
        bdb <- lapply(bdb, function(v) v[ v %in% x$targets_annotation$pdb ])
        bdb <- lst_clear0(bdb)
        if (!length(bdb)) {
          stop("No candidates found in BindingDB.")
        }
        x$dock_layout <- bdb
        .add_internal_job(.job(method = "Database `BindingDB` used for pre-filtering of docking candidates.",
            cite = "[@BindingdbIn20Gilson2016]"))
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
  function(x, try_cluster_random = T, nGroup = 30, nMember = 3, cl = 10)
  {
    step_message("Download sdf files and convert as pdbqt for ligands.")
    sdfFile <- query_sdfs(unique(names(x$dock_layout)), curl_cl = cl)
    Show_filter <- F
    if (try_cluster_random) {
      if (length(object(x)$cids) > nGroup) {
        message("To reduce docking candidates, clustering the molecules, and random sample each group to get `n`")
        set.seed(100)
        sdfset <- e(ChemmineR::read.SDFset(sdfFile))
        message("Read ", length(sdfset), " molecules.")
        apset <- e(ChemmineR::sdf2ap(sdfset))
        cluster <- e(ChemmineR::cmp.cluster(db = apset, cutoff = seq(0.9, 0.4, by = -0.1)))
        cluster <- dplyr::select(cluster, 1, dplyr::starts_with("CLID"))
        nCl <- apply(cluster[, -1], 2, function(x) length(unique(x)))
        useWhich <- which( nCl < 30 )[1] + 1
        if (is.na(useWhich)) {
          useWhich <- menu(paste0(names(nCl), "__", unname(nCl)), title = "Need custom select for cutoff.")
          useWhich <- useWhich + 1
        }
        message("Use cutoff: ", colnames(cluster)[useWhich])
        groups <- split(seq_along(sdfset), cluster[useWhich])
        message("Now, get group number: ", length(groups))
        groups <- lapply(groups,
          function(x) {
            if (length(x) > nMember) {
              sample(x, nMember)
            } else x
          })
        sdfset <- sdfset[unlist(groups, use.names = F)]
        MaybeError <- vapply(ChemmineR::cid(sdfset), FUN.VALUE = logical(1),
          function(x) {
            dim(sdfset[[x]][[2]])[2] < 3
        })
        sdfset <- sdfset[which(!MaybeError)]
        sdfFile <- gs(sdfFile, "\\.sdf$", "_random.sdf")
        e(ChemmineR::write.SDF(sdfset, sdfFile))
        Show_filter <- T
        .add_internal_job(.job(method = "R package `ChemmineR` used for similar chemical compounds clustering",
            cite = "[@ChemminerACoCaoY2008]"
            ))
        x$chem_lich <- new_lich(list(
            "Clustering method:" = "Binning Clustering",
            "Use Cut-off:" = colnames(cluster)[useWhich],
            "Cluster number:" = length(groups),
            "Each sampling for next step (docking):" = nMember
            ))
      }
    }
    sdfFile <- cal_3d_sdf(sdfFile)
    if (Show_filter) {
      message("Filter out (Due to Chemmine bug): ", length(which(MaybeError)))
      message("Now, Total docking molecules: ", length(sdfset))
    }
    res.pdbqt <- mk_prepare_ligand.sdf(sdfFile)
    x$res.ligand <- nl(res.pdbqt$pdbqt.cid, res.pdbqt$pdbqt)
    message("Got (filter out in `mk_prepare_ligand.sdf`): ", length(x$res.ligand))
    alls <- unique(as.character(object(x)$cid))
    message("Not got:", paste0(alls[!alls %in% names(x$res.ligand)], collapse = ", "))
    message("Filter the `x$dock_layout`")
    x$dock_layout <- x$dock_layout[ names(x$dock_layout) %in% names(x$res.ligand) ]
    return(x)
  })

setMethod("step3", signature = c(x = "job_vina"),
  function(x, cl = 10, pattern = NULL,
    extra_pdb.files = NULL, extra_layouts = NULL, extra_symbols = NULL, filter = T)
  {
    step_message("Dowload pdb files for Receptors.")
    ids <- unique(unlist(x$dock_layout, use.names = F))
    if (length(ids)) {
      pdb.files <- get_pdb(ids, cl = cl)
    } else {
      pdb.files <- NULL
    }
    if (!is.null(extra_pdb.files)) {
      names(extra_pdb.files) <- tolower(names(extra_pdb.files))
      pdb.files <- c(pdb.files, extra_pdb.files)
      if (is.null(extra_symbols)) {
        extra_symbols <- nl(toupper(names(extra_pdb.files)), names(extra_pdb.files), F)
      }
      if (!nrow(x@params$targets_annotation)) {
        x@params$targets_annotation <- dplyr::mutate(x@params$targets_annotation, pdb = character(0))
      }
      x@params$targets_annotation %<>%
        tibble::add_row(hgnc_symbol = names(extra_symbols), pdb = unname(extra_symbols))
      layoutNeedRevise <- T
    } else {
      layoutNeedRevise <- F
    }
    if (!is.null(extra_layouts)) {
      if (is(extra_layouts, "character")) {
        extra_layouts <- as.list(extra_layouts)
      }
      x$dock_layout <- c(x$dock_layout, extra_layouts)
    } else if (layoutNeedRevise) {
      x$dock_layout %<>% lapply(function(x) c(x, names(extra_pdb.files)))
    }
    if (!is.null(pattern)) {
      pdb.files <- filter_pdbs(pdb.files, pattern)
    }
    x$res.receptor <- prepare_receptor(pdb.files)
    names <- names(x$res.receptor)
    gotSymbols <- x$targets_annotation$hgnc_symbol[ match(tolower(names), tolower(x$targets_annotation$pdb)) ]
    x$res.receptor.symbol <- gotSymbols
    fun <- function(x) x[ !x %in% gotSymbols ]
    message("Not got: ", paste0(fun(unique(object(x)$hgnc_symbol)), collapse = ", "))
    if (filter) {
      if (!is.null(x$.layout)) {
        message("Customize using `x$.layout` columns: ", paste0(colnames(x$.layout)[1:2], collapse = ", "))
        x <- filter(x, x$.layout[[ 1 ]], x$.layout[[ 2 ]])
      }
    }
    return(x)
  })

setMethod("filter", signature = c(x = "job_vina"),
  function(x, cpd, symbol, cid = NULL)
  {
    message("Custom specified docking.")
    if (x@step != 3L) {
      stop("x@step != 3L")
    }
    if (is.null(cid)) {
      cid <- unname(object(x)$cids[ match(cpd, names(object(x)$cids)) ])
    }
    ref <- dplyr::distinct(x$targets_annotation, hgnc_symbol, .keep_all = T)
    pdb <- ref$pdb[ match(symbol, ref$hgnc_symbol) ]
    layout <- nl(cid, pdb)
    x$dock_layout <- layout
    return(x)
  })

setMethod("step4", signature = c(x = "job_vina"),
  function(x, time = 3600 * 2, savedir = "vina_space", log = "/tmp/res.log", save.object = "vn3.rds")
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
        if (!is.null(lig) & !is.null(recep)) {
          vina_limit(lig, recep, time, dir = savedir, x = x)
        }
      }
    )
    if (.Platform$OS.type == "unix") {
      system("notify-send 'AutoDock vina' 'All job complete'")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_vina"),
  function(x, compounds, by.y, facet = "Ingredient_name", excludes = NULL, top = 5, cutoff.af = NULL)
  {
    step_message("Summary and visualization for results.")
    x$summary_vina <- summary_vina(x$savedir)
    x$summary_vina <- dplyr::filter(x$summary_vina,
      PubChem_id %in% !!names(x$dock_layout),
      PDB_ID %in% tolower(unlist(x$dock_layout, use.names = F))
    )
    if (!is.null(top)) {
      x$summary_vina <- split_lapply_rbind(x$summary_vina, ~ PDB_ID,
        function(x) {
          dplyr::slice_min(x, Affinity, n = top)
        }
      )
    }
    res_dock <- dplyr::mutate(x$summary_vina, PubChem_id = as.integer(PubChem_id))
    anno <- dplyr::mutate(x$targets_annotation, pdb = tolower(pdb))
    res_dock <- map(res_dock, "PDB_ID", anno, "pdb", "hgnc_symbol", col = "hgnc_symbol")
    if (!is.null(x$from_job_herb)) {
      res_dock <- map(res_dock, "PubChem_id",
        x$compounds, "PubChem_id", "Ingredient_name", col = "Ingredient_name")
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
      res_dock <- map(res_dock, "PubChem_id",
        compounds, "PubChem_id", "Ingredient_name", col = "Ingredient_name")
    }
    res_dock <- dplyr::arrange(res_dock, Affinity)
    res_dock <- dplyr::filter(res_dock, !is.na(!!rlang::sym(facet)))
    data <- dplyr::distinct(res_dock, PubChem_id, hgnc_symbol, .keep_all = T)
    if (!is.null(cutoff.af)) {
      data <- dplyr::filter(data, Affinity < !!cutoff.af)
    }
    if (!is.null(excludes)) {
      data <- dplyr::filter(data, !hgnc_symbol %in% !!excludes)
    }
    if (nrow(data) > 20) {
      data <- dplyr::filter(data, Affinity < 0)
    }
    if (nrow(data) > 20) {
      data <- head(data, n = 20)
    }
    data <- dplyr::mutate(data, dplyr::across(!!rlang::sym(facet), function(x) stringr::str_trunc(x, 30)))
    p.res_vina <- ggplot(data) + 
      geom_col(aes(x = reorder(hgnc_symbol, Affinity, decreasing = T), y = Affinity, fill = Affinity), width = .7) +
      geom_text(data = dplyr::filter(data, Affinity <= 0),
        aes(x = hgnc_symbol, y = Affinity - .5, label = round(Affinity, 1)), hjust = 1) +
      labs(x = "", y = "Affinity (kcal/mol)") +
      coord_flip() +
      ylim(zoRange(c(-1, data$Affinity, 1), 1.4)) +
      facet_wrap(as.formula(paste0("~ Hmisc::capitalize(paste0(", facet, ", \" (CID: \", PubChem_id, \")\"))")),
        ncol = 1, scales = "free_y") +
      theme()
    p.res_vina <- wrap(p.res_vina, 7, nrow(data) + 1)
    p.res_vina <- .set_lab(p.res_vina, sig(x), "Overall combining Affinity")
    x@tables[[ 5 ]] <- namel(res_dock, unique_tops = data)
    x@plots[[ 5 ]] <- namel(p.res_vina)
    return(x)
  })

setMethod("step6", signature = c(x = "job_vina"),
  function(x, time = 15, top = NULL, save = T, unique = F, symbol = NULL, cpd = NULL)
  {
    step_message("Use pymol for all visualization.")
    data <- dplyr::filter(x@tables$step5$res_dock, Affinity < 0)
    if (!is.null(symbol)) {
      data <- dplyr::filter(data, hgnc_symbol %in% !!symbol)
    }
    if (!is.null(cpd)) {
      data <- dplyr::filter(data, Ingredient_name %in% !!cpd)
    }
    if (unique) {
      data <- dplyr::distinct(data, hgnc_symbol, .keep_all = T)
    }
    if (!is.null(top)) {
      data <- head(data, n = top)
    }
    figs <- pbapply::pbapply(data, 1, simplify = F,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], timeLimit = time, save = save)
      }
    )
    figs <- unlist(figs, recursive = F)
    lab(figs) <- paste(sig(x), "docking visualization")
    names(figs) <- paste0("Top", seq_along(figs), "_", names(figs))
    x@plots[[ 6 ]] <- figs
    data <- .set_lab(data, sig(x), "Metadata of visualized Docking")
    data <- dplyr::select(data, -dir, -file)
    x@tables[[ 6 ]] <- namel(data)
    return(x)
  })

setMethod("step7", signature = c(x = "job_vina"),
  function(x, save = T){
    step_message("Show docking results in deep and in detail.")
    data <- x@tables$step6$data
    figs <- pbapply::pbapply(data, 1, simplify = F,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], save = save, detail = T)
      }
    )
    figs <- unlist(figs, recursive = F)
    lab(figs) <- paste(sig(x), "docking interaction details")
    names(figs) <- paste0("Top", seq_along(figs), "_", names(figs))
    x@plots[[ 7 ]] <- figs
    return(x)
  })

pretty_docking <- function(protein, ligand, path,
  save = "annotation.png",
  script = paste0(.expath, "/pretty_docking.pymol"))
{
  script <- readLines(script)
  script <- gs(script, "{{protein.pdbqt}}", protein, fixed = T)
  script <- gs(script, "{{ligand.pdbqt}}", ligand, fixed = T)
  temp <- tempfile("Pymol", fileext = ".pml")
  writeLines(script, temp)
  cli::cli_alert_info(paste0("Pymol run script: ", temp))
  output <- paste0(path, "/", save)
  message("Save png to ", output)
  expr <- paste0(" png ", save, ",2000,2000,dpi=300")
  gett(expr)
  cdRun(pg("pymol"), " ",
    " -d \"run ", temp, "\"",
    " -d \"ray; ", expr, "\"",
    path = path)
  return(output)
}

vina_limit <- function(lig, recep, timeLimit = 120, dir = "vina_space", ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit, dir = dir), T)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoring = "ad4", stout = "/tmp/res.log", timeLimit = 60,
  x)
{
  remote <- F
  if (!missing(x)) {
    if (is.remote(x)) {
      remote <- T
    }
  }
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
    .cdRun(pg("prepare_gpf.py"), " -l ", files[1], " -r ", files[2], " -y")
    .cdRun(pg("autogrid4"), " -p ", reals[2], ".gpf ", " -l ", reals[2], ".glg")
    if (remote) {
      message("Run in remote server.")
      cdRun("scp -r ", subdir, " ", x$remote, ":", x$wd, path = dir)
      remoteWD <- x$wd
      x$wd <- paste0(x$wd, "/", subdir)
      output <- paste0(subdir, "_out.pdbqt")
      rem_run("timeout ", timeLimit, 
        " vina  --ligand ", files[1],
        " --maps ", reals[2],
        " --scoring ", scoring,
        " --exhaustiveness ", exhaustiveness,
        " --out ", output, 
        " >> ", stout)
      try(get_file_from_remote(output, x$wd, paste0(wd, "/", output)))
    } else {
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
}

vinaShow <- function(Combn, recep, subdir = Combn, dir = "vina_space",
  timeLimit = 3, backup = NULL, save = T, detail = F)
{
  if (!file.exists(path <- paste0(dir, "/", subdir))) {
    stop("file.exists(path <- paste0(dir, \"/\", subdir))")
  } 
  wd <- paste0(dir, "/", subdir)
  out <- paste0(Combn, "_out.pdbqt")
  recep <- paste0(recep, ".pdbqt")
  if (!file.exists(recep)) {
    maybethat <- list.files(wd, recep, ignore.case = T)
    if (length(maybethat)) {
      recep <- maybethat[1]
    }
  }
  res <- paste0(Combn, ".png")
  if (detail) {
    res <- paste0("detail_", res)
  }
  .cdRun <- function(...) cdRun(..., path = wd)
  img <- paste0(wd, "/", res)
  if (file.exists(img)) {
    file.remove(img)
  }
  if (detail) {
    pretty_docking(recep, out, wd, save = res)
  } else {
    expr <- paste0(" png ", res, ",2000,2000,dpi=300")
    gett(expr)
    if (!save) {
      expr <- ""
    }
    try(.cdRun("timeout ", timeLimit, 
        " pymol ",
        " -d \"load ", out, ";",
        " load ", recep, ";",
        " ray;", expr, "\" "), T)
  }
  if (is.character(backup)) {
    dir.create(backup, F)
    file.copy(img, backup, T)
  }
  fig <- file_fig(Combn, img)
  lab(fig[[1]]) <- paste0("docking ", Combn, if (!detail) NULL else " detail")
  return(fig)
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
    PDB_ID = tolower(stringr::str_extract(Combn, "[^_]{1,}$")),
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
    message("Find columns of: ", paste0(names[ names %in% cols ], collapse = ", "))
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

mk_prepare_ligand.sdf <- function(sdf_file, mkdir.pdbqt = "pdbqt", check = F) {
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
  cdRun(pg("mk_prepare_ligand.py"), " -i ", sdf_file,
    " --multimol_outdir ", mkdir.pdbqt)
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

setMethod("set_remote", signature = c(x = "job_vina"),
  function(x, wd = "/data/hlc/vina", remote = "remote")
  {
    x$wd <- wd
    x$set_remote <- T
    x$remote <- remote
    return(x)
  })

cal_3d_sdf <- function(sdf) {
  output <- gs(sdf, "\\.sdf$", "_3D.sdf")
  isThat <- T
  if (file.exists(output)) {
    isThat <- usethis::ui_yeah("Overwrite the exists file?")
  }
  if (isThat) {
    cdRun(pg("obgen"), " ", sdf, " -ff UFF > ", output)
  }
  output
}

setMethod("res", signature = c(x = "job_vina"),
  function(x, meta,
    use = "Ingredient_name", target = "Ingredient.name", get = "Herb_pinyin_name")
  {
    data <- dplyr::select(x@tables$step5$res_dock, -dir, -file)
    data <- dplyr::relocate(data, hgnc_symbol, Ingredient_name, Affinity)
    if (!missing(meta)) {
      meta <- dplyr::filter(meta, !!rlang::sym(target) %in% !!data[[ use ]])
      meta <- dplyr::distinct(meta, !!rlang::sym(target), !!rlang::sym(get))
      meta <- dplyr::group_by(meta, !!rlang::sym(target))
      fun <- function(x) paste0(x, collapse = "; ")
      meta <- dplyr::reframe(meta, .get = fun(!!rlang::sym(get)))
      data <- map(data, use, meta, target, ".get", col = get)
    }
    data
  })
