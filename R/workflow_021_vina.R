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
    method = "The CLI tools of `AutoDock vina` and `ADFR` software used for auto molecular docking",
    tag = "dock:vina",
    analysis = "AutoDock vina 分子对接"
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

setGeneric("asjob_vina", group = list("asjob_series"),
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
      from_job_herb <- TRUE
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
    step_message("Prepare your data with function `job_vina`. ")
  })

setMethod("step1", signature = c(x = "job_vina"),
  function(x, order = TRUE, each_target = 1, custom_pdbs = NULL, 
    bdb_file = .prefix("BindingDB_All_202401.tsv", "db"))
  {
    step_message("Prepare Docking Combination.")
    x <- .find_proper_pdb(x, order, each_target, custom_pdbs)
    x$dock_layout <- sapply(object(x)$cids, function(cid) x$used_pdbs, simplify = FALSE)
    names(x$dock_layout) <- object(x)$cids
    return(x)
  })

.find_proper_pdb <- function(x, order = TRUE, each_target = 1, custom_pdbs = NULL) {
  if (!is(x, "job")) {
    stop('!is(x, "job").')
  }
  if (!is(object(x)$hgnc_symbols, "character")) {
    stop('!is(object(x)$hgnc_symbols, "character").')
  }
  if (is.null(x$mart)) {
    mart <- new_biomart()
  } else {
    mart <- x$mart
  }
  x$targets_annotation <- filter_biomart(mart, c("hgnc_symbol", "pdb"), "hgnc_symbol",
    object(x)$hgnc_symbols, distinct = FALSE)
  if (nrow(x$targets_annotation)) {
    x <- methodAdd(x, "以 R 包 `biomaRt` ({packageVersion('biomaRt')}) {cite_show('MappingIdentifDurinc2009')} 获取基因 Symbol 对应的蛋白结构 PDB (<https://www.rcsb.org/>) 数据库 ID。")
    x <- snapAdd(x, "以 `biomaRt` 获取基因 Symbol 对应的蛋白结构 (PDB，详见方法章节)。")
    x$targets_annotation <- dplyr::filter(x$targets_annotation, pdb != "")
    x$targets_annotation <- dplyr::distinct(x$targets_annotation, pdb, .keep_all = TRUE)
  }
  if (!is.null(custom_pdbs)) {
    custom_pdbs <- list(hgnc_symbol = names(custom_pdbs), pdb = unname(custom_pdbs))
    if (!is.character(x$targets_annotation$pdb)) {
      x$targets_annotation <- dplyr::mutate(x$targets_annotation, pdb = as.character(pdb))
    }
    x$targets_annotation <- tibble::add_row(x$targets_annotation, !!!custom_pdbs)
  }
  if (any(isThat <- !object(x)$hgnc_symbols %in% x$targets_annotation$hgnc_symbol)) {
    x$pdb_notGot <- unique(object(x)$hgnc_symbols[ isThat ])
    message("PDB not found:\n\t", paste0(x$pdb_notGot, collapse = ", "))
    if (!sureThat("Continue? (Suggests: `.get_pdb_files`)"))
    {
      stop("Consider other ways to found PDB files for docking.")
    }
  } else {
    message("Got PDB for all `hgnc_symbol`.")
  }
  if (nrow(x$targets_annotation)) {
    x$annoPdbs <- as_tibble(e(bio3d::pdb.annotate(x$targets_annotation$pdb)))
    x <- methodAdd(
      x, "以 R 包 `bio3d` ({packageVersion('bio3d')}) 获取 PDB ID 对应的注释 (蛋白结构分辨率, resolution) 。"
    )
    if (order) {
      annoPdbs <- dplyr::distinct(x$annoPdbs, structureId, resolution)
      x$targets_annotation <- map(
        x$targets_annotation, "pdb", annoPdbs, "structureId", 
        "resolution", col = "resolution"
      )
      x$targets_annotation <- dplyr::arrange(x$targets_annotation, resolution)
      x <- methodAdd(x, "首要以 resolution 选取用于分子对接的蛋白结构 (resolution 越小，分辨率越高) 。")
      x <- snapAdd(x, "选取分辨率最高 (即，resolution 值最小) 的 PDB 作为分子对接的蛋白结构。")
    }
    used_pdbs <- dplyr::distinct(
      x$targets_annotation, hgnc_symbol, .keep_all = TRUE
    )
    isMultiChains <- vapply(used_pdbs$pdb,
      function(pdb) {
        length(which(x$annoPdbs$structureId == pdb)) > 1
      }, logical(1))
    if (any(isMultiChains)) {
      x$pdb_MultiChains <- used_pdbs$pdb[ isMultiChains ]
      message(glue::glue("Some 'pdb' has multiple chains ({bind(x$pdb_MultiChains)})"))
    }
    used_pdbs <- setLegend(used_pdbs, "基因 Symbol 所用的蛋白结构，以及对应的 PDB ID 和分辨率。")
    x <- tablesAdd(x, t.proteins_used_PDB = used_pdbs)
    x$used_pdbs <- used_pdbs <- nl(
      used_pdbs$hgnc_symbol, used_pdbs$pdb, FALSE
    )
    if (!length(used_pdbs)) {
      stop('!length(used_pdbs), no any protein need docking?')
    }
  } else {
    x$used_pdbs <- NULL
  }
  return(x)
}

setMethod("step2", signature = c(x = "job_vina"),
  function(x, try_cluster_random = TRUE, nGroup = 100, nMember = 3, cl = 10, sdf.3d = NULL)
  {
    step_message("Download sdf files and convert as pdbqt for ligands.")
    sdfFile <- query_sdfs(unique(names(x$dock_layout)), curl_cl = cl)
    x <- methodAdd(x, "以 PubChem API
      (<https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>) 获取化合物 SDF
      结构文件。"
    )
    x <- snapAdd(x, "从 PubChem 获取化合物 SDF 结构文件(2D)。")
    Show_filter <- FALSE
    if (try_cluster_random && length(object(x)$cids) > nGroup) {
      message("To reduce docking candidates, clustering the molecules, and random sample each group to get `n`")
      set.seed(100)
      sdfset <- e(ChemmineR::read.SDFset(sdfFile))
      message("Read ", length(sdfset), " molecules.")
      apset <- e(ChemmineR::sdf2ap(sdfset))
      cluster <- e(ChemmineR::cmp.cluster(db = apset, cutoff = seq(0.9, 0.4, by = -0.1)))
      x <- methodAdd(x, "以 R 包 `ChemmineR` ({packageVersion('ChemmineR')}) {cite_show('ChemminerACoCaoY2008')} 计算化合物结构相似度。")
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
      sdfset <- sdfset[unlist(groups, use.names = FALSE)]
      MaybeError <- vapply(ChemmineR::cid(sdfset), FUN.VALUE = logical(1),
        function(x) {
          dim(sdfset[[x]][[2]])[2] < 3
        })
      sdfset <- sdfset[which(!MaybeError)]
      sdfFile <- gs(sdfFile, "\\.sdf$", "_random.sdf")
      e(ChemmineR::write.SDF(sdfset, sdfFile))
      Show_filter <- TRUE
      .add_internal_job(.job(method = "R package `ChemmineR` used for similar chemical compounds clustering",
          cite = "[@ChemminerACoCaoY2008]"
          ))
      x$chem_lich <- new_lich(list(
          "Clustering method:" = "Binning Clustering",
          "Use Cut-off:" = colnames(cluster)[useWhich],
          "Cluster number:" = length(groups),
          "Each sampling for next step (docking):" = nMember
          )
      )
    }
    if (is.null(sdf.3d)) {
      message("Overwrite exists 3D SDF file.")
      guess_sdf.3d <- gs(sdfFile, "\\.sdf$", "_3D.sdf")
      if (file.exists(guess_sdf.3d)) {
        if (sureThat("File exists: {guess_sdf.3d}, use it?")) {
          sdfFile <- guess_sdf.3d
        } else {
          sdfFile <- cal_3d_sdf(sdfFile)
        }
      } else {
        sdfFile <- cal_3d_sdf(sdfFile)
      }
      x <- methodAdd(x, "使用 `openbabel` 的工具 (`obgen`) 计算 SDF 文件的 3D 构象 (转化为 3D SDF文件)。")
      x <- snapAdd(x, "以 `openbabel` 计算化合物的 3D 构象。")
    } else if (file.exists(sdf.3d)) {
      message("Use '", sdf.3d, "'")
      sdfFile <- sdf.3d
    } else {
      stop("file.exists(sdf.3d) == FALSE")
    }
    if (Show_filter) {
      message("Filter out (Due to Chemmine bug): ", length(which(MaybeError)))
      message("Now, Total docking molecules: ", length(sdfset))
    }
    res.pdbqt <- mk_prepare_ligand.sdf(sdfFile)
    x <- methodAdd(x, "以 Python `meeko` 包 (`mk_prepare_ligand.py`) 转化 SDF 文件获取配体 PDBQT 用于分子对接。")
    x <- snapAdd(x, "以 `meeko` 从 SDF 转化得到配体的 PDBQT 文件。")
    x$res.ligand <- nl(res.pdbqt$pdbqt.cid, res.pdbqt$pdbqt)
    message("Got (filter out in `mk_prepare_ligand.sdf`): ", length(x$res.ligand))
    alls <- unique(as.character(object(x)$cid))
    notGot <- alls[!alls %in% names(x$res.ligand)]
    message("Not got:", paste0(notGot, collapse = ", "))
    x$ligand_notGot <- notGot
    message("Filter the `x$dock_layout`")
    x$dock_layout <- x$dock_layout[ names(x$dock_layout) %in% names(x$res.ligand) ]
    return(x)
  })

setMethod("step3", signature = c(x = "job_vina"),
  function(x, cl = 10, pattern = NULL, extra_pdb.files = NULL, extra_layouts = NULL,
    filter = TRUE, use_complex = TRUE, exclude_nonStd = c(
      "NAG", "BMA", "FUL"
    ), tryGetFromAlphaFold = TRUE)
  {
    step_message("Dowload pdb files for Receptors.")
    res <- .get_pdb_files(
      x, cl = cl, pattern = pattern, extra_pdb.files = extra_pdb.files, 
      extra_layouts = extra_layouts,
      tryGetFromAlphaFold = tryGetFromAlphaFold
    )
    x <- res$x
    pdb.files <- res$pdb.files
    used_pdbs <- x$used_pdbs
    if (!is.null(exclude_nonStd)) {
      pattern <- paste0(
        paste0("\\b", exclude_nonStd, "\\b"), collapse = "|"
      )
      files_nonStds <- pdb.files[select_files_by_grep(pdb.files, pattern)]
      cmdRmd <- paste0(
        "remove ", paste0(
          glue::glue("resn {exclude_nonStd}"), collapse = " or "
        )
      )
      lapply(files_nonStds, 
        function(file) {
          cdRun(glue::glue("{pg('pymol')} -c -Q -d 'load {file}; {cmdRmd}; save {file}; quit'"))
        })
      writeLines(exclude_nonStd, file.path("protein_pdb", ".pymol_cleaned"))
    }
    if (file.exists(file.path("protein_pdb", ".pymol_cleaned"))) {
      x <- methodAdd(x, "以 `pymol` {cite_show('LigandDockingSeelig2010')} 删除受体 PDB 文件中的非标准残基 (例如 {bind(exclude_nonStd)})。")
      x <- snapAdd(x, "以 `pymol` 去除非标准残基 ({bind(exclude_nonStd)})。")
    }
    dir.create(cpdir <- "protein_clean_pdb", FALSE)
    pdb.files <- lapply(pdb.files, 
      function(file) {
        newfile <- file.path(cpdir, basename(file))
        .pymol_select_polymer.protein(file, newfile)
        newfile
      })
    x <- snapAdd(x, "随后，以 `pymol` 仅保留蛋白结构 (polymer.protein)。")
    x <- methodAdd(x, "以 `pymol` 仅保留蛋白结构 (polymer.protein) (去除了原 PDB 中的配体等其他结构)。")
    x$res.receptor <- prepare_receptor(pdb.files)
    x <- methodAdd(x, "以 `ADFR` {cite_show('AutogridfrImpZhang2019')} 工具组的准备受体蛋白的 PDBQT 文件 (以 `prepare_receptor` 添加氢原子，并转化为 PDBQT 文件) 。请参考 <https://autodock-vina.readthedocs.io/en/latest/docking_basic.html>。")
    x <- snapAdd(x, "以 `ADFR` 工具给受体添加氢原子，转化为 PDBQT 文件。")
    names <- names(x$res.receptor)
    gotSymbols <- names(x$used_pdbs)[ match(tolower(names(x$res.receptor)), tolower(x$used_pdbs)) ]
    x$res.receptor.symbol <- gotSymbols
    fun <- function(x) x[ !x %in% gotSymbols ]
    message("Not got: ", paste0(fun(unique(object(x)$hgnc_symbol)), collapse = ", "))
    if (filter) {
      if (!is.null(x$.layout)) {
        message("Customize using `x$.layout` columns: ", paste0(colnames(x$.layout)[1:2], collapse = ", "))
        x <- filter(x, x$.layout[[ 1 ]], x$.layout[[ 2 ]])
      }
    }
    if (any(lengths(x$dock_layout) > 1)) {
      names <- rep(names(x$dock_layout), lengths(x$dock_layout))
      dock_layout <- unlist(x$dock_layout, use.names = FALSE)
      names(dock_layout) <- names
      x$dock_layout <- dock_layout
    }
    if (use_complex && any(isThats <- grpl(names(used_pdbs), "\\+"))) {
      isUsedChains <- vapply(x$dock_layout, 
        function(x) {
          grpl(x, "_")
        }, logical(1))
      cpdsUsedChains <- names(x$dock_layout)[isUsedChains]
      dat <- dplyr::distinct(
        data.frame(
          cpd = cpdsUsedChains, pdb = gs(
            unlist(x$dock_layout[isUsedChains]), "_[a-z]$", ""
          )
        )
      )
      x$dock_layout <- c(x$dock_layout, nl(dat$cpd, dat$pdb))
    }
    if (!length(x$res.receptor) || !length(x$res.ligand)) {
      stop(
        '!length(x$res.receptor) || !length(x$res.ligand), no any ligand or receptor?'
      )
    }
    return(x)
  })

.get_pdb_files <- function(x, cl = 10, pattern = NULL, 
  extra_pdb.files = NULL, extra_layouts = NULL,
  tryGetFromAlphaFold = TRUE)
{
  if (!is(x, "job")) {
    stop('!is(x, "job").')
  }
  if (is.null(x$dock_layout)) {
    stop('is.null(x$dock_layout).')
  }
  ids <- rm.no(unlist(x$dock_layout, use.names = FALSE))
  if (length(ids)) {
    # The returned files is got by list.files
    pdb.files <- get_pdb(ids, cl = cl, mkdir.pdb = "protein_pdb")
    pdb.files <- pdb.files[ names(pdb.files) %in% tolower(x$used_pdbs) ]
    x <- methodAdd(x, "以 RCSB API  (<https://www.rcsb.org/docs/programmatic-access/web-apis-overview>) 获取蛋白 PDB 文件。")
    x <- snapAdd(x, "从 RCSB PDB 获取 PDB 文件。")
  } else {
    pdb.files <- NULL
  }
  if (length(x$pdb_notGot) && tryGetFromAlphaFold) {
    res_af <- get_pdb_from_alphaFold(x$pdb_notGot, "protein_pdb")
    x <- methodAdd(x, "以 R 包 `UniProt.ws` ({packageVersion('UniProt.ws')}) 获取基因 (symbol) 的 `UniProtKB-Swiss-Prot` ID (Entry ID)，随后，以 Entry ID 从数据库 `AlphaFold` (<https://alphafold.ebi.ac.uk/>) 获取数据库 `PDB` 中不包含的蛋白结构 (预测的结构)。")
    x <- snapAdd(x, "对于未从 `PDB` 数据库找到结构文件的，从数据库 `AlphaFold` 获取 {less(x$pdb_notGot)} 预测的蛋白结构 (根据 `UniProtKB-Swiss-Prot` ID，详见方法章节)。")
    x$pdb_notGot_uniprot <- res_af$info
    # extra_pdb.files: hgnc_symbol = file
    if (!is.null(extra_pdb.files)) {
      # symbol = pdb_ID
      customInput <- basename(tools::file_path_sans_ext(extra_pdb.files))
      customFiles <- extra_pdb.files
      # pdb_ID = files
      names(customFiles) <- unname(customInput)
    } else {
      customInput <- NULL
      customFiles <- NULL
    }
    extra_pdb.files <- c(
      res_af$files, customFiles
    )
    used_pdbs_extra <- c(res_af$used_pdbs, customInput)
    # need revision: x$dock_layout, pdb.files, x$used_pdbs
  }
  if (!is.null(extra_pdb.files)) {
    pdb.files <- c(pdb.files, extra_pdb.files)
    x$used_pdbs <- c(x$used_pdbs, used_pdbs_extra)
    layoutNeedRevise <- TRUE
  } else {
    layoutNeedRevise <- FALSE
  }
  if (!is.null(extra_layouts)) {
    if (is(extra_layouts, "character")) {
      extra_layouts <- as.list(extra_layouts)
    }
    x$dock_layout <- c(x$dock_layout, extra_layouts)
  } else if (layoutNeedRevise) {
    x$dock_layout %<>% lapply(function(x) c(x, unname(used_pdbs_extra)))
  }
  if (!is.null(x$pdb_MultiChains)) {
    pdb_MultiChains <- tolower(unique(x$pdb_MultiChains))
    fun_extract_cpdName <- function(lines) {
      chains <- stringr::str_extract(lines, "(?<=CHAIN: ).*(?=;)")
      names <- stringr::str_extract(lines, "(?<=MOLECULE: ).*(?=;)")
      names <- names[ !is.na(names) ]
      names <- if (length(names)) names else "__Omit__"
      names(names) <- chains[ !is.na(chains) ]
      return(names)
    }
    used_pdbs <- x$used_pdbs
    pdb_files_MultiChains <- lapply(pdb_MultiChains,
      function(pdb) {
        file <- pdb.files[[ pdb ]]
        contentPdb <- readLines(file)
        contentPdb <- contentPdb[ grpl(contentPdb, "^COMPND") ]
        contentPdb <- sep_list(contentPdb, "^COMPND.*MOL_ID", 0)
        if (length(contentPdb) > 1) {
          names <- unlist(lapply(unname(contentPdb), fun_extract_cpdName))
          message(glue::glue("PDB: {pdb.files[[ pdb ]]} has multiple compounds: {bind(names)}."))
          if (any(nchar(names(names)) > 1)) {
            message(
              glue::glue('Chain names not in expection: {bind(names(names), co = " | ")}')
            )
            return(file)
          }
          chains <- names(names)
          newfiles <- add_filename_suffix(file, tolower(chains))
          .pymol_select_chains(file, chains, newfiles)
          gene_pdbChains <- nl(names, tools::file_path_sans_ext(basename(newfiles)), FALSE)
          if (any(hasWhich <- names(used_pdbs) %in% names(gene_pdbChains))) {
            genes <- names(used_pdbs)[hasWhich]
            message(
              glue::glue("Detected input gene ({bind(genes)}) in input chains: bind(gene_pdbChains).")
            )
            used_pdbs[ hasWhich ] <- gene_pdbChains[ match(genes, names(gene_pdbChains)) ]
            complexPdb <- pdb
            names(complexPdb) <- paste0(
              names(gene_pdbChains), collapse = "+"
            )
            used_pdbs <<- c(used_pdbs, complexPdb)
          }
          c(file, newfiles)
        } else {
          file
        }
      })
    lessUsedChain <- head(
      used_pdbs[ grpl(used_pdbs, "_[a-e]$") ], n = 2
    )
    if (length(lessUsedChain)) {
      x <- snapAdd(x, "对于复合体 PDB (文件中包含支链分子信息)，将使用对应的支链 (以 `pymol` 获取支链) 进行分子对接 (例如 {bind(names(lessUsedChain))}，使用 {bind(lessUsedChain)}) 。")
    }
    x$used_pdbs <- used_pdbs
    x$dock_layout <- lapply(x$dock_layout,
      function(vec) {
        vec[ names(vec) %in% names(used_pdbs) ] <- used_pdbs[ match(
          names(vec), names(used_pdbs)
          ) ]
        vec
      })
    pdb_files_MultiChains <- unlist(
      pdb_files_MultiChains, use.names = FALSE
    )
    names(pdb_files_MultiChains) <- tools::file_path_sans_ext(basename(pdb_files_MultiChains))
    x$pdb_files_MultiChains <- pdb_files_MultiChains
    pdb.files <- c(pdb.files[!names(pdb.files) %in% pdb_MultiChains], pdb_files_MultiChains)
  }
  if (!is.null(pattern)) {
    pdb.files <- filter_pdbs(pdb.files, pattern)
  }
  namel(x, pdb.files)
}

get_pdb_from_alphaFold <- function(symbols, dir = "protein_pdb")
{
  message("Try get from alphaFold database (taxId: 9606, Human)")
  info <- UniProt.ws::mapUniProt(
    from = "Gene_Name",
    to = "UniProtKB-Swiss-Prot",
    columns = c("accession", "id"),
    query = list(taxId = 9606, ids = symbols)
  )
  if (any(which <- !symbols %in% info$From)) {
    message(glue::glue("Not Got in alphaFold: {bind(symbols[which])}"))
  }
  if (any(duplicated(info$From))) {
    message('any(duplicated(info$From)), distinct herein.')
    info <- dplyr::distinct(info, From, .keep_all = TRUE)
  }
  files <- apply(info, 1, 
    function(vec) {
      id <- vec[[ "Entry" ]]
      if (!is.na(id)) {
        url <- glue::glue("https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v4.pdb")
        save <- file.path(dir, paste0(id, ".pdb"))
        if (!file.exists(save)) {
          res <- try(download.file(url, save))
          if (inherits(res, "try-error")) {
            message(glue::glue("Download Failed of {id}, skip."))
            return(NULL)
          }
        }
        nl(vec[[ "Entry" ]], save, FALSE)
      } else {
        NULL
      }
    }, simplify = FALSE)
  files <- unlist(files)
  info <- dplyr::filter(info, Entry %in% names(files))
  lst <- list(
    files = files, 
    info = info, used_pdbs = nl(info$From, info$Entry, FALSE)
  )
  lst
}

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
    if (is.null(x$used_pdbs)) {
      stop('is.null(x$used_pdbs).')
    }
    layout <- pdb <- x$used_pdbs[ match(symbol, names(x$used_pdbs)) ]
    names(layout) <- unlist(cid)
    x$dock_layout <- as.list(layout)
    return(x)
  })

setMethod("step4", signature = c(x = "job_vina"),
  function(x, time = 3600 * 2, savedir = "vina_space", log = "/tmp/res.log", save.object = "vn3.rds")
  {
    step_message("Run vina ...")
    runs <- tibble::tibble(
      Ligand = rep(names(x$dock_layout), lengths(x$dock_layout)),
      Receptor = tolower(unlist(x$dock_layout, use.names = FALSE))
    )
    x$show_layout <- runs
    runs <- apply(runs, 1, unname, simplify = FALSE)
    x$runs <- runs
    x$savedir <- savedir
    n <- 0
    if (file.exists(log)) {
      file.remove("/tmp/res.log")
    }
    saveRDS(x, save.object)
    message(glue::glue("Save this 'job_vina' as: {save.object}."))
    x <- methodAdd(x, "以 `AutoDock-Vina` 提供的工具 (`prepare_gpf.py`) (<https://github.com/ccsb-scripps/AutoDock-Vina>) 创建 GPF (grid parameter file)。")
    x <- methodAdd(x, "以 `ADFR` {cite_show('AutogridfrImpZhang2019')} 工具 `autogrid4` 计算亲和图谱 (Affinity Maps)。")
    x <- snapAdd(x, "以 `ADFR` 创建 Affinity Maps (详见方法章节) 。")
    if (any(grpl(names(x$res.receptor), "[A-Z]"))) {
      message('any(grpl(names(x$res.receptor), "[A-Z]")), convert to lower case.')
      names(x$res.receptor) <- tolower(names(x$res.receptor) )
    }
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
    x <- methodAdd(
      x, "运行 AutoDock-Vina {cite_show('AutodockVina1Eberha2021')} (parameters: scoring = ad4; exhaustiveness = 32)。"
    )
    x <- snapAdd(x, "以 `Autodock-Vina` 进行自动分子对接。")
    if (.Platform$OS.type == "unix" && Sys.which("notify-send") != "") {
      system("notify-send 'AutoDock vina' 'All job complete'")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_vina"),
  function(x, compounds, by.y, facet = "Ingredient_name", excludes = NULL, top = NULL, cutoff.af = NULL)
  {
    step_message("Summary and visualization for results.")
    x$summary_vina <- summary_vina(x$savedir)
    x$summary_vina <- dplyr::filter(x$summary_vina,
      PubChem_id %in% !!names(x$dock_layout),
      PDB_ID %in% tolower(unlist(x$dock_layout, use.names = FALSE))
    )
    if (!is.null(top)) {
      x$summary_vina <- split_lapply_rbind(x$summary_vina, ~ PDB_ID,
        function(x) {
          dplyr::slice_min(x, Affinity, n = top)
        }
      )
    }
    res_dock <- dplyr::mutate(
      x$summary_vina, PubChem_id = as.integer(PubChem_id),
      hgnc_symbol = names(x$used_pdbs)[ match(PDB_ID, tolower(x$used_pdbs)) ]
    )
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
    res_dock <- .set_lab(res_dock, sig(x), "All combining Affinity data")
    res_dock <- setLegend(res_dock, "分子对接得分 (亲和度) 附表。")
    data <- dplyr::distinct(res_dock, PubChem_id, hgnc_symbol, .keep_all = TRUE)
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
      geom_col(aes(x = reorder(hgnc_symbol, Affinity, decreasing = TRUE), y = Affinity, fill = Affinity), width = .7) +
      geom_text(data = dplyr::filter(data, Affinity <= 0),
        aes(x = hgnc_symbol, y = Affinity - .5, label = round(Affinity, 1)), hjust = 1) +
      labs(x = "", y = "Affinity (kcal/mol)") +
      coord_flip() +
      ylim(zoRange(c(-1, data$Affinity, 1), 1.4)) +
      facet_wrap(as.formula(paste0("~ Hmisc::capitalize(paste0(", facet, ", \" (CID: \", PubChem_id, \")\"))")),
        ncol = 1, scales = "free_y") +
      theme()
    height <- nrow(data) + 1
    p.res_vina <- wrap(p.res_vina, 8, if (height > 9.5) 9.5 else height)
    p.res_vina <- .set_lab(p.res_vina, sig(x), "Overall combining Affinity")
    p.res_vina <- setLegend(p.res_vina, "为分子对接亲和度得分可视化，能量越低，代表亲和度越高。")
    x@tables[[ 5 ]] <- namel(res_dock, unique_tops = data)
    x@plots[[ 5 ]] <- namel(p.res_vina)
    return(x)
  })

setMethod("step6", signature = c(x = "job_vina"),
  function(x, time = 30, top = NULL, save = TRUE, unique = FALSE, symbol = NULL, cpd = NULL)
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
      data <- dplyr::distinct(data, hgnc_symbol, .keep_all = TRUE)
    }
    if (!is.null(top)) {
      data <- head(data, n = top)
    }
    figs <- pbapply::pbapply(data, 1, simplify = FALSE,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], timeLimit = time, save = save)
      }
    )
    figs <- unlist(figs, recursive = FALSE)
    lab(figs) <- paste(sig(x), "docking visualization")
    names(figs) <- paste0("Top", seq_along(figs), "_", names(figs))
    x <- snapAdd(x, "使用 `pymol` 将分子对接结果可视化。")
    figs <- setLegend(
      figs, glue::glue("分子对接结果。蛋白(Symbol: {data$hgnc_symbol}) (PDB: {data$PDB_ID}) 与化合物 (name: {data$Ingredient_name}) (PubChem CID: {data$PubChem_id})，亲和度为 {data$Affinity}。")
    )
    x@plots[[ 6 ]] <- figs
    data <- .set_lab(data, sig(x), "Metadata of visualized Docking")
    data <- dplyr::select(data, -dir, -file)
    x@tables[[ 6 ]] <- namel(data)
    return(x)
  })

setMethod("step7", signature = c(x = "job_vina"),
  function(x, save = TRUE){
    step_message("Show docking results in deep and in detail.")
    data <- x@tables$step6$data
    figs <- pbapply::pbapply(data, 1, simplify = FALSE,
      function(v) {
        vinaShow(v[[ "Combn" ]], v[[ "PDB_ID" ]], save = save, detail = TRUE)
      }
    )
    figs <- unlist(figs, recursive = FALSE)
    lab(figs) <- paste(sig(x), "docking interaction details")
    names(figs) <- paste0("Top", seq_along(figs), "_", names(figs))
    figs <- setLegend(
      figs, glue::glue("分子对接结果局部细节。图中展示了蛋白(Symbol: {data$hgnc_symbol}) (PDB: {data$PDB_ID}) 与化合物 (name: {data$Ingredient_name}) (PubChem CID: {data$PubChem_id}) 之间的氢键结合。")
    )
    x@plots[[ 7 ]] <- figs
    return(x)
  })

pretty_docking <- function(protein, ligand, path,
  save = "annotation.png",
  script = paste0(.expath, "/pretty_docking.pymol"))
{
  script <- readLines(script)
  script <- gs(script, "{{protein.pdbqt}}", protein, fixed = TRUE)
  script <- gs(script, "{{ligand.pdbqt}}", ligand, fixed = TRUE)
  temp <- tempfile("Pymol", fileext = ".pml")
  writeLines(script, temp)
  cli::cli_alert_info(paste0("Pymol run script: ", temp))
  output <- paste0(path, "/", save)
  message("Save png to ", output)
  expr <- paste0(" png ", save, ",2500,2000,dpi=300")
  gett(expr)
  cdRun(pg("pymol"), " ",
    " -d \"run ", temp, "\"",
    " -d \"ray; ", expr, "\"",
    path = path)
  return(output)
}

vina_limit <- function(lig, recep, timeLimit = 120, dir = "vina_space", ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit, dir = dir), TRUE)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoring = "ad4", stout = "/tmp/res.log", timeLimit = 60,
  excludes.atom = c("G0", "CG0"), x)
{
  remote <- FALSE
  if (!missing(x)) {
    if (is.remote(x)) {
      remote <- TRUE
    }
  }
  if (!file.exists(dir)) {
    dir.create(dir, FALSE)
  } 
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
  if (!file.exists(paste0(wd, "/", subdir, "_out.pdbqt"))) {
    dir.create(wd, FALSE)
    file.copy(c(recep, lig), wd)
    .message_info("Generating affinity maps", subdir)
    .cdRun <- function(...) cdRun(..., path = wd)
    files <- basename(c(lig, recep))
    .cdRun(pg("prepare_gpf.py"), " -l ", files[1], " -r ", files[2], " -y")
    if (!is.null(excludes.atom)) {
      message("Excludes atom type: ", paste0(excludes.atom, collapse = ", "))
      ## ligand type and map file
      .cdRun("sed -i",
        " -e '/", paste0(paste0("^map.*", excludes.atom), collapse = "\\|"), "/d'",
        " -e 's/", paste0(paste0("\\b", excludes.atom, "\\b"), collapse = "\\|"), "//g'",
        " ", reals[2], ".gpf")
    }
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
      cat("\n$$$$\n", date(), "\n", subdir, "\n\n", file = stout, append = TRUE)
      try(.cdRun("timeout ", timeLimit, 
          " ", pg("vina"), "  --ligand ", files[1],
          " --maps ", reals[2],
          " --scoring ", scoring,
          " --exhaustiveness ", exhaustiveness,
          " --out ", subdir, "_out.pdbqt",
          " >> ", stout), TRUE)
    }
  }
}

vinaShow <- function(Combn, recep, subdir = Combn, dir = "vina_space",
  timeLimit = 3, backup = NULL, save = TRUE, detail = FALSE)
{
  if (!file.exists(path <- paste0(dir, "/", subdir))) {
    stop("file.exists(path <- paste0(dir, \"/\", subdir))")
  } 
  wd <- paste0(dir, "/", subdir)
  out <- paste0(Combn, "_out.pdbqt")
  recep <- paste0(recep, ".pdbqt")
  if (!file.exists(recep)) {
    maybethat <- list.files(wd, recep, ignore.case = TRUE)
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
    expr <- paste0(" png ", res, ",2500,2000,dpi=300")
    gett(expr)
    if (!save) {
      expr <- ""
    }
    try(.cdRun("timeout ", timeLimit, 
        " ", pg("pymol"), " ",
        " -d \"load ", out, ";",
        " load ", recep, ";",
        " ray;", expr, "\" "), TRUE)
  }
  if (is.character(backup)) {
    dir.create(backup, FALSE)
    file.copy(img, backup, TRUE)
  }
  fig <- file_fig(Combn, img)
  lab(fig[[1]]) <- paste0("docking ", Combn, if (!detail) NULL else " detail")
  return(fig)
}

summary_vina <- function(space = "vina_space", pattern = "_out\\.pdbqt$")
{
  files <- list.files(space, pattern, recursive = TRUE, full.names = TRUE)
  res_dock <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (length(lines) >= 1) {
        name <- gsub("_out\\.pdbqt", "", basename(file))
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
    res_dock, PubChem_id = gs(Combn, "(.*)_into_.*", "\\1"),
    PDB_ID = tolower(gs(Combn, ".*_into_(.*)", "\\1")),
    dir = paste0(space, "/", Combn),
    file = paste0(dir, "/", Combn, "_out.pdbqt")
  )
  dplyr::select(res_dock, PubChem_id, PDB_ID, Affinity, dir, file, Combn)
}

smiles_as_sdfs.obabel <- function(smiles) {
  lst.sdf <- pbapply::pbsapply(smiles, simplify = FALSE,
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
  lst[ vapply(lst, function(v) if (length(v) > len) TRUE else FALSE, logical(1)) ]
}

ld_cutRead <- function(file, cols, abnum = TRUE, sep = "\t", tmp = "/tmp/ldtmp.txt") {
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

mk_prepare_ligand.sdf <- function(sdf_file, mkdir.pdbqt = "pdbqt", check = FALSE) {
  dir.create(mkdir.pdbqt, FALSE)
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
  lst$pdbqt <- list.files(mkdir.pdbqt, "\\.pdbqt$", full.names = TRUE)
  lst$pdbqt.num <- length(lst$pdbqt)
  lst$pdbqt.cid <- stringr::str_extract(lst$pdbqt, "(?<=/|^)[0-9]{1,}")
  return(lst)
}

select_files_by_grep <- function(files, pattern){
  vapply(files,
    function(file) {
      any(grepl(pattern, readLines(file), TRUE))
    }, logical(1))
}

filter_pdbs <- function(files, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS") {
  files[ select_files_by_grep(files, pattern) ]
}

prepare_receptor <- function(files, mkdir.pdbqt = "protein_pdbqt") {
  dir.create(mkdir.pdbqt, FALSE)
  file <- lapply(files,
    function(file) {
      if (!is.null(file)) {
        newfile <- gsub("\\.pdb$", ".pdbqt", basename(file))
        newfile <- paste0(mkdir.pdbqt, "/", newfile)
        system(paste0("prepare_receptor -r ", file, " -o ", newfile, " -A hydrogens"))
        return(newfile)
      }
    })
  file <- file[ !vapply(file, is.null, logical(1)) ]
  file[ vapply(file, file.exists, logical(1)) ]
}

setMethod("set_remote", signature = c(x = "job_vina"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })

cal_3d_sdf <- function(sdf) {
  output <- gs(sdf, "\\.sdf$", "_3D.sdf")
  isThat <- TRUE
  if (file.exists(output)) {
    isThat <- usethis::ui_yeah("Overwrite the exists file?")
  }
  if (isThat) {
    cdRun(pg("obgen"), " ", sdf, " -ff UFF > ", output)
  }
  output
}

.pymol_select_polymer.protein <- function(file, newfile) {
  cmd <- glue::glue(
    "select protein, polymer.protein; save {newfile}, protein"
  )
  cdRun(glue::glue("{pg('pymol')} -c -Q -d 'load {file}; {cmd}; quit'"))
}

.pymol_select_chains <- function(file, chains, newfiles = add_filename_suffix(file, chains))
{
  cmd <- glue::glue(
    "select chain{chains}, chain {chains}; save {newfiles}, chain{chains}"
  )
  cmd <- paste0(cmd, collapse = "; ")
  cdRun(glue::glue("{pg('pymol')} -c -Q -d 'load {file}; {cmd}; quit'"))
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

# deprecated:
# if (bdb) {
#   if (file.exists(bdb_file)) {
#     message("Database `BindingDB` used for pre-filtering of docking candidates.")
#     bdb <- ld_cutRead(bdb_file, c("PubChem CID", "PDB ID(s) of Target Chain"))
#     bdb <- dplyr::filter(bdb, `PubChem CID` %in% object(x)$cids)
#     colnames(bdb) <- c("PubChem_id", "pdb_id")
#     x@params$bdb_compounds_targets <- bdb
#     bdb <- nl(bdb$PubChem_id, strsplit(bdb$pdb_id, ","))
#     bdb <- lst_clear0(bdb)
#     bdb <- lapply(bdb, function(v) v[ v %in% x$targets_annotation$pdb ])
#     bdb <- lst_clear0(bdb)
#     if (!length(bdb)) {
#       stop("No candidates found in BindingDB.")
#     }
#     x$dock_layout <- bdb
#     .add_internal_job(.job(method = "Database `BindingDB` used for pre-filtering of docking candidates.",
#         cite = "[@BindingdbIn20Gilson2016]"))
#   }
# }
