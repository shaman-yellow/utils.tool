# ==========================================================================
# workflow of batman
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_batman <- setClass("job_batman", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "batman",
    info = c("http://bionet.ncpsb.org.cn/batman-tcm/#/home"),
    cite = "[@BatmanTcm20Kong2024]",
    method = "Database `BATMAN-TCM` was used as source data of TCM ingredients and target proteins",
    tag = "herb:BATMAN, pharm",
    analysis = "BATMAN 网络药理学"
    ))

job_batman <- function(herbs, db = get_batman_data())
{
  type <- c(herb_browse = "herbs_ingredients",
    known_browse_by_ingredients = "targets",
    predicted_browse_by_ingredients = "predicted_targets"
  )
  names(db) <- dplyr::recode(names(db), !!!as.list(type))
  herbs_info <- dplyr::filter(db$herbs_ingredients, Chinese.Name %in% !!herbs)
  print(herbs_info)
  message("Got the herbs:\n\n\t", paste0(herbs_info$Chinese.Name, collapse = ", "),
    "\n\n", "\tTotal: ", colSum(herbs_info$Chinese.Name))
  x <- .job_batman(object = db)
  x$herbs_info <- herbs_info
  x$herbs <- herbs
  return(x)
}

setMethod("step0", signature = c(x = "job_batman"),
  function(x){
    step_message("Prepare your data with function `job_batman`.")
  })

setMethod("step1", signature = c(x = "job_batman"),
  function(x, anno = FALSE, filter.hob = TRUE, filter.dl = FALSE, test = FALSE)
  {
    step_message("Filter and format the compounds and targets.")
    x$herbs_info <- dplyr::mutate(x$herbs_info,
      Ingredients = strsplit(Ingredients, "\\|"),
      CID = lapply(Ingredients,
        function(x) as.integer(strx(x, "(?<=\\()[0-9]+(?=\\)$)")))
    )
    x$compounds_info <- tibble::tibble(
      cids = unlist(x$herbs_info$CID),
      name = gs(unlist(x$herbs_info$Ingredients), "\\([0-9]+\\)$", "")
    )
    x$compounds_info <- dplyr::distinct(x$compounds_info)
    x <- methodAdd(x, "从数据库 `BATMAN-TCM` ({cite_show('BatmanTcm20Kong2024')}) 中获取 {bind(x$herbs_info$Pinyin.Name)} 等中药的成分、靶点数据。(即中药：{bind(x$herbs)})。")
    x <- snapAdd(x, "从数据库 `BATMAN-TCM` 中药的成分、靶点数据 (详见方法章节) 。")
    if (anno) {
      x <- anno(x, "LiteratureCount")
      x$compounds_info <- dplyr::arrange(x$compounds_info, dplyr::desc(LiteratureCount))
      x <- methodAdd(x, "以 `PubChemR` 获取 `PubChemR` 数据库中，化合物的文献记录数量。", TRUE)
    }
    cids <- rm.no(unlist(x$herbs_info$CID))
    compounds_targets <- dplyr::filter(
      object(x)$targets, PubChem_CID %in% !!cids
    )
    predicted <- dplyr::filter(
      object(x)$predicted_targets, PubChem_CID %in% !!cids
    )
    object(x) <- NULL
    herbs_compounds <- reframe_col(
      dplyr::select(x$herbs_info, -Ingredients), "CID", unlist
    )
    herbs_compounds <- dplyr::distinct(herbs_compounds, Pinyin.Name, Latin.Name, CID)
    lst <- filter_hob_dl(unique(herbs_compounds$CID), filter.hob, filter.dl, test, symbol = sig(x))
    if (filter.hob) {
      x <- snapAdd(x, "通过 `PubChemR` 获取化合物的结构式 (SMILES)。")
      x <- methodAdd(x, lst$hob)
      x <- snapAdd(x, lst$hob)
    }
    if (filter.dl) {
      x <- methodAdd(x, lst$dl)
      x <- snapAdd(x, lst$dl)
    }
    if (length(lst)) {
      p.upset <- .set_lab(lst$p.upset, sig(x), "Prediction of HOB and Drug-Likeness")
      p.upset <- setLegend(p.upset, "UpSet 图展示了经过 HOB 和 Drug-Likeness 预测前后的彼此交集。")
      x@plots[[ 1 ]] <- namel(p.upset)
      message("Filter the `herbs_compounds` and `compounds_targets`.")
      cids <- lst$ids
      herbs_compounds <- dplyr::filter(herbs_compounds, CID %in% !!cids)
      compounds_targets <- dplyr::filter(compounds_targets, PubChem_CID %in% !!cids)
      predicted <- dplyr::filter(predicted, PubChem_CID %in% !!cids)
      x <- snapAdd(x, "经过筛选，")
    }
    x@tables[[ 1 ]] <- namel(herbs_compounds, compounds_targets, predicted)
    s.com <- try_snap(herbs_compounds, "Pinyin.Name", "CID")
    s.n <- length(unique(compounds_targets$PubChem_CID))
    x <- snapAdd(x, "各中药的化合物组成统计 (可能有交叉涵盖)：{s.com}。
      共 {length(unique(herbs_compounds$CID))} 个化合物 (注：根据唯一 PubChem CID 统计)，
      其中，含有靶点信息记录的化合物共 {s.n} 个 (非重复)。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_batman"),
  function(x, use.predicted = TRUE, cutoff = .9)
  {
    step_message("Format target genes symbol")
    compounds_targets <- x@tables$step1$compounds_targets
    compounds_targets <- reframe_split(compounds_targets, "known_target_proteins", "\\|")
    x <- methodAdd(x, "使用 `BATMAN-TCM` 数据库中的 `known_target_proteins` 作为成分靶点。", TRUE)
    if (use.predicted) {
      predicted <- x@tables$step1$predicted
      predicted <- reframe_split(predicted, "predicted_target_proteins", "\\|")
      predicted <- dplyr::mutate(predicted,
        entrez_id = as.integer(strx(predicted_target_proteins, "[0-9]+")),
        score = as.double(strx(predicted_target_proteins, "(?<=\\()[0-9.]*(?=\\)$)"))
      )
      predicted <- dplyr::filter(predicted, score > cutoff)
      predicted <- dplyr::select(predicted, -predicted_target_proteins, -score)
      entrez_ids <- rm.no(predicted$entrez_id)
      mart <- new_biomart("hsa")
      anno <- filter_biomart(mart,
        c("entrezgene_id", "hgnc_symbol"), "entrezgene_id", entrez_ids
      )
      predicted <- map(predicted, "entrez_id",
        anno, "entrezgene_id", "hgnc_symbol", col = "hgnc_symbol")
      lst <- list(
        known = dplyr::select(compounds_targets, PubChem_CID, IUPAC_name, symbol = known_target_proteins),
        predicted = dplyr::select(predicted, PubChem_CID, IUPAC_name, symbol = hgnc_symbol)
      )
      compounds_targets <- frbind(lst, idcol = "target_type")
      compounds_targets <- dplyr::distinct(compounds_targets,
        PubChem_CID, IUPAC_name, symbol, .keep_all = TRUE)
      compounds_targets <- dplyr::relocate(compounds_targets, IUPAC_name, symbol)
      x <- methodAdd(x, "此外，还使用了 `BATMAN-TCM` 数据库中的 `predicted_target_proteins` 作为成分靶点，并设定 分数 cut-off 为 {cutoff}。合并靶点数据。", TRUE)
      x <- methodAdd(x, "以 `BiomaRt` ({cite_show('MappingIdentifDurinc2009')}) 对靶点信息的 entrez_id 转化为基因 Symbol (hgnc_symbol) 。", TRUE)
    }
    compounds_targets <- .set_lab(compounds_targets, sig(x), "Compounds and targets")
    x@tables[[ 2 ]] <- namel(compounds_targets)
    x <- snapAdd(x, "共包含靶点 {length(unique(compounds_targets$symbol))} 个 (非重复)。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_batman"),
  function(x, HLs = NULL, disease = NULL, shortName = TRUE, test = FALSE)
  {
    step_message("Network pharmacology.")
    ########################
    ########################
    cli::cli_alert_info("Run: job_herb")
    hb <- .job_herb(step = 2L)
    herbs_compounds <- map(x@tables$step1$herbs_compounds, "CID",
      x$compounds_info, "cids", "name", col = "Ingredient.name"
    )
    if (shortName) {
      message("Get Synonyms (short name) for compounds.")
      if (test) {
        herbs_compounds <- head(herbs_compounds, 300)
      }
      syns <- try_get_syn(unique(herbs_compounds$CID))
      herbs_compounds <- map(herbs_compounds, "CID", syns, "CID", "Synonym", col = "Synonym")
      herbs_compounds <- dplyr::mutate(herbs_compounds,
        Ingredient.name = ifelse(is.na(Synonym), Ingredient.name, Synonym))
      if (test) {
        stop_debug(namel(herbs_compounds))
      }
      x <- methodAdd(x, "以 `PubChemR` 获取化合物同义名 (Synonym)，按正则表达式 (Regex) 匹配化合物简短的同义名用以化合物注释。", TRUE)
    }
    hb@tables$step1$herbs_compounds <- dplyr::select(
      herbs_compounds, herb_id = Pinyin.Name, Ingredient.id = CID,
      Ingredient.name
    )
    hb@tables$step2$compounds_targets <- dplyr::select(
      x@tables$step2$compounds_targets,
      Ingredient_id = PubChem_CID, Target.name = symbol
    )
    hb@params$herbs_info <- dplyr::select(
      dplyr::mutate(x@params$herbs_info, Herb_ = Pinyin.Name),
      ## Herb_, it is, herb_id
      Herb_ = Pinyin.Name, Herb_pinyin_name = Pinyin.Name, Herb_cn_name = Chinese.Name
    )
    hb@object$herb <- hb@params$herbs_info
    hb <- step3(hb, disease = disease, HLs = HLs)
    x@plots[[ 3 ]] <- hb@plots$step3
    x$easyRead <- hb@params$easyRead
    x@tables[[ 3 ]] <- namel(disease_targets_annotation = hb@tables$step3$disease_targets_annotation)
    x$data.allu <- hb$data.allu
    return(x)
  })


reframe_split <- function(data, col, sep) {
  reframe_col(data, col, function(x) strsplit(x, sep)[[1]])
}

reframe_col <- function(data, col, fun) {
  groups <- colnames(data)[ colnames(data) != col ]
  data <- dplyr::group_by(data, !!!rlang::syms(groups))
  data <- dplyr::reframe(data, dplyr::across(!!rlang::sym(col), fun))
  dplyr::ungroup(data)
}

get_batman_data <- function(savedir = .prefix("batman", "db"), reload = FALSE) {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  rdata <- file.path(savedir, "all.rds")
  if (!file.exists(rdata) || reload) {
    url_base <- "http://batman2.cloudna.cn/downloadApiFile/data/browser/"
    files <- c(
      "herb_browse.txt",
      "known_browse_by_ingredients.txt.gz",
      "predicted_browse_by_ingredients.txt.gz"
    )
    lst <- pbapply::pblapply(files,
      function(file) {
        url <- paste0(url_base, file)
        target <- file.path(savedir, file)
        if (!file.exists(target)) {
          x <- RCurl::getURLContent(url)
          writeBin(as.vector(x), target)
        }
        data <- ftibble(target, sep = "\t")
        if (ncol(data) == 1 & file == "predicted_browse_by_ingredients.txt.gz") {
          colnames(data) <- "Name"
          data <- dplyr::mutate(data,
            PubChem_CID = as.integer(strx(Name, "^[0-9]+")),
            IUPAC_name = gs(Name, "^[^ ]* (.*) [^ ]*$", "\\1"),
            predicted_target_proteins = strx(Name, "[^ ]+$")
          )
          data <- dplyr::select(data, - Name)
        }
        data
      })
    names(lst) <- get_realname(files)
    saveRDS(lst, rdata)
    return(lst)
  } else {
    readRDS(rdata)
  }
}

setMethod("anno", signature = c(x = "job_batman"),
  function(x, use = "LiteratureCount")
  {
    message("Add annotation into `x@params$compounds_info`")
    use <- match.arg(use)
    if (use == "LiteratureCount") {
      .add_properties.PubChemR()
      querys <- grouping_vec2list(x@params$compounds_info$cids, 100, TRUE)
      anno <- pbapply::pblapply(querys,
        function(query) {
          props <- PubChemR::get_properties(
            properties = "LiteratureCount",
            identifier = query,
            namespace = "cid"
          )
          PubChemR::retrieve(props, .combine.all = TRUE)
        })
      anno <- unlist(anno, recursive = FALSE)
      anno <- dplyr::bind_rows(anno)
      anno <- dplyr::mutate_all(anno, as.integer)
      x@params$compounds_info <- map(
        x@params$compounds_info, "cids",
        anno, "CID", "LiteratureCount", col = "LiteratureCount"
      )
      x@params$compounds_info <- .set_lab(x@params$compounds_info, sig(x), "Compounds information")
    }
    return(x)
  })

setMethod("asjob_vina", signature = c(x = "job_batman"),
  function(x, ref, ...){
    ref <- resolve_feature_snapAdd_onExit("x", ref)
    herComTar <- x$easyRead
    # Target.name, symbol, hgnc_symbol or Others?
    herComTar <- dplyr::filter(herComTar, Target.name %in% ref)
    herComTar <- .set_lab(herComTar, sig(x), "herbs compounds and targets for docking")
    herComTar <- setLegend(herComTar, "含靶点基因的化合物和对应中药的附表 (Ingredient.id 为 PubChem CID)。")
    s.com <- try_snap(herComTar, "Herb_pinyin_name", "Ingredient.name")
    snapAdd_onExit("x", "含靶点基因的化合物与对应的中药统计：{s.com}。")
    comTar <- dplyr::distinct(
      herComTar, Ingredient.id, Ingredient.name, Target.name
    )
    herComTar <- .set_lab(herComTar, sig(x), "compounds and targets for docking")
    comTar <- setLegend(comTar, "将用于分子对接的化合物和靶点组合。")
    s.com <- try_snap(comTar, "Target.name", "Ingredient.name")
    snapAdd_onExit("x", "用于分子对接的靶点对应化合物统计：{s.com}。")
    object <- dplyr::distinct(
      # cpd names, hgnc_symbols, cids
      herComTar, Ingredient.name, Target.name, Ingredient.id
    )
    x <- job_vina(.layout = object)
    x$herComTar <- herComTar
    x$comTar <- comTar
    return(x)
  })

.add_properties.PubChemR <- function(name = c(literature_count = "LiteratureCount"))
{
  stop("As the PubchemR Version update, This function is deprecated.")
  maps <- PubChemR:::property_map 
  if (any(names(name) %in% names(maps)) || any(name %in% unlist(maps))) {
    message("The properties already in `PubChemR:::property_map`")
  } else {
    maps <- c(maps, as.list(name))
    replaceFunInPackage("property_map", maps, "PubChemR")
  }
}

try_get_syn <- function(cids, db_dir = .prefix("synonyms", "db")) {
  cids <- unique(cids)
  if (!is.null(db_dir)) {
    dir.create(db_dir, FALSE)
    db <- new_db(file.path(db_dir, "synonyms.rdata"), "CID")
    db <- not(db, cids)
    query <- db@query
    if (length(query)) {
      syns <- e(PubChemR::get_synonyms(query))
      syns <- e(PubChemR::synonyms(syns))
      db <- upd(db, syns)
    }
    syns <- dplyr::filter(db@db, CID %in% !!cids)
  }
  syns <- dplyr::rename(syns, syno = Synonyms)
  syns <- dplyr::filter(syns, !!!.filter_pick.general)
  syns <- dplyr::group_by(syns, CID)
  syns <- dplyr::reframe(syns, Synonym = PickGeneral(syno))
  dplyr::ungroup(syns)
}
