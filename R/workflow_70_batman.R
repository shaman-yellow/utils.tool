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
  return(x)
}

setMethod("step0", signature = c(x = "job_batman"),
  function(x){
    step_message("Prepare your data with function `job_batman`.")
  })

setMethod("step1", signature = c(x = "job_batman"),
  function(x, anno = F, filter.hob = T, filter.dl = T, test = F)
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
    if (anno) {
      x <- anno(x, "LiteratureCount")
      x$compounds_info <- dplyr::arrange(x$compounds_info, dplyr::desc(LiteratureCount))
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
    lst <- filter_hob_dl(unique(herbs_compounds$CID), filter.hob, filter.dl, test)
    if (length(lst)) {
      p.upset <- .set_lab(lst$p.upset, sig(x), "Prediction of HOB and Drug-Likeness")
      x@plots[[ 1 ]] <- namel(p.upset)
      message("Filter the `herbs_compounds` and `compounds_targets`.")
      cids <- lst$ids
      herbs_compounds <- dplyr::filter(herbs_compounds, CID %in% !!cids)
      compounds_targets <- dplyr::filter(compounds_targets, PubChem_CID %in% !!cids)
      predicted <- dplyr::filter(predicted, PubChem_CID %in% !!cids)
    }
    x@tables[[ 1 ]] <- namel(herbs_compounds, compounds_targets, predicted)
    return(x)
  })

setMethod("step2", signature = c(x = "job_batman"),
  function(x, use.predicted = T, cutoff = .9)
  {
    step_message("Format target genes symbol")
    compounds_targets <- x@tables$step1$compounds_targets
    compounds_targets <- reframe_split(compounds_targets, "known_target_proteins", "\\|")
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
        PubChem_CID, IUPAC_name, symbol, .keep_all = T)
      compounds_targets <- dplyr::relocate(compounds_targets, IUPAC_name, symbol)
    }
    compounds_targets <- .set_lab(compounds_targets, sig(x), "Compounds and targets")
    x@tables[[ 2 ]] <- namel(compounds_targets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_batman"),
  function(x, HLs = NULL, disease = NULL, shortName = T, test = F)
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

get_batman_data <- function(savedir = .prefix("batman", "db"), reload = F) {
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
      querys <- grouping_vec2list(x@params$compounds_info$cids, 100, T)
      anno <- pbapply::pblapply(querys,
        function(query) {
          PubChemR::get_properties(
            properties = "LiteratureCount",
            identifier = query,
            namespace = "cid"
          )
        })
      anno <- unlist(anno, recursive = F)
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

.add_properties.PubChemR <- function(name = c(literature_count = "LiteratureCount"))
{
  maps <- PubChemR:::property_map 
  if (any(names(name) %in% names(maps)) || any(name %in% unlist(maps))) {
    message("The properties already in `PubChemR:::property_map`")
  } else {
    maps <- c(maps, as.list(name))
    replaceFunInPackage("property_map", maps, "PubChemR")
  }
}

try_get_syn <- function(cids) {
  syns <- e(PubChemR::get_synonyms(unique(cids)))
  syns <- dplyr::bind_rows(syns)
  syns <- dplyr::rename(syns, syno = Synonym)
  syns <- dplyr::filter(syns, !!!.filter_pick.general)
  syns <- dplyr::group_by(syns, CID)
  syns <- dplyr::reframe(syns, Synonym = PickGeneral(syno))
  dplyr::ungroup(syns)
}

get_smiles_batch <- function(cids, n = 100, sleep = .5) {
  res <- .get_properties_batch(cids, n = 100, as_dataframe = T, properties = "IsomericSMILES", sleep = sleep)
  frbind(res, fill = T)
}

.get_properties_batch <- function(ids, ..., n = 100, sleep = .5) {
  groups <- grouping_vec2list(unique(ids), n, T)
  pbapply::pblapply(groups,
    function(ids) {
      Sys.sleep(sleep)
      PubChemR::get_properties(identifier = ids, ...)
    })
}
