# ==========================================================================
# workflow of ctd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ctd <- setClass("job_ctd", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "ctd",
    info = c("https://ctdbase.org/downloads/"),
    cite = "[@ComparativeToxDavis2023]",
    method = "The Comparative Toxicogenomics Database (CTD) used for finding relationship between chemicals and disease",
    tag = "ctd",
    analysis = "CTD 化合物和疾病关联数据"
    ))

# https://github.com/r-dbi/RMariaDB

job_ctd <- function(db = get_ctd_data())
{
  .job_ctd(object = db)
}

setMethod("step0", signature = c(x = "job_ctd"),
  function(x){
    step_message("Prepare your data with function `job_ctd`.")
  })

setMethod("step1", signature = c(x = "job_ctd"),
  function(x, pattern)
  {
    step_message("Search disease in `CTD_diseases`")
    dis <- dplyr::filter(object(x)$CTD_diseases, grpl(DiseaseName, pattern, T))
    x$dis <- .set_lab(dis, sig(x), "matched disease")
    message("Get disease:\n\t", paste0(x$dis$DiseaseName, collapse = ",\n\t"))
    return(x)
  })

setMethod("step2", signature = c(x = "job_ctd"),
  function(x, rm.na = T){
    step_message("Get relationship table.")
    t.chemical <- dplyr::filter(
      object(x)$CTD_chemicals_diseases,
      DiseaseID %in% x@params$dis$DiseaseID
    )
    object(x)$CTD_chemicals_diseases <- NULL
    t.chemical <- dplyr::mutate(t.chemical,
      ChemicalID = paste0("MESH:", ChemicalID)
    )
    t.chemical <- map(t.chemical, "ChemicalID", object(x)$CTD_chemicals,
      "ChemicalID", "CID")
    t.chemical <- dplyr::filter(t.chemical, !is.na(CID))
    t.chemical <- .set_lab(t.chemical, sig(x), "Disease related Compounds")
    x@tables[[ 2 ]] <- namel(t.chemical)
    return(x)
  })

get_ctd_data <- function(name = "chemical.rds",
  savedir = .prefix("ctd", "db"), reload = F)
{
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  name <- match.arg(name)
  rdata <- paste0(savedir, "/", name)
  if (!file.exists(rdata) || reload) {
    url_base <- "https://ctdbase.org/reports/"
    if (name == "chemical.rds") {
      files <- c(
        "CTD_chemicals_diseases.tsv.gz",
        "CTD_chemicals.tsv.gz",
        "CTD_diseases.tsv.gz"
      )
      selects <- list(NULL, 1:2, 1:2)
      col.names <- list(
        c("ChemicalName", "ChemicalID", "CasRN",
          "DiseaseName", "DiseaseID", "DirectEvidence", "InferenceGeneSymbol",
          "InferenceScore", "OmimIDs", "PubMedIDs"),
        c("ChemicalName", "ChemicalID"),
        c("DiseaseName", "DiseaseID")
      )
    }
    lst <- pbapply::pbmapply(files, selects, col.names, SIMPLIFY = F,
      FUN = function(file, select, col.names) {
        url <- paste0(url_base, file)
        target <- paste0(savedir, "/", file)
        if (!file.exists(target)) {
          x <- RCurl::getURLContent(url)
          writeBin(as.vector(x), target)
        }
        if (file == "CTD_chemicals.tsv.gz") {
          # particular error
          data <- gs(readLines(target), "^([^\t]+\t[^\t]+).*", "\\1")
          data <- tibble::as_tibble(
            data.table::fread(text = data, select = select, col.names = col.names)
          )
        } else {
          data <- ftibble(target, sep = "\t", select = select, col.names = col.names)
        }
        data
      })
    names(lst) <- get_realname(files)
    if (any(names(lst) == "CTD_chemicals")) {
      ## add CID annotation into Chemicals data.frame
      anno <- get_ctd_compoundInfo()
      lst$CTD_chemicals <- dplyr::mutate(lst$CTD_chemicals,
        meshID = gs(ChemicalID, "^MESH:", ""))
      lst$CTD_chemicals <- map(
        lst$CTD_chemicals, "meshID",
        anno, "RegistryID", "CID"
      )
    }
    saveRDS(lst, rdata)
    return(lst)
  } else {
    readRDS(rdata)
  }
}

get_ctd_compoundInfo <- function(savedir = .prefix("ctd", "db"), reload = F)
{
  file <- paste0(savedir, "/meshID2sid.rds")
  if (!file.exists(file) || reload) {
    url <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceall/",
      "Comparative%20Toxicogenomics%20Database",
      "/xrefs/RegistryID/JSON"
    )
    res <- e(RCurl::getURL(url))
    res <- e(rjson::fromJSON(res))
    while (length(res) == 1) {
      res <- res[[1]]
    }
    res <- dplyr::bind_rows(lapply(res, unlist))
    ## get CID
    file_sid2cid <- paste0(savedir, "/sid2cid.rds")
    if (!file.exists(file_sid2cid)) {
      anno <- batch_get_cids.sids(res$SID, sleep = .2)
      saveRDS(anno, file_sid2cid)
    } else {
      anno <- readRDS(file_sid2cid)
    }
    res <- map(res, "SID", anno, "SID", "CID", col = "CID")
    saveRDS(res, file)
    res
  } else {
    readRDS(file)
  }
}
