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
    method = "The Comparative Toxicogenomics Database (CTD) used for finding relationship between chemicals and disease"
    ))

job_ctd <- function()
{
  .job_ctd()
}

setMethod("step0", signature = c(x = "job_ctd"),
  function(x){
    step_message("Prepare your data with function `job_ctd`.")
  })

setMethod("step1", signature = c(x = "job_ctd"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

get_ctd_data <- function(savedir = .prefix("ctd", "db"), reload = F) {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  rdata <- paste0(savedir, "/all.rds")
  if (!file.exists(rdata) || reload) {
    url_base <- "https://ctdbase.org/reports/"
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
    saveRDS(lst, rdata)
    return(lst)
  } else {
    readRDS(rdata)
  }
}
