# ==========================================================================
# workflow of esearch
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_esearch <- setClass("job_esearch", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://www.ncbi.nlm.nih.gov/books/NBK179288/"),
    cite = ""
    ))

job_esearch <- function(key, jour = .jour_bioinf(), in_title_and_abstract = T)
{
  if (in_title_and_abstract) {
    key <- paste0(key, " [TIAB] ")
  }
  data <- esearch.mj(key, jour)
  .job_esearch(object = data)
}

setMethod("vis", signature = c(x = "job_esearch"),
  function(x, ref = NULL, n = 100, col = c(".id", "Title")){
    data <- object(x)
    if (!is.null(ref)) {
      data <- dplyr::filter(data, grpl(.id, !!ref))
    }
    data <- dplyr::select(data, dplyr::all_of(col))
    print(data, n = n)
  })

setMethod("step0", signature = c(x = "job_esearch"),
  function(x){
    step_message("Prepare your data with function `job_esearch`.
      "
    )
  })

# map_from <- function(lst, db_names, db_values){}

.jour_bioinf <- function() {
  c("Nature Genetics", "Genome Biology", "Genome Research", "Nucleic Acids Res",
    "Briefings in Bioinformatics", "BMC genomics", "Nature Methods", "Nature Biotechnology")
}

esearch.mj <- function(key, jour = .jour_bioinf(), rbind = T)
{
  query <- paste(jour, "[JOUR] AND", key)
  sear <- pbapply::pblapply(query, esearch)
  names(sear) <- jour
  if (rbind) {
    sear <- tibble::as_tibble(data.table::rbindlist(sear, idcol = T, fill = T))
    if (nrow(sear) > 0)
      sear <- dplyr::arrange(sear, dplyr::desc(SortPubDate))
  }
  sear
}

esearch <- function(query = NULL, fetch.save = paste0(gsub(" ", "_", query), ".xml"),
  path = "search", tract.save = "res.tsv",
  fields = c("SortPubDate", "Title", "FullJournalName", "Name", "Id"))
{
  if (!dir.exists(path)) {
    dir.create(path)
  }
  if (!is.null(query)) {
    if (!file.exists(paste0(path, "/", fetch.save))) {
      cdRun("esearch -db pubmed -query \"", query, "\"",
        " | efetch -format docsum > ", fetch.save,
        path = path)
    }
  } else if (fetch.save == ".xml") {
    stop("fetch.save == \".xml\"")
  }
  cdRun(" cat ", fetch.save, " | xtract -pattern DocumentSummary ",
    " -sep \"|\" -element ", paste0(fields, collapse = " "),
    " > ", tract.save, path = path)
  file <- paste0(path, "/", tract.save)
  res <- ftibble(file, sep = "\t", quote = "", fill = T, header = F)
  if (nrow(res > 0)) {
    colnames(res) <- fields
    if (any("SortPubDate" == fields)) {
      res <- dplyr::mutate(res, SortPubDate = as.Date(SortPubDate))
      res <- dplyr::arrange(res, dplyr::desc(SortPubDate))
    }
    if (any("Id" == fields)) {
      res <- dplyr::mutate(res, Id = as.character(Id))
    }
  }
  res
}


