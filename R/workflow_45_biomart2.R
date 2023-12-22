# ==========================================================================
# workflow of biomart2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_biomart2 <- setClass("job_biomart2", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/grimbough/biomaRt/issues/61"),
    cite = "[@MappingIdentifDurinc2009]",
    method = "The `biomart` was used for mapping genes between organism (e.g., mgi_symbol to hgnc_symbol)"
    ))

job_biomart2 <- function(values,
  from = c("hsapiens_gene_ensembl", "sscrofa_gene_ensembl", "mmusculus_gene_ensembl"),
  to = c("hsapiens_gene_ensembl", "sscrofa_gene_ensembl", "mmusculus_gene_ensembl"))
{
  dataset_map <- list("hsapiens_gene_ensembl" = "hgnc_symbol", "mmusculus_gene_ensembl" = "mgi_symbol")
  args <- list()
  args$from <- match.arg(from)
  args$to <- match.arg(to)
  args$fromAttri <- dplyr::recode(args$from, !!!dataset_map)
  args$toAttri <- dplyr::recode(args$to, !!!dataset_map)
  x <- .job_biomart2()
  x$args <- args
  object(x) <- rm.no(values)
  return(x)
}

setMethod("step0", signature = c(x = "job_biomart2"),
  function(x){
    step_message("Prepare your data with function `job_biomart2`.")
  })

setMethod("step1", signature = c(x = "job_biomart2"),
  function(x, group_number = 500, cl = 10,
    savepath = paste0("biomart2_", x$args$fromAttri, "_AS_", x$args$toAttri))
  {
    step_message("Query ...")
    if (!dir.exists(savepath)) {
      dir.create(savepath)
    }
    rdata <- paste0(savepath, "/", "mapping.rdata")
    db <- extract_rdata_list(rdata)
    if (!is.null(db)) {
      db <- frbind(db, fill = T)
      values <- object(x) %>% .[!. %in% db[[1]]]
    } else {
      values <- object(x)
    }
    while (length(values)) {
      groups <- .multi_tasks_urls.biomart2(values, group_number, x$args)
      urls <- groups$urls
      res <- pbapply::pblapply(names(urls), cl = cl,
        function(gn) {
          file <- paste0(savepath, "/", gn, ".tsv")
          cdRun("wget -t 0 -c -O ", file, " ",
            "'", urls[[ gn ]], "'")
        })
      db <- list(extra = db)
      packing_as_rdata_list(savepath, "^G.*\\.tsv$", get_filename(rdata), dedup = F, extra = db)
      db <- extract_rdata_list(rdata)
      db <- frbind(db, fill = T)
      # values <- object(x) %>% .[ !. %in% db[[1]] ]
      message("Check integrity...")
      notGet <- lapply(groups$theValues,
        function(values) {
          if (!any(values %in% db[[1]])) {
            values
          }
        })
      values <- unlist(notGet, use.names = F)
      message("Target of length ", length(values), " can retry.")
      isThat <- usethis::ui_yeah("Re-try?")
      if (!isThat) {
        values <- character(0)
      }
    }
    mapped <- data.frame(db[ db[[1]] %in% object(x), ])
    colnames(mapped) <- c(x$args$fromAttri, x$args$toAttri)
    x$mapped <- as_tibble(mapped)
    return(x)
  })

setMethod("step2", signature = c(x = "job_biomart2"),
  function(x, tops, use = c("P.Value", "adj.P.Val")){
    step_message("Add column in `tops`.")
    by <- colnames(x$mapped)[[1]]
    to <- colnames(x$mapped)[[2]]
    if (any(to == colnames(tops))) {
      tops <- dplyr::select(tops, -!!rlang::sym(to))
    }
    data <- tbmerge(tops, x$mapped, by = by)
    data <- dplyr::distinct(data, !!rlang::sym(to), .keep_all = T)
    use <- match.arg(use)
    data <- dplyr::arrange(data, !!rlang::sym(use))
    x$tops_mapped <- data
    return(x)
  })

.multi_tasks_urls.biomart2 <- function(values, group_number, args) {
  theValues <- grouping_vec2list(values, group_number, T)
  urls <- lapply(theValues,
    function(values) {
      theValue <- paste0(values, collapse = ",")
      url <- paste0('https://dec2021.archive.ensembl.org/biomart/martservice?query=<?xml version="1.0"
        encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"
        uniqueRows = "1" count = "0" datasetConfigVersion = "0.6"
        header="1" formatter = "TSV" requestid= "biomaRt">
        <Dataset name = "', args$from, '"><Attribute name = "', args$fromAttri, '"/>
        <Filter name = "', args$fromAttri, '" value = "', theValue, '" />
        </Dataset><Dataset name = "', args$to, '" ><Attribute name = "', args$toAttri, '"/>
        </Dataset></Query>'
      )
      gs(url, "\n|\\s+", " ")
    })
  names(urls) <- vapply(theValues, function(x) attr(x, "name"), character(1))
  namel(urls, theValues)
}
