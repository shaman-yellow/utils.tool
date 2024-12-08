# ==========================================================================
# workflow of gds
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gds <- setClass("job_gds", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "gds",
    info = c("https://www.ncbi.nlm.nih.gov/books/NBK3837/"),
    cite = "",
    method = "",
    tag = "gds",
    analysis = "GSE 数据搜索"
    ))

job_gds <- function(keys,
  n = "6:1000", org = "Homo Sapiens",
  elements = c("GSE", "title", "summary", "taxon", "gdsType",
    "n_samples", "PubMedIds", "BioProject"),
  ...)
{
  query <- add_field(keys, "[Description]")
  if (is.integer(n)) {
    n <- paste0(min(n), ":", max(n))
  }
  extra <- glue::glue(
    "({n}[Number of Samples]) AND ({org}[Organism]) AND (GSE[Entry Type])"
  )
  object <- edirect_db("gds", c(query, extra), elements, ...)
  object <- dplyr::mutate(object, GSE = paste0("GSE", GSE))
  object <- .set_lab(object, keys[1], "EDirect query")
  x <- .job_gds(object = object, params = namel(query, elements, org, n))
  x <- methodAdd(x, "使用 Entrez Direct (EDirect) <https://www.ncbi.nlm.nih.gov/books/NBK3837/> 搜索 GEO 数据库 (`esearch -db gds`)，查询信息为: {query}。")
  x
}

add_field <- function(keys, field) {
  if (any(which <- !grpl(keys, "\\["))) {
    keys[which] <- vapply(keys[which],
      function(x) paste0(x, field), character(1)
    )
  }
  keys
}

edirect_db <- function(db, query, elements,
  cache = .prefix(paste0("query_", db), "db"),
  force = F)
{
  dir.create(cache, F)
  if (length(query)) {
    query <- paste0("(", query, ")")
    query <- paste0(query, collapse = " AND ")
  }
  target <- file.path(cache, make.names(query))
  if (!file.exists(target) || force) {
    script <- glue::glue("esearch -db {db} -query '{query}' |\\
      efetch -format docsum |\\
      xtract -pattern DocumentSummary \\
      -sep '|' -element {paste0(elements, collapse = ' ')} \\
      > {target}"
    )
    cdRun(script)  
  }
  object <- ftibble(target, header = F, sep = "\t")
  colnames(object) <- elements
  object
}

setMethod("active", signature = c(x = "job_gds"),
  function(x, n = 10){
    ids <- head(object(x)$GSE, n = n)
    urls <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ids)
    lapply(urls, browseURL)
  })

setMethod("vis", signature = c(x = "job_gds"),
  function(x, ref = NULL, n = 10){
    data <- object(x)
    if (!is.null(ref)) {
      data <- dplyr::filter(data, grpl(summary, !!ref))
    }
    print(data, n = n)
  })

setMethod("step0", signature = c(x = "job_gds"),
  function(x){
    step_message("Prepare your data with function `job_gds`.")
  })
