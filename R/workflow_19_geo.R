# ==========================================================================
# workflow of geo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_geo <- setClass("job_geo", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    method = "GEO <https://www.ncbi.nlm.nih.gov/geo/> used for expression dataset aquisition"
    ))

job_geo <- function(id)
{
  .job_geo(object = id)
}

setMethod("step0", signature = c(x = "job_geo"),
  function(x){
    step_message("Prepare your data with function `job_geo`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_geo"),
  function(x){
    step_message("Get GEO metadata and information.")
    about <- e(GEOquery::getGEO(object(x)))
    metas <- get_metadata.geo(about)
    prods <- get_prod.geo(metas)
    x@params$about <- about
    x@params$metas <- metas
    x@params$prods <- prods
    guess <- metas$res[[1]]
    guess <- dplyr::rename_all(guess, make.names)
    guess <- dplyr::select(guess, 1:2, dplyr::ends_with(".ch1"))
    x@params$guess <- guess
    return(x)
  })

setMethod("step2", signature = c(x = "job_geo"),
  function(x, filter_regex = NULL){
    step_message("Download geo datasets.")
    e(GEOquery::getGEOSuppFiles(object(x), filter_regex = filter_regex))
    return(x)
  })

setMethod("meta", signature = c(x = "job_geo"),
  function(x, use = 1L){
    counts <- as_tibble(x@params$about[[ use ]]@assayData$exprs)
    genes <- as_tibble(x@params$about[[ use ]]@featureData@data)
    namel(counts, genes)
  })

setMethod("asjob_limma", signature = c(x = "job_geo"),
  function(x, metadata, use = 1L, normed = F){
    counts <- as_tibble(x@params$about[[ use ]]@assayData$exprs)
    genes <- as_tibble(x@params$about[[ use ]]@featureData@data)
    if (any(colnames(genes) == "rownames")) {
      if (grpl(colnames(genes), "Gene Symbol")) {
        genes <- dplyr::relocate(genes, rownames, hgnc_symbol = `Gene Symbol`)
      }
    } else {
      guess <- grpf(colnames(genes), "^gene", ignore.case = T)
      message("All available:\n\t", paste0(guess, collapse = ", "))
      message("Use the firist.")
      guess <- guess[[ 1 ]]
      keep <- !is.na(genes[[ guess ]]) & genes[[ guess ]] != "" & !duplicated(genes[[ guess ]])
      ## format
      genes <- dplyr::mutate(genes, rownames = !!rlang::sym(guess))
      genes <- dplyr::relocate(genes, rownames)
      message("Col.names of the data.frame:\n\t",
        stringr::str_trunc(paste0(colnames(counts), collapse = ", "), 40))
      counts <- dplyr::mutate(counts, rownames = !!genes[[ guess ]])
      counts <- dplyr::relocate(counts, rownames)
      genes <- genes[ keep, ]
      counts <- counts[ keep,  ]
    }
    if (normed) {
      counts <- dplyr::select(counts, -1)
      counts <- data.frame(counts)
      rownames(counts) <- genes$rownames
      job_limma_normed(counts, metadata)
    } else {
      job_limma(new_dge(metadata, counts, genes))
    }
  })
