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

.job_publish <- setClass("job_publish", 
  contains = c("job"),
  representation = representation(),
  prototype = prototype(
    method = "Supplementary file from article refer to"))

setClassUnion("job_PUBLISH", c("job_geo", "job_publish"))

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
    x@params$test <- list(
      genes = try(as_tibble(about[[ 1 ]]@featureData@data), T),
      counts = try(as_tibble(about[[ 1 ]]@assayData$exprs), T),
      anno = about[[ 1 ]]@annotation,
      anno.db = try(.get_biocPackage.gpl(about[[ 1 ]]@annotation), T)
    )
    return(x)
  })

.get_biocPackage.gpl <- function(gpl) {
  if (!is.character(gpl) | identical(gpl, character(0))) {
    stop("Incorrect character of gpl ID.")
  } else {
    data <- RCurl::getURL("https://gist.githubusercontent.com/seandavi/bc6b1b82dc65c47510c7/raw/b2027d7938896dce6145a1ebbcea75826813f6e1/platformMap.txt")
    data <- ftibble(data)
    dplyr::filter(data, gpl == !!gpl)$bioc_package
  }
}

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

get_metadata.geo <- function(lst,
  select = rlang::quos(rownames, title),
  pattern = c("diagnosis", "Sex", "^age", "^time point", "data_processing"),
  abbrev = c("data_processing"))
{
  res <- lapply(lst,
    function(eset){
      as_tibble(eset@phenoData@data)
    })
  if (!is.null(select)) {
    main <- lapply(res,
      function(data){
        cols <- colnames(data)
        extra <- unlist(.find_and_sort_strings(cols, pattern), use.names = F)
        dplyr::select(data, !!!select, dplyr::all_of(extra))
      })
    if (!is.null(abbrev)) {
      abbrev <- lapply(res,
        function(data){
          cols <- colnames(data)
          extra <- unlist(.find_and_sort_strings(cols, abbrev), use.names = F)
          dplyr::distinct(data, dplyr::pick(extra))
        })
    } else abbrev <- NULL
    res <- namel(main, abbrev, res)
  }
  lst_clear0(res)
}

get_prod.geo <- function(lst) {
  res <- as.list(dplyr::distinct(do.call(rbind, lst$abbrev)))
  .lich(res)
}
