# ==========================================================================
# workflow of tcga
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tcga <- setClass("job_tcga", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("Tutorial: https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html",
      "\nhttps://portal.gdc.cancer.gov/"),
    cite = "[@TcgabiolinksAColapr2015]"
    ))

job_tcga <- function(project)
{
  object <- list(
    query_1 = list(
      project = project, 
      data.category = "Transcriptome Profiling", 
      data.type = "Gene Expression Quantification", 
      workflow.type = "STAR - Counts"
    ),
    query_2 = list(
      project = project,
      data.category = "Clinical",
      data.type = "Clinical Supplement",
      data.format = "BCR XML"
    )
  )
  .job_tcga(object = object)
}

setMethod("step0", signature = c(x = "job_tcga"),
  function(x){
    step_message("Prepare your data with function `job_tcga`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_tcga"),
  function(x, keep_consensus = T){
    step_message("Get information in TCGA.
      "
    )
    pblapply <- pbapply::pblapply
    if (is.null(x@tables$step1)) {
      object(x) <- e(pblapply(object(x),
          function(args) {
            do.call(TCGAbiolinks::GDCquery, args)
          }))
      res_query <- sapply(names(object(x)), simplify = F,
        function(name) {
          as_tibble(object(x)[[ name ]]$results[[1]])
        })
    } else {
      res_query <- x@tables$step1[ grepl("^query_", names(x@tables$step1)) ]
    }
    n <- 0
    data_id <- lapply(res_query,
      function(data) {
        n <<- n + 1
        data <- dplyr::select(data, .id = cases)
        data <- dplyr::mutate(data, .id = stringr::str_extract(.id, "[^\\-]*-[^\\-]*-[^\\-]*"))
        data[[ paste0(".row", "_", n) ]] <- 1:nrow(data)
        data
      })
    cons <- data_id[[1]]
    for (i in 2:length(data_id)) {
      cons <- merge(cons, data_id[[ i ]], by = ".id")
    }
    cons <- distinct(cons, .id, .keep_all = T)
    if (keep_consensus) {
      for (i in 1:length(object(x))) {
        object(x)[[i]]$results[[1]] %<>%
          dplyr::slice(dplyr::all_of(cons[[ paste0(".row_", i) ]]))
      }
    }
    x@tables[[ 1 ]] <- c(res_query, namel(cons))
    return(x)
  })

setMethod("step2", signature = c(x = "job_tcga"),
  function(x){
    step_message("Download data from TCGA.")
    e(lapply(object(x),
        function(query) {
          TCGAbiolinks::GDCdownload(query)
        }))
    return(x)
  })

setMethod("step3", signature = c(x = "job_tcga"),
  function(x, use = "vital_status"){
    step_message("Prepare 'Transcriptome Profiling' data (red{{only one query}}).")
    query <- lapply(object(x),
      function(query) {
        if (query$data.category == "Transcriptome Profiling")
          query
      })
    query <- lst_clear0(query)[[1]]
    x@params$queries <- object(x)
    object(x) <- e(TCGAbiolinks::GDCprepare(query = query))
    p.vital <- new_pie(SummarizedExperiment::colData(object(x))[[ use ]])
    x@plots[[ 3 ]] <- namel(p.vital)
    return(x)
  })

setGeneric("asjob_limma", 
  function(x, ...) standardGeneric("asjob_limma"))

setMethod("asjob_limma", signature = c(x = "job_tcga"),
  function(x, col_id = "sample", row_id = "gene_id", group = "vital_status",
    get_treatment = T)
  {
    step_message("Use `object(x)@assays@data$unstranded` converted as job_limma.")
    metadata <- data.frame(object(x)@colData)
    metadata <- dplyr::relocate(metadata, !!rlang::sym(col_id))
    metadata <- dplyr::mutate(metadata, group = vital_status)
    genes <- data.frame(object(x)@rowRanges)
    genes <- dplyr::relocate(genes, !!rlang::sym(row_id))
    counts <- object(x)@assays@data$unstranded
    colnames(counts) <- metadata[[ col_id ]]
    rownames(counts) <- genes[[ row_id ]]
    object <- e(edgeR::DGEList(counts, samples = metadata, genes = genes))
    x <- job_limma(object)
    if (get_treatment) {
      x <- .get_treatment.lm.tc(x)
      x <- .get_treatment.lm.tc(x, "R", name_suffix = ".radiation")
    }
    x <- meta(x, group)
    return(x)
  })

.get_treatment.lm.tc <- function(x,
  type = c("Pharmaceutical Therapy, NOS", "Radiation Therapy, NOS"),
  attr = c("treatment_or_therapy"), add_into = T, name_suffix = NULL)
{
  type <- match.arg(type)
  treats <- select(x@object$samples, dplyr::contains("treat"))
  data <- lapply(treats$treatments,
    function(data) {
      data <- filter(data, treatment_type == !!type)
      data <- select(data, treatment_type, dplyr::all_of(attr))
      data
    })
  data <- do.call(dplyr::bind_rows, data)
  if (!is.null(name_suffix)) {
    data <- dplyr::rename_all(data, function(x) paste0(x, name_suffix))
  }
  if (!add_into)
    return(data)
  x@object$samples <- dplyr::bind_cols(x@object$samples, data)
  return(x)
}


