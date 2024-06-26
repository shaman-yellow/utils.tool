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
    cite = "[@TcgabiolinksAColapr2015]",
    method = "R package `TCGAbiolinks` used for abtain TCGA dataset",
    tag = "db:TCGA",
    analysis = "TCGA 数据获取"
    ))

job_tcga <- function(project)
{
  object <- list(
    RNA = list(
      project = project, 
      data.category = "Transcriptome Profiling", 
      data.type = "Gene Expression Quantification", 
      workflow.type = "STAR - Counts"
      ),
    mutation = list(
      project = project,
      data.category = "Simple Nucleotide Variation", 
      access = "open",
      data.type = "Masked Somatic Mutation", 
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
      ),
    protein = list(
      project = project,
      data.category = "Proteome Profiling",
      data.type = "Protein Expression Quantification"
      ),
    clinical = list(
      project = project,
      data.category = "Clinical",
      data.type = "Clinical Supplement",
      data.format = "BCR XML"
    )
  )
  .job_tcga(object = object, params = list(project = project))
}

setMethod("step0", signature = c(x = "job_tcga"),
  function(x){
    step_message("Prepare your data with function `job_tcga`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_tcga"),
  function(x, query = c("RNA", "clinical"), keep_consensus = T){
    step_message("Get information in TCGA.")
    object(x) <- object(x)[ names(object(x)) %in% query ]
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
    if (length(query) > 1) {
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
    } else {
      cons <- NULL
    }
    res_query <- .set_lab(res_query, paste(sig(x), x$object), names(res_query), "metadata")
    x@tables[[ 1 ]] <- c(res_query, namel(cons))
    return(x)
  })

setMethod("filter", signature = c(x = "job_tcga"),
  function(x, ids, type = "RNA", use.tnbc = F){
    message("Filter the 'query' before downloading the data.")
    if (x@step != 1) {
      stop("x@step != 1")
    }
    if (use.tnbc) {
      if (!identical(x$project, "TCGA-BRCA")) {
        stop("The param `use.tnbc` only used for TCGA-BRCA.")
      }
      pb <- get_data.rot2016()
      ids <- dplyr::filter(object(pb), TNBC == "YES")$BARCODE
      pb@object <- NULL
      .add_internal_job(pb)
    }
    if (!missing(ids)) {
      ids <- substr(unique(ids), 1, 12)
    }
    for (i in type) {
      object(x)[[ type ]]$results[[ 1 ]] %<>%
        dplyr::filter(substr(cases, 1, 12) %in% !!ids)
    }
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
  function(x, use = "vital_status", query = "RNA",
    clinical.info = c("patient", "drug", "admin", "follow_up", "radiation"))
  {
    step_message("Prepare data for next step analysis.")
    lst.query <- object(x)[[ query ]]
    if (is.null(lst.query))
      stop("is.null(lst.query)")
    x@params$queries <- object(x)
    obj <- e(TCGAbiolinks::GDCprepare(query = lst.query))
    if (query == "RNA") {
      p.vital <- new_pie(SummarizedExperiment::colData(obj)[[ use ]])
      x@plots[[ 3 ]] <- namel(p.vital)
    } else if (query == "protein") {
      obj <- as_tibble(obj)
      if (!is.null(object(x)[[ "clinical" ]])) {
        clinical.info <- match.arg(clinical.info)
        x$metadata <- e(TCGAbiolinks::GDCprepare_clinic(
            query = object(x)[[ "clinical" ]],
            clinical.info = clinical.info))
        x$metadata <- as_tibble(x$metadata)
      }
    }
    object(x) <- obj
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
    # tumor, 01~09; normal, 10~19
    metadata <- dplyr::mutate(metadata, isTumor = ifelse(
        as.numeric(substr(rownames(metadata), 14, 15)) < 10,
        'tumor', 'normal'))
    metadata <- dplyr::mutate(metadata, group = !!rlang::sym(group))
    p.isTumor <- new_pie(metadata$isTumor)
    genes <- data.frame(object(x)@rowRanges)
    genes <- dplyr::relocate(genes, !!rlang::sym(row_id))
    counts <- object(x)@assays@data$unstranded
    colnames(counts) <- metadata[[ col_id ]]
    rownames(counts) <- genes[[ row_id ]]
    object <- e(edgeR::DGEList(counts, samples = metadata, genes = genes))
    project <- x$project
    x <- job_limma(object)
    if (get_treatment) {
      x <- .get_treatment.lm.tc(x)
      x <- .get_treatment.lm.tc(x, "R", name_suffix = ".radiation")
    }
    x <- meta(x, group)
    x$p.isTumor <- p.isTumor
    x$isTcga <- T
    x$project <- project
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


