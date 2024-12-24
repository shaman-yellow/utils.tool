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
  x <- .job_tcga(object = object, params = list(project = project))
  x <- snapAdd(x, "获取 {project} 数据。")
  x
}

setMethod("step0", signature = c(x = "job_tcga"),
  function(x){
    step_message("Prepare your data with function `job_tcga`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_tcga"),
  function(x, query = c("RNA", "clinical"), keep_consensus = F){
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
          data[[ paste0(".row", "_", n) ]] <- seq_len(nrow(data))
          data
        })
      cons <- data_id[[1]]
      for (i in 2:length(data_id)) {
        cons <- merge(cons, data_id[[ i ]], by = ".id")
      }
      cons <- distinct(cons, .id, .keep_all = T)
      if (keep_consensus) {
        for (i in seq_along(object(x))) {
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
  function(x, dir = .prefix("GDCdata", "db")){
    step_message("Download data from TCGA.")
    e(lapply(object(x),
        function(query) {
          TCGAbiolinks::GDCdownload(query, directory = dir)
        }))
    x$dir <- dir
    meth(x)$step2 <- glue::glue("以 R 包 `TCGAbiolinks` ({packageVersion('TCGAbiolinks')}) {cite_show('TcgabiolinksAColapr2015')} 获取 {x$project} 数据集。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_tcga"),
  function(x, use = "vital_status", query = "RNA",
    clinical.info = c("patient", "drug", "admin", "follow_up", "radiation"),
    type = "unstranded")
  {
    step_message("Prepare data for next step analysis.")
    lst.query <- object(x)[[ query ]]
    if (is.null(lst.query))
      stop("is.null(lst.query)")
    x@params$queries <- object(x)
    obj <- e(TCGAbiolinks::GDCprepare(query = lst.query, directory = x$dir))
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
    obj <- object@assays@data[ names(object@assays@data %in% type) ]
    x$dim <- dim(obj)
    object(x) <- obj
    return(x)
  })

setMethod("merge", signature = c(x = "job_tcga", y = "job_tcga"),
  function(x, y, ...){
    if (!identical(object(x)@rowRanges, object(y)@rowRanges)) {
      stop("object(x) and object(y) do not has the same feature data.")
    }
    projects <- c(x$project, y$project)
    object(x)@colData$batch <- make.names(x$project)
    object(y)@colData$batch <- make.names(y$project)
    cols <- intersect(colnames(object(x)@colData), colnames(object(y)@colData))
    object(x)@colData <- rbind(object(x)@colData[, cols], object(y)@colData[, cols])
    object(x)@assays@data <- object(x)@assays@data[names(object(x)@assays@data) == "unstranded"]
    object(x)@assays@data$unstranded <- cbind(
      object(x)@assays@data$unstranded, object(y)@assays@data$unstranded
    )
    print(table(object(x)@colData$batch))
    x$dim <- dim(object(x))
    s.com <- try_snap(data.frame(object(x)@colData, check.names = F), "batch", "sample")
    x <- methodAdd(x, "将 {s.com} 数据集合并。")
    x <- snapAdd(x, "将 {s.com} 数据集合并。", step = "merge")
    x$project <- bind(projects)
    return(x)
  })

setGeneric("asjob_limma", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_limma"))

setMethod("asjob_limma", signature = c(x = "job_tcga"),
  function(x, ..., col_id = "sample", row_id = "gene_id", group = "vital_status",
    get_treatment = T)
  {
    step_message("Use `object(x)@assays@data$unstranded` converted as job_limma.")
    filter_sample <- rlang::enquos(...)
    metadata <- data.frame(object(x)@colData)
    metadata <- dplyr::relocate(metadata, !!rlang::sym(col_id))
    # tumor, 01~09; normal, 10~19
    metadata <- dplyr::mutate(metadata, isTumor = ifelse(
        as.numeric(
          gs(rownames(metadata), "[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-([0-9]+).*", "\\1")
          ) < 10, 'tumor', 'normal'))
    metadata <- dplyr::mutate(metadata, group = !!rlang::sym(group))
    genes <- data.frame(object(x)@rowRanges)
    genes <- dplyr::relocate(genes, !!rlang::sym(row_id))
    counts <- object(x)@assays@data$unstranded
    colnames(counts) <- metadata[[ col_id ]]
    rownames(counts) <- genes[[ row_id ]]
    if (length(filter_sample)) {
      message("Filter by: ", paste0(filter_sample, collapse = ", "))
      onrow <- nrow(metadata)
      metadata <- trace_filter(metadata, quosures = filter_sample)
      meta.snap <- snap(metadata)
      nnrow <- nrow(metadata)
      message(glue::glue("Filter out: {onrow - nnrow}"))
      counts <- counts[, colnames(counts) %in% metadata[[ col_id]] ]
    }
    if (any(duplicated(colnames(counts)))) {
      message('any(duplicated(colnames(counts))), Duplicated patients found, ',
        'remove thats herein (Prioritize retaining normal sample).')
      print(table(duplicated(colnames(counts))))
      if (!identical(colnames(counts), metadata$sample)) {
        stop('identical(colnames(counts), metadata$sample), ')
      }
      orders <- unique(c(which(metadata$isTumor == "normal"), seq_len(nrow(metadata))))
      counts <- counts[, orders]
      metadata <- metadata[orders, ]
      counts <- counts[ , !duplicated(colnames(counts)) ]
      metadata <- dplyr::distinct(metadata, sample, .keep_all = T)
    }
    object <- e(edgeR::DGEList(counts, samples = metadata, genes = genes))
    project <- x$project
    x <- job_limma(object)
    if (get_treatment && grpl(project, "^TCGA")) {
      x <- .get_treatment.lm.tc(x)
      x <- .get_treatment.lm.tc(x, "R", name_suffix = ".radiation")
    }
    x <- meta(x, group)
    p.group <- new_pie(metadata[[ group ]], "Group distribution")
    p.group <- .set_lab(p.group, sig(x), "Group distribution")
    p.isTumor <- new_pie(metadata$isTumor, "Tumor distribution")
    p.isTumor <- .set_lab(p.isTumor, sig(x), "Tumor distribution")
    if (grpl(project, "^TCGA")) {
      t.common <- dplyr::select(metadata, age_at_index, vital_status, gender, tumor_grade,
        ajcc_pathologic_stage, classification_of_tumor)
      x$t.common <- t.common
    }
    x$p.group <- p.group
    x$p.isTumor <- p.isTumor
    x$isTcga <- T
    x$project <- project
    if (length(filter_sample)) {
      x <- methodAdd(x, "{meta.snap}")
      x <- snapAdd(x, "{meta.snap}")
    }
    # if (!is.null(filter_follow_up)) {
    # meth(x)$step0 <- glue::glue("以 days_to_last_follow_up 大于 {filter_follow_up} 用于分析。")
    # }
    return(x)
  })

summary_tibble <- function(data, markdown = knitr::is_latex_output(), ...) {
  e(summarytools::st_options(use.x11 = FALSE))
  e(summarytools::dfSummary(data, style = if (markdown) "grid" else "multiline", ...))
}

summary_tibble2 <- function(data, group) {
  tabone <- tableone::CreateTableOne(data = data, strata = group)
  res <- as_tibble(tableone:::print.TableOne(tabone))
  colnames(res)[1] <- "Type"
  return(res)
}

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


