# ==========================================================================
# workflow of reactome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_reactome <- setClass("job_reactome", 
  contains = c("job"),
  prototype = prototype(
    pg = "reactome",
    info = c("https://bioconductor.posit.co/packages/devel/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html"),
    cite = "",
    method = "",
    tag = "reactome",
    analysis = "ReactomeGSA 富集分析"
    ))

setGeneric("asjob_reactome",
  function(x, ...) standardGeneric("asjob_reactome"))

setMethod("asjob_reactome", signature = c(x = "job_seurat"),
  function(x, ref = NULL, compare.by = "group", sample.by = "orig.ident",
    cellType = x$group.by, ids = c(compare.by, sample.by, cellType))
  {
    object <- object(x)
    if (is.null(object@meta.data[[compare.by]])) {
      stop('is.null(object@meta.data[[compare.by]]).')
    }
    if (any(colnames(object@meta.data) == "ids")) {
      stop('any(colnames(object@meta.data) == "ids").')
    }
    object@meta.data <- tidyr::unite(
      object@meta.data, "ids", dplyr::all_of(ids), 
      sep = "___", remove = FALSE
    )
    metadata <- as_tibble(object@meta.data, idcol = "cell")
    e(Seurat::Idents(object) <- "ids")
    message(
      glue::glue("Clusters for enrichment: {less(levels(Seurat::Idents(object)), 10)}")
    )
    if (is.null(ref)) {
      fea <- as_feature(levels(Seurat::Idents(object)), "关键细胞", "cell")
    }
    snapAdd_onExit("x", "将{snap(fea)}进行通路富集分析。")
    message("Replace aggregation with a sparse matrix algorithm")
    cli::cli_alert_info(".reactome_prepare")
    object <- .reactome_prepare(object, verbose = TRUE)
    x <- .job_reactome(object = object)
    x$metadata <- metadata
    x$compare.by <- compare.by
    x$sample.by <- sample.by
    x$cellType <- cellType
    x$ids <- ids
    x$from <- "Seurat"
    return(x)
  })

setMethod("step0", signature = c(x = "job_reactome"),
  function(x){
    step_message("Prepare your data with function `job_reactome`.")
  })

setMethod("step1", signature = c(x = "job_reactome"),
  function(x){
    step_message("Send to reactome server.")
    cli::cli_alert_info("ReactomeGSA::perform_reactome_analysis")
    object(x) <- ReactomeGSA::perform_reactome_analysis(object(x), verbose = TRUE)
    x <- methodAdd(x, "反应组基因集分析（Reactome Gene Set Analysis，ReactomeGSA）是 Reactome 知识库旗下的多组学通路分析工具，可通过 Camera、PADOG 等方法将单细胞数据转化为通路水平信息，通过计算细胞集群平均表达量实现通路富集分析，识别细胞功能特征，解析细胞异质性。")
    x <- methodAdd(x, "使用 R 包 `ReactomeGSA` ({packageVersion('ReactomeGSA')}) 基于基因表达量对 Reactome 通路进行“活性评分”。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_reactome"),
  function(x, rerun = FALSE)
  {
    step_message("Gather for test.")
    compare.by <- x$compare.by
    if (!is.null(x$from) && x$from == "Seurat") {
      if (length(unique(x$metadata[[compare.by]])) != 2) {
        stop('length(unique(x$metadata[[gather.by]])) != 2, only comparison of two group allowed.')
      }
      paths <- ReactomeGSA::get_result(object(x), "pathways", "Seurat")
      genes <- ReactomeGSA::get_result(object(x), "fold_changes", "Seurat")
      fun_test <- function(x, y) {
        pbapply::pbvapply(seq_along(x),
          function(n) {
            wilcox.test(x[[n]], y[[n]])$p.value
          }, double(1))
      }
      lst <- setNames(list(paths, genes), c("pathways", "genes"))
      lapply_test <- function(lst) {
        lapply(lst,
          function(data) {
            whichIdCols <- which(vapply(data, is.character, logical(1)))
            if (!all(whichIdCols %in% 1:2)) {
              stop('!all(whichIdCols %in% 1:3).')
            }
            idCols <- colnames(data)[whichIdCols]
            if (!all(idCols %in% c("Pathway", "Name", "Identifier"))) {
              stop('!all(idCols %in% c("Pathway", "Name", "Identifier")).')
            }
            data <- dplyr::select(data, -dplyr::contains("Pathway"))
            data <- dplyr::rename(data, Name = 1)
            data <- tidyr::pivot_longer(
              data, dplyr::where(is.double), names_to = "type", values_to = "value"
            )
            data <- tidyr::separate(data, type, into = x$ids, sep = "___")
            # drop the 'sample' columns, for pivot_wider gathered to generate 'list'
            data <- dplyr::select(data, -!!rlang::sym(x$sample.by))
            data <- tidyr::pivot_wider(data, names_from = !!rlang::sym(compare.by), values_from = value)
            groups <- unique(x$metadata[[compare.by]])
            cli::cli_alert_info("wilcox.test")
            data <- dplyr::mutate(
              data, p.value = fun_test(!!rlang::sym(groups[1]), !!rlang::sym(groups[2]))
            )
            data <- dplyr::group_by(data, !!rlang::sym(x$cellType))
            data <- dplyr::mutate(
              data, p.adjust = p.adjust(p.value, method = "BH")
            )
            data <- dplyr::ungroup(data)
            dplyr::arrange(data, p.value)
          })
      }
      t.wilcox <- expect_local_data(
        "tmp", "reactome_wilcox", lapply_test, list(lst = lst), rerun = rerun
      )
      x <- tablesAdd(x, t.wilcox)
      x <- methodAdd(x, "以 wilcox.test 检验评估每个通路和基因在组间的表达差异显著性，并以 P value 对其排序。")
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_reactome"),
  function(x, top = 20){
    step_message("")
    lst <- x@tables$step2$t.wilcox
    vapply_mean <- function(x) {
      vapply(x, mean, double(1))
    }
    groups <- unique(x$metadata[[x$compare.by]])
    p.hps <- sapply(names(lst), simplify = FALSE,
      function(type) {
        data <- lst[[ type ]]
        fea <- head(unique(data$Name), n = top)
        data <- dplyr::filter(data, Name %in% !!fea)
        data <- dplyr::select(data, -p.value, -p.adjust)
        data <- dplyr::mutate(
          data, dplyr::across(dplyr::all_of(groups), vapply_mean)
        )
        data <- tidyr::pivot_longer(
          data, dplyr::all_of(groups), 
          names_to = "group", values_to = "value"
        )
        data <- dplyr::mutate(
          data, Cell = paste0(group, "_", !!rlang::sym(x$cellType)),
          value = scale(log2(value + 1), center = TRUE)[, 1]
        )
        args <- list(
          .data = data, .row = quote(Name), .column = quote(Cell),
          .value = quote(value), cluster_columns = FALSE, column_names_rot = 45,
          cluster_rows = TRUE, group_by = quote(group),
          row_names_max_width = grobWidth(textGrob(rownames(data), gpar(fontsize = 10, fontface = 1)))
        )
        fun_heatmap <- function(.data, .row, .column, .value, ..., group_by)
        {
          tidyHeatmap::annotation_tile(
            tidyHeatmap::heatmap(
              .data, .row = {{ .row }}, .column = {{ .column }}, 
              .value = {{ .value }}, ...
            ), .column = {{ group_by }}
          )
        }
        wrap_scale_heatmap(
          # ComplexHeatmap::Heatmap
          funPlot(fun_heatmap, args),
          data$Cell, data$Name, pre_width = max(nchar(data$Name)) * .1 + 2,
          pre_height = max(nchar(data$Name)) * .05 + 1,
          w.size = .3
        )
      })
    cn <- dplyr::recode(names(p.hps), genes = "基因", pathways = "通路")
    p.hps <- set_lab_legend(
      p.hps,
      glue::glue("{x@sig} top {names(p.hps)} ReactomeGSA heatmap"),
      glue::glue("ReactomeGSA 富集结果中各细胞类型按显著性排列靠前的{cn}")
    )
    x <- plotsAdd(x, p.hps)
    return(x)
  })

.reactome_prepare <- function (object, use_interactors = TRUE, include_disease_pathways = FALSE, 
  create_reactome_visualization = FALSE, create_reports = FALSE, 
  report_email = NULL, verbose = FALSE, assay = "RNA", 
  slot = "counts", ...)
{
  if (!assay %in% Seurat::Assays(object)) {
    stop("Error: Assay '", assay, "' does not exist in passed Seurat object.", 
      call. = FALSE)
  }
  raw_data <- Seurat::GetAssayData(object, assay = assay, layer = slot)
  cell_ids <- as.character(Seurat::Idents(object))
  if (length(unique(cell_ids)) < 2) {
    stop("Only one identification found: '", cell_ids[1], 
      "'. Please ensure that cell / cluster ids are stored as the primary identification (Ident) of your Seurat object.",
      " Clustering has to be performed prior to this pathway analysis.",
      call. = FALSE)
  }
  if (verbose) {
    message("Calculating average cluster expression...")
  }
  groups <- factor(cell_ids)
  n_groups <- nlevels(groups)
  require(Matrix)
  G <- Matrix::sparseMatrix(
    i = seq_along(groups),
    j = as.integer(groups),
    x = 1,
    dims = c(length(groups), n_groups)
  )
  sums <- raw_data %*% G
  group_sizes <- tabulate(groups, n_groups)
  av_counts <- sweep(sums, 2, group_sizes, "/")
  rownames(av_counts) <- rownames(raw_data)
  colnames(av_counts) <- levels(groups)
  av_counts <- data.frame(av_counts)
  request <- e(ReactomeGSA::ReactomeAnalysisRequest(method = "ssGSEA"))
  request <- e(ReactomeGSA::set_parameters(request, use_interactors = use_interactors, 
      include_disease_pathways = include_disease_pathways, 
      create_reactome_visualization = create_reactome_visualization, 
      create_reports = create_reports))
  if (!is.null(report_email)) {
    request <- ReactomeGSA::set_parameters(request, email = report_email)
  }
  cell_groups <- colnames(av_counts)
  request <- ReactomeGSA::add_dataset(request, expression_values = av_counts, 
    name = "Seurat", type = "rnaseq_counts", comparison_factor = "Cluster", 
    comparison_group_1 = unique(cell_groups)[1], comparison_group_2 = unique(cell_groups)[2], 
    sample_data = data.frame(row.names = cell_groups, 
      Cluster = cell_groups))
  request
}


