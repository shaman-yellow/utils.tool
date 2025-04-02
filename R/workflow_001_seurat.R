# ==========================================================================
# workflow of seurat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat <- setClass("job_seurat", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html"),
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019]",
    method = "The R package `Seurat` used for scRNA-seq processing",
    tag = "scrna:anno",
    analysis = "Seurat 单细胞数据分析"
    ))

job_seurat <- function(dir = NULL, project = basename(sub("/$", "", dir)),
  min.cells = 3, min.features = 200, file_h5 = NULL, ...)
{
  if (!is.null(file_h5)) {
    data <- e(Seurat::Read10X_h5(file_h5))
  } else {
    data <- e(Seurat::Read10X(dir))
  }
  object <- e(Seurat::CreateSeuratObject(counts = data, project = project,
      min.cells = min.cells, min.features = min.features, ...))
  x <- .job_seurat(object = object)
  meth(x)$step0 <- glue::glue("")
  x
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    step_message("Prepare your data with function `job_seurat`. ",
      crayon::red("Required: Seurat v5. "),
      "The `dir` to read would passed to `Seurat::Read10X` to ", 
      "load single-cell data. The directory must contains three ",
      "files: ", crayon::red("barcodes.tsv"), "; ",
      crayon::red("genes.tsv"), "; ", crayon::red("matrix.mtx"), "."
    )
  })

setMethod("step1", signature = c(x = "job_seurat"),
  function(x){
    step_message("QC and selecting cells for further analysis.",
      "This do:",
      "Generate column in `object(x)@meta.data`; ",
      "Plots in `x@plots[[ 1 ]]`.",
      "Check the plots then decided the ",
      crayon::red("min nFeature_RNA"), " and ",
      crayon::red("max nFeature_RNA"), ", ",
      "as well as ", crayon::red("percent.mt"), "."
    )
    object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(
      object(x), pattern = "^MT-"
    ))
    p.qc_pre <- plot_qc.seurat(x)
    p.qc_pre <- .set_lab(p.qc_pre, sig(x), "Pre-Quality control")
    p.qc_pre <- setLegend(p.qc_pre, "为 QC (质量控制) 图 (数据过滤前) 。")
    x$p.qc_pre <- p.qc_pre
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    x <- methodAdd(x, "使用 Seurat R 包 ({packageVersion('Seurat')}) 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000,
    use = "nFeature_RNA", sct = FALSE, ndims = 20)
  {
    step_message("This contains several execution: Subset the data...
      NOTE: inspect the plots and red{{Determined dims}} for downstream analysis. "
    )
    if (missing(min.features) | missing(max.features)) {
      stop("missing(min.features) | missing(max.features)")
    }
    object(x) <- e(SeuratObject:::subset.Seurat(
        object(x), subset = !!rlang::sym(use) > min.features &
          !!rlang::sym(use) < max.features & percent.mt < max.percent.mt
        ))
    x <- methodAdd(x, "一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。")
    x <- snapAdd(x, "前期质量控制，一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。过滤后，数据集共包含 {dim(x@object)[1]} 个基因，{dim(x@object)[2]} 个细胞。")
    p.qc_aft <- plot_qc.seurat(x)
    p.qc_aft <- .set_lab(p.qc_aft, sig(x), "After Quality control")
    p.qc_aft <- setLegend(p.qc_aft, "为数据过滤后的 QC 图。")
    x$p.qc_aft <- p.qc_aft
    if (sct) {
      object(x) <- e(Seurat::SCTransform(
          object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE,
          assay = SeuratObject::DefaultAssay(object(x))
          ))
    } else {
      object(x) <- e(Seurat::NormalizeData(object(x)))
      object(x) <- e(Seurat::FindVariableFeatures(object(x)))
      object(x) <- e(Seurat::ScaleData(object(x)))
    }
    object(x) <- e(Seurat::RunPCA(
        object(x), 
        features = SeuratObject::VariableFeatures(object(x)), reduction.name = "pca"
        ))
    p.pca_rank <- e(Seurat::ElbowPlot(object(x), ndims))
    p.pca_rank <- wrap(pretty_elbowplot(p.pca_rank), 4, 4)
    p.pca_rank <- .set_lab(p.pca_rank, sig(x), "Standard deviations of PCs")
    p.pca_rank <- setLegend(p.pca_rank, "为主成分 (PC) 的 Standard deviations。")
    x@plots[[ 2 ]] <- namel(p.pca_rank)
    x <- methodAdd(x, "执行标准 Seurat 分析工作流 (`NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`)。以 `ElbowPlot` 判断后续分析的 PC 维度。")
    x <- snapAdd(x, "数据归一化，PCA 聚类 (Seurat 标准工作流，见方法章节) 后，绘制 PC standard deviations 图。")
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat"),
  function(x, dims = 1:15, resolution = 1.2, features = NULL, reset = FALSE,
    reduction = if (is.null(features)) "pca" else "features_pca", force = FALSE)
  {
    step_message("This contains several execution: Cluster the cells, UMAP...")
    if (!force) {
      if (ncol(object(x)) > 3000 & resolution < 1.2) {
        stop("ncol(object(x)) > 3000 but resolution < 1.2. ",
          "Please consider higher `resolution`. ")
      }
    }
    if (reset) {
      object(x) <- e(Seurat::FindVariableFeatures(object(x)))
      object(x) <- e(Seurat::ScaleData(object(x)))
    }
    if (!is.null(features)) {
      message(glue::glue("Run PCA in specified features."))
      object(x) <- e(
        Seurat::RunPCA(object(x), features = unique(features), reduction.name = reduction)
      )
    }
    object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = reduction))
    object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution))
    object(x) <- e(Seurat::RunUMAP(object(x), dims = dims, reduction = reduction))
    p.umap <- e(Seurat::DimPlot(object(x), cols = color_set(TRUE)))
    p.umap <- wrap(p.umap, 6, 5)
    p.umap <- .set_lab(p.umap, sig(x), "UMAP", "Clustering")
    x@plots[[ 3 ]] <- namel(p.umap)

    x <- snapAdd(x, "在 1-{max(dims)} PC 维度，{resolution} 分辨率下，对细胞群 UMAP 聚类。")
    x <- methodAdd(x, "在 1-{max(dims)} PC 维度下，以 Seurat::FindNeighbors 构建 Nearest-neighbor Graph。随后在 {resolution} 分辨率下，以 Seurat::FindClusters 函数识别细胞群并以 Seurat::RunUMAP 进行 UMAP 聚类。")
    return(x)
  })

setMethod("step4", signature = c(x = "job_seurat"),
  function(x, use = "SingleR", use.level = c("label.main", "label.fine"),
    ref = celldex::HumanPrimaryCellAtlasData())
  {
    if (use == "scAnno") {
      ## PMID: 37183449
      ## https://github.com/liuhong-jia/scAnno
      stop("Deprecated. Too many bugs in that package.")
    } else if (use == "SingleR") {
      step_message("Use `SingleR` and `celldex` to annotate cell types.
        By default, red{{`celldex::HumanPrimaryCellAtlasData`}} was used
        as red{{`ref`}} dataset. This annotation would generate red{{'SingleR_cell'}}
        column in `object(x)@meta.data`. Plots were generated in `x@plots[[ 4 ]]`;
        tables in `x@tables[[ 4 ]]`.
        "
      )
      ref <- e(celldex::HumanPrimaryCellAtlasData())
      clusters <- object(x)@meta.data$seurat_clusters
      use.level <- match.arg(use.level)
      anno_SingleR <- e(SingleR::SingleR(object(x)@assays$SCT@scale.data,
          ref = ref, labels = ref[[ use.level ]], clusters = clusters
          ))
      score <- as.matrix(anno_SingleR$scores)
      rownames(score) <- rownames(anno_SingleR)
      p.score_SingleR <- callheatmap(
        new_heatdata(as_data_long(score, 
            row_var = "Cluster", col_var = "Cell_type"))
      )
      p.score_SingleR <- wrap(p.score_SingleR, 30, 8)
      anno_SingleR <- tibble::tibble(
        seurat_clusters = rownames(anno_SingleR),
        SingleR_cell = anno_SingleR$labels
      )
      object(x)@meta.data$SingleR_cell <-
        anno_SingleR$SingleR_cell[match(clusters,
          anno_SingleR$seurat_clusters)]
      p.map_SingleR <- e(Seurat::DimPlot(
          object(x), reduction = "umap", label = FALSE, pt.size = .7,
          group.by = "SingleR_cell", cols = color_set()
          ))
      p.map_SingleR <- p.map_SingleR
      p.map_SingleR <- wrap(p.map_SingleR, 10, 7)
      x@tables[[ 4 ]] <- namel(anno_SingleR)
      x@plots[[ 4 ]] <- namel(p.score_SingleR, p.map_SingleR)
      x@params$group.by <- "SingleR_cell"
      x <- methodAdd(x, "以 R 包 `SingleR` ({packageVersion('SingleR')}) 注释细胞群。")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_seurat"),
  function(x, workers = NULL, min.pct = .25, logfc.threshold = .25, 
    force = FALSE, topn = 5)
  {
    step_message("Find all Marders for Cell Cluster.")
    if (is.remote(x)) {
      file_markers <- file.path(x$map_local, filename_markers <- "markers.csv")
      if (!file.exists(file_markers) || force) {
        if (!rem_file.exists(filename_markers) || force) {
          if (is.null(workers)) {
            stop('is.null(workers).')
          }
          file_obj <- file.path(
            x$map_local, filename_obj <- paste0("seurat_5.rds")
          )
          if (!file.exists(file_obj) || force) {
            message(glue::glue("Save Seurat object as {file_obj}..."))
            saveRDS(object(x), file_obj)
          }
          if (!rem_file.exists(filename_obj) || force) {
            cdRun(glue::glue("scp {file_obj} {x$remote}:{x$wd}"))
          }
          script <- expression({
            require(Seurat)
            require(future)
            args <- commandArgs(trailingOnly = TRUE)
            ## args
            file_rds <- args[1]
            output <- args[2]
            workers <- as.integer(args[3])
            min.pct <- as.double(args[4])
            logfc.threshold <- as.double(args[5])
            ## run
            object <- readRDS(file_rds)
            options(future.globals.maxSize = Inf)
            future::plan(future::multicore, workers = workers)
            markers <- Seurat::FindAllMarkers(object, min.pct = min.pct,
              logfc.threshold = logfc.threshold, only.pos = TRUE)
            write.table(
              markers, output, sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE
            )
          })
          script <- as.character(script)
          if (force) {
            message(glue::glue("force == TRUE, remove remote and local file."))
            rem_file.remove(filename_markers)
            file.remove(file_markers)
          }
          rem_run(glue::glue(
              "{pg('Rscript', is.remote(x))} -e '{script}' \\
              {filename_obj} {filename_markers} {workers} {min.pct} {logfc.threshold}"
              ))
          testRem_file.exists(x, filename_markers, 1, later = FALSE)
        }
        file_markers <- get_file_from_remote(
          filename_markers, x$wd, file_markers, x$remote
        )
      }
      markers <- ftibble(file_markers)
      markers <- dplyr::mutate(markers, rownames = gene, .before = 1)
    } else {
      fun <- function(x) {
        markers <- as_tibble(
          e(Seurat::FindAllMarkers(object(x), min.pct = min.pct,
              logfc.threshold = logfc.threshold, only.pos = TRUE))
        )
      }
      if (!is.null(workers)) {
        markers <- parallel(x, fun, workers)
      } else {
        markers <- fun(x)
      }
    }
    all_markers_no_filter <- markers
    markers <- dplyr::filter(markers, p_val_adj < .05)
    markers <- .set_lab(markers, sig(x), "significant markers of cell clusters")
    markers <- setLegend(markers, "为所有细胞群的 Marker (LogFC 阈值 {logfc.threshold}; 最小检出率 {min.pct}; 矫正 P 值阈值 {0.05})")
    if (TRUE) {
      tops <- dplyr::slice_max(dplyr::group_by(markers, cluster), avg_log2FC, n = topn)
      p.toph <- e(Seurat::DoHeatmap(object(x), features = tops$gene, raster = TRUE))
      p.toph <- wrap(p.toph, 14, 12)
      x@plots[[ 5 ]] <- namel(p.toph)
    }
    x@tables[[ 5 ]] <- list(all_markers = markers, all_markers_no_filter = all_markers_no_filter)
    x <- snapAdd(x, "计算所有细胞群的 Marker。")
    x <- methodAdd(x, "以 `Seurat::FindAllMarkers` (LogFC 阈值 {logfc.threshold}; 最小检出率 {min.pct}) 为所有细胞群寻找 Markers。")
    return(x)
  })

setMethod("step6", signature = c(x = "job_seurat"),
  function(x, tissue, ref.markers = NULL, method = c("cellMarker", "gpt", "scsa"),
    # scsa
    org = c("Human", "Mouse"), filter.p = 0.01, filter.fc = .5, res.col = "scsa_cell",
    # cellMarker
    type = c("Normal cell"),
    # plot markers
    exclude_pattern = "derived|progenitor|Transitional|Memory|switch|white blood cell",
    exclude = NULL, include = NULL, show = NULL, notShow = NULL,
    reset = NULL, keep_markers = 3,
    # chatGPT 
    filter.pct = .25, toClipboard = TRUE, post_modify = FALSE,
    n = 30, variable = FALSE, hp_type = c("pretty", "seurat"), ...)
  {
    args <- c(as.list(environment()), list(...))
    args$init <- FALSE
    do.call(anno, args)
  })

setMethod("anno", signature = c(x = "job_seurat"),
  function(x, tissue, ref.markers = NULL, method = c("cellMarker", "gpt", "scsa"),
    # scsa
    org = c("Human", "Mouse"), filter.p = 0.01, filter.fc = .5, res.col = "scsa_cell",
    # cellMarker
    type = c("Normal cell"),
    # plot markers
    exclude_pattern = "derived|progenitor|Transitional|Memory|switch|white blood cell",
    exclude = NULL, include = NULL, show = NULL, notShow = NULL,
    reset = NULL, keep_markers = 3,
    # chatGPT 
    filter.pct = .25, toClipboard = TRUE, post_modify = FALSE,
    n = 30, variable = FALSE, hp_type = c("pretty", "seurat"), 
    ..., init = TRUE)
  {
    if (x@step < 5L) {
      stop('x@step < 5L.')
    }
    if (init) {
      x <- init(x)
      message(glue::glue("Set step as 6."))
      x@step <- 6L
    }
    method <- match.arg(method)
    if (method == "gpt") {
      step_message("Use ChatGPT (ChatGPT-4o) to annotate cell types.")
      all_markers <- x@tables$step5$all_markers
      if (!is.null(filter.pct)) {
        all_markers <- dplyr::filter(all_markers, pct.1 >= filter.pct)
      }
      if (variable) {
        all_markers <- dplyr::filter(
          all_markers, gene %in% SeuratObject::VariableFeatures(object(x))
        )
      }
      query <- prepare_GPTmessage_for_celltypes(
        tissue, all_markers, n = n, toClipboard = toClipboard
      )
      if (sureThat("Get results from clipboard?")) {
        feedback <- get_clipboard()
        if (length(feedback) != query$ncluster) {
          stop(glue::glue("The res number ({nrow(res)}) not match cluster number ({query$ncluster})."))
        }
        res <- parse_GPTfeedback(feedback)
        seqs <- match(object(x)@meta.data[[ "seurat_clusters"]], 0:(query$ncluster - 1))
        object(x)@meta.data[[ "ChatGPT_cell"]] <- res$cells[ seqs ]
        if (post_modify) {
          x <- .modify_gpt_marker_and_cell(
            x, all_markers, query, res, keep_markers
          )
          validMarkers <- x$gpt_excluOverMarkers
        } else {
          validMarkers <- unique(reframe_col(res$data, "markers", unlist)$markers)
        }
        ## heatmap
        hp_type <- match.arg(hp_type)
        if (hp_type == "seurat") {
          p.markers <- e(Seurat::DoHeatmap(object(x), features = validMarkers,
              group.by = "ChatGPT_cell", raster = TRUE, group.colors = color_set(), label = FALSE))
        } else if (hp_type == "pretty") {
          p.markers <- .plot_marker_heatmap(
            x, validMarkers, "ChatGPT_cell"
          )
        }
        p.markers <- .set_lab(
          wrap(p.markers, 7, 7.5), sig(x), "Markers in cell types"
        )
        p.markers <- setLegend(p.markers, "为 ChatGPT 注释细胞群使用的首要 Marker 热图。")
        attr(p.markers, "lich") <- new_lich(
          namel(
            ChatGPT_Query = gs(query$query, "\n", "\t\t\n"), 
            ChatGPT_feedback = bind(feedback, co = "\t\t\n")
          ), sep = "\n"
        )
        t.gptRes <- setLegend(res$data, "为 ChatGPT-4 对细胞类型的注释。")
        x <- tablesAdd(x, t.ChatGPT4_cell_types_annotation = t.gptRes)
        ## dim plot
        p.map_gpt <- e(Seurat::DimPlot(
            object(x), reduction = "umap", label = FALSE, pt.size = .7,
            group.by = "ChatGPT_cell", cols = color_set()
            ))
        p.map_gpt <- wrap(as_grob(p.map_gpt), 7, 4)
        p.map_gpt <- .set_lab(p.map_gpt, sig(x), "ChatGPT", "Cell type annotation")
        p.map_gpt <- setLegend(p.map_gpt, "ChatGPT 对细胞群的注释结果。")
        p.props_gpt <- plot_cells_proportion(
          object(x)@meta.data, "orig.ident", "ChatGPT_cell"
        )
        p.props_gpt <- .set_lab(p.props_gpt, sig(x), "ChatGPT4 Cell Proportions in each sample")
        p.props_gpt <- setLegend(p.props_gpt, "为 ChatGPT-4 注释的细胞群在各个样本中的占比。")
        x <- plotsAdd(x, p.map_gpt, p.markers, p.props_gpt)
        x <- snapAdd(
          x, "根据细胞群 Markers (检出率至少为 {filter.pct}，选取 Top {n}) ，让 ChatGPT-4 对细胞类型注释。"
        )
        x <- methodAdd(
          x, "以 ChatGPT-4 注释细胞类型 {cite_show('Assessing_GPT_4_Hou_W_2024')}。将每一个细胞群的 Top {n} 基因提供给 ChatGPT，使其注释细胞类型。询问信息为：\n\n{text_roundrect(query$message, 'docx')}"
        )
        x@params$group.by <- "ChatGPT_cell"
      } else {
        stop("Terminated.")
      }
    } else if (method == "scsa" || method == "cellMarker") {
      step_message("Use SCSA to annotate cell types (<https://github.com/bioinfo-ibms-pumc/SCSA>).")
      # lst <- do.call(scsa_annotation, as.list(environment()))
      if (method == "cellMarker" && is.null(ref.markers)) {
        job_markers <- job_markers(
          tissue, type = type, use_numbering = FALSE, 
          org = org
        )
        tissue <- bind(unique(job_markers$db$tissueType))
        ref.markers <- ref(job_markers)
        if (!is.null(include)) {
          ref.markers <- dplyr::filter(
            ref.markers, grpl(cell, include, TRUE)
          )
        }
        if (!is.null(exclude_pattern) || !is.null(exclude)) {
          exclude_pattern <- paste0(
            c(exclude_pattern, exclude), collapse = "|"
          )
          ref.markers <- dplyr::filter(
            ref.markers, !grpl(cell, exclude_pattern, TRUE)
          )
        }
      } else {
        job_markers <- NULL
      }
      lst <- scsa_annotation(
        x = x, tissue = if (method == "scsa") tissue else "All",
        ref.markers = ref.markers, filter.p = filter.p,
        filter.fc = filter.fc, org = org, reset = reset, res.col = res.col
      )
      x <- lst$x
      x <- tablesAdd(x, scsa_res_all = lst$scsa_res_all)
      x@params$group.by <- lst$res.col
      p.props_scsa <- plot_cells_proportion(
        object(x)@meta.data, "orig.ident", "scsa_cell"
      )
      p.props_scsa <- .set_lab(p.props_scsa, sig(x), "SCSA Cell Proportions in each sample")
      p.props_scsa <- setLegend(p.props_scsa, "为 SCSA 注释的细胞群在各个样本中的占比。")
      p.markers <- NULL
      if (method == "cellMarker" && is.null(ref.markers)) {
        ## .plot_marker_heatmap
        validMarkers <- dplyr::filter(ref.markers, cell %in% lst$scsa_res$Cell.Type)
        x <- map(
          x, job_markers, markers = split(
            validMarkers$markers, validMarkers$cell
          ), group.by = "scsa_cell", max = keep_markers, show = show, notShow = notShow
        )
        p.markers <- x$p.cellMarker
      } else if (!is.null(ref.markers)) {
        p.markers <- .plot_marker_heatmap(
          x, split(ref.markers$markers, ref.markers$cell), "scsa_cell",
          show = show, max = keep_markers, notShow = notShow, ...
        )
        p.markers <- set_lab_legend(
          wrap(p.markers, 7, 7),
          paste(sig(x), "Marker Validation"),
          "使用特异性 Marker 对细胞注释结果的验证热图。"
        )
      }
      x <- plotsAdd(
        x, p.map_scsa = lst$p.map_scsa, p.props_scsa = p.props_scsa, 
        p.markers = p.markers
      )
      if (method == "scsa" || !is.null(ref.markers)) {
        x <- snapAdd(x, "{if (!is.null(ref.markers)) '使用特异性 Marker，' else ''}以 SCSA 对细胞群注释。")
        x <- methodAdd(
          x, "以 Python 工具 `SCSA` ({cite_show('ScsaACellTyCaoY2020')}) (<https://github.com/bioinfo-ibms-pumc/SCSA>) 对细胞群注释。"
        )
      } else if (method == "cellMarker") {
        x$job_markers <- job_markers
        x <- snapAdd(x, "使用 CellMarker 数据库的 Marker 基因 (Tissue: {tissue})，以 SCSA 对细胞群注释。")
        x <- methodAdd(
          x, "以 Python 工具 `SCSA` ({cite_show('ScsaACellTyCaoY2020')}) (<https://github.com/bioinfo-ibms-pumc/SCSA>)，结合 CellMarker 数据库 ({cite_show('CellMarker_a_m_Zhang_2019')})，对细胞群注释。"
        )
      }
    }
    return(x)
  })

setMethod("step7", signature = c(x = "job_seurat"),
  function(x, classifier, db = org.Hs.eg.db::org.Hs.eg.db){
    step_message("
      Use `garnett::classify_cells` to anntate cells.
      Prarameter red{{`classifier`}} specify the pre-difined classifier.
      "
    )
    object <- e(SeuratWrappers::as.cell_data_set(object(x),
        group.by = x@params$group.by, ...))
    object <- e(garnett::classify_cells(object, classifier, db = db, 
      cds_gene_id_type = "SYMBOL"))
    object(x) <- e(SeuratWrappers:::as.Seurat.cell_data_set(object))
    x@params$group.by <- "cell_type"
    return(x)
  })

as_type_group <- function(type, group) {
  group <- paste0(type, "_", gs(group, "[0-9]*$", ""))
  gs(make.names(group), "[.]+", "_")
}

setMethod("diff", signature = c(x = "job_seurat"),
  function(x, group.by, contrasts, name = "contrasts",
    cut.fc = .3, cut.p = .5, min.pct = .25, force = FALSE)
  {
    message("Differential analysis for cells.")
    if (is.null(object(x)@meta.data[[ group.by ]])) {
      stop('is.null(object(x)@meta.data[[ group.by ]]).')
    }
    if (is(contrasts, "character")) {
      contrasts <- .get_versus_cell(contrasts, NULL, unlist = FALSE)
    }
    if (!all(unlist(contrasts) %in% ids(x, group.by))) {
      stop('!all(unlist(contrasts) %in% ids(x, group.by)).')
    }
    if (is.data.frame(contrasts)) {
      contrasts <- apply(contrasts, 1, c, simplify = FALSE)
    }
    if (is(contrasts, "list") && !all(lengths(contrasts) == 2)) {
      stop('is(contrasts, "list") && !all(lengths(contrasts) == 2).')
    }
    numCells <- table(x@object@meta.data[[ group.by ]])
    cellsFew <- numCells[ numCells <= 3 ]
    excluThat <- vapply(contrasts, function(x) any(x %in% names(cellsFew)), logical(1))
    if (any(excluThat)) {
      message(glue::glue("Some cells too few, exclude from contrasts ({bind(which(excluThat))})."))
      contrasts <- contrasts[ !excluThat ]
    }
    if (is.null(x@params[[ glue::glue("{name}_contrast") ]]) || force) {
      res <- e(lapply(contrasts,
          function(con) {
            data <- Seurat::FindMarkers(object(x),
              ident.1 = con[1], ident.2 = con[2], min.pct = min.pct,
              group.by = group.by
            )
            data$gene <- rownames(data)
            data
          }))
      names(res) <- vapply(contrasts, function(x) paste0(x[1], "_vs_", x[2]), character(1))
      res <- dplyr::as_tibble(data.table::rbindlist(res, idcol = TRUE))
      res <- dplyr::rename(res, contrast = .id)
      res <- dplyr::filter(res, p_val_adj < cut.p, abs(avg_log2FC) > cut.fc)
      res <- .set_lab(res, sig(x), "DEGs of the contrasts")
      res <- setLegend(
        res, "细胞群差异表达基因附表 (其中 'contrast' 列为比较的两类细胞) (|log~2~(FC)| &gt; {cut.fc}, P-Adjust &lt; {cut.p})。"
      )
      x@params[[ glue::glue("{name}_contrast") ]] <- res
    } else {
      res <- x@params[[ glue::glue("{name}_contrast") ]]
    }
    p.volcano <- plot_volcano(
      res, "gene", use = "p_val_adj", use.fc = "avg_log2FC", fc = cut.fc, keep_cols = TRUE
    )
    p.volcano <- p.volcano + facet_wrap(~ contrast)
    p.volcano <- set_lab_legend(
      wrap(p.volcano),
      glue::glue("{x@sig} cell differential expression volcano plot"),
      glue::glue("细胞群差异表达基因火山图。")
    )
    x[[ glue::glue("{name}_volcano") ]] <- p.volcano
    init(snap(x)) <- TRUE
    snap <- vapply(contrasts, bind, character(1), co = " vs ")
    x <- snapAdd(
      x, "对细胞群差异分析 (依据 {group.by}, 分析 {less(snap)})，筛选差异表达基因。",
      step = paste0("diff_", name)
    )
    ## contrast intersection
    tops <- split(res, ~ contrast)
    use.gene <- "gene"
    fun_filter <- function(x) rm.no(x)
    tops <- lapply(tops,
      function(data){
        up <- dplyr::filter(data, avg_log2FC > 0)[[ use.gene ]]
        down <- dplyr::filter(data, avg_log2FC < 0)[[ use.gene ]]
        lst <- list(up = up, down = down)
        lapply(lst, fun_filter)
      })
    feature(x) <- as_feature(tops, x, analysis = "Seurat 细胞群差异表达分析")
    tops <- unlist(tops, recursive = FALSE)
    x[[ paste0(name, "_intersection") ]] <- tops
    p.sets_intersection <- new_upset(lst = tops, trunc = NULL)
    p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "contrasts-DEGs-intersection")
    p.sets_intersection <- setLegend(p.sets_intersection, "细胞群差异表达基因的 UpSet 交集图。")
    x[[ paste0("p.", name, "_intersection") ]] <- p.sets_intersection
    if (identical(parent.frame(2), .GlobalEnv)) {
      job_append_heading(x, heading = "细胞差异表达分析")
      job_append_method(
        x, paste0("diff_", name), TRUE, FALSE
      )
    }
    return(x)
  })

diff_group <- function(x, group.by, patterns) {
  group <- combn(as.character(ids(x, group.by)), 2, simplify = FALSE)
  try_contrast(group, patterns)
}

try_contrast <- function(combn, patterns) {
  lapply(combn,
    function(x) {
      res <- lapply(patterns,
        function(pattern) {
          matched <- grepl(pattern, x)
          if (any(matched))
            c(x[matched], x[!matched])
        })
      res <- lst_clear0(res)
      res[[1]]
    })
}

setMethod("map", signature = c(x = "job_seurat", ref = "marker_list"),
  function(x, ref, p.adjust = .05, log2fc = 1, use_all = TRUE, prop = .6, print = TRUE){
    if (x@step < 5L)
      stop("x@step < 5L")
    all_markers <- tables(x, 5)$all_markers
    all_markers <- filter(all_markers, p_val_adj < p.adjust, avg_log2FC > log2fc)
    all_markers <- select(all_markers, cluster, gene)
    all_markers <- split(all_markers$gene, all_markers$cluster)
    refs <- lapply(ref@level, function(n) ref@marker[[ n ]])
    stats <- lapply(all_markers,
      function(markers) {
        alls <- lapply(refs,
          function(ref) {
            all <- vapply(ref, FUN.VALUE = numeric(1),
              function(ref.) {
                res <- prop.table(table(ref. %in% markers))
                res <- unname(res[ names(res) == "TRUE" ])
                if (!length(res)) 0
                else res
              })
            tibble::tibble(cell_type = names(all), prop = unname(all))
          })
        names(alls) <- ref@level
        alls <- data.table::rbindlist(alls, fill = TRUE, idcol = TRUE)
        rename(alls, unique_level = .id)
      })
    map_prop <- data.table::rbindlist(stats, fill = TRUE, idcol = TRUE)
    map_prop <- rename(map_prop, cluster = .id)
    x@tables$map_prop <- map_prop
    map_prop <- mutate(map_prop, unique_level = as.integer(unique_level))
    map_prop <- arrange(map_prop, dplyr::desc(prop))
    if (use_all) {
      mapu <- distinct(map_prop, cell_type, .keep_all = TRUE)
      mapu <- rbind(mapu, map_prop)
      map_res <- distinct(mapu, cluster, .keep_all = TRUE)
    } else {
      map_res <- filter(map_prop, prop > !!prop)
      map_res <- distinct(map_res, cluster, .keep_all = TRUE)
    }
    if (print)
      print(map_res, n = 100)
    object(x)@meta.data$map_cell_type <- do.call(
      dplyr::recode,
      c(list(as.character(object(x)@meta.data$seurat_clusters)),
        nl(map_res$cluster, map_res$cell_type))
    )
    x@params$group.by <- "map_cell_type"
    p.map_marker_list <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = TRUE, pt.size = 0.5,
        group.by = "map_cell_type", cols = color_set()
        ))
    x@plots[[ 4 ]]$p.map_marker_list <- wrap(as_grob(p.map_marker_list))
    message(crayon::yellow(paste0("Use `plots(",
          as.character(substitute(x, env = parent.frame(1))),
          ", 4)$p.map_marker_list` to show result.")))
    return(x)
  })

setMethod("ids", signature = c(x = "job_seurat"),
  function(x, id = x@params$group.by, unique = TRUE){
    ids <- object(x)@meta.data[[ id ]]
    if (unique)
      ids <- unique(ids)
    ids
  })

plot_var2000 <- function(x) {
  top20 <- head(SeuratObject::VariableFeatures(x), 20)
  p.var2000 <- e(Seurat::VariableFeaturePlot(x))
  p.var2000 <- e(Seurat::LabelPoints(p.var2000, points = top20, repel = TRUE))
  p.var2000 <- wrap(p.var2000, 12, 8)
  p.var2000
}

plot_pca.seurat <- function(x) {
  p.pca_pcComponents <- e(Seurat::VizDimLoadings(
      x, dims = 1:2, reduction = "pca",
      col = ggsci::pal_npg()(2)[2], combine = FALSE
      ))
  p.pca_pcComponents <- lapply(p.pca_pcComponents,
    function(p) {
      p$layers[[1]]$aes_params$size <- 4
      p$layers[[1]]$aes_params$alpha <- .5
      p + theme_minimal()
    })
  p.pca_pcComponents <- patchwork::wrap_plots(p.pca_pcComponents, ncol = 2)
  p.pca_1v2 <- e(Seurat::DimPlot(x, reduction = "pca", cols = color_set()))
  p.pca_1v2 <- wrap(p.pca_1v2, 6, 5)
  p.pca_heatmap <- e(Seurat::DimHeatmap(
      x, dims = 1:10, cells = 500,
      balanced = TRUE, fast = FALSE, combine = FALSE
      ))
  n <- 0
  p.pca_heatmap <- lapply(p.pca_heatmap,
    function(p) {
      n <<- n + 1
      # p <- suppressMessages(p + scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF"))
      p <- p + ggtitle(paste0("PC_", n))
      p <- p + theme(plot.title = element_text(hjust = .5, face = "bold"))
      if (n < length(p.pca_heatmap))
        p <- p + theme(legend.position = "none")
      p
    })
  p.pca_heatmap <- patchwork::wrap_plots(p.pca_heatmap)
  p.pca_heatmap <- wrap(p.pca_heatmap, 12, 6)
  p.pca_rank <- pretty_elbowplot(e(Seurat::ElbowPlot(x)))
  p.pca_rank <- wrap(p.pca_rank, 4, 4)
  namel(
    p.pca_pcComponents,
    p.pca_1v2, p.pca_heatmap, p.pca_rank
  )
}

pretty_elbowplot <- function(plot) {
  plot$layers[[1]]$aes_params$size <- 5
  plot$layers[[1]]$aes_params$alpha <- .5
  plot$layers[[1]]$aes_params$colour <- "brown"
  plot
}

setMethod("getsub", signature = c(x = "job_seurat"),
  function(x, ...){
    object(x) <- e(SeuratObject:::subset.Seurat(object(x), ...))
    return(x)
  })

setMethod("active", signature = c(x = "job_seurat"),
  function(x, assay = "RNA"){
    object(x)@active.assay <- assay
    x
  })

plot_qc.seurat <- function(x) {
  require(patchwork)
  if (is(x, "job_seuratSp")) {
    nFeature <- "nFeature_Spatial"
    nCount <- "nCount_Spatial"
  } else if (is(x, "job_seurat")) {
    nFeature <- "nFeature_RNA"
    nCount <- "nCount_RNA"
  }
  x <- object(x)
  p.feature_count_mt <- e(Seurat::VlnPlot(x, features = c(nFeature, nCount,
        "percent.mt"), ncol = 3, pt.size = 0, alpha = .3, cols = color_set()))
  p.qcv1 <- e(Seurat::FeatureScatter(
      x, feature1 = nCount, feature2 = "percent.mt",
      cols = color_set()))
  p.qcv2 <- e(Seurat::FeatureScatter(
      x, feature1 = nCount, feature2 = nFeature,
      cols = color_set()))
  p.qc <- p.feature_count_mt / (p.qcv1 + p.qcv2)
  p.qc <- wrap(p.qc, 13, 10)
  p.qc
}

setMethod("vis", signature = c(x = "job_seurat"),
  function(x, group.by = x@params$group.by,
    pt.size = .7, mode = c("cell", "sample", "type"),
    palette = x$palette, reduction = "umap", type_pattern = "_[^_]+$",
    orig.ident = "orig.ident", ...)
  {
    mode <- match.arg(mode)
    if (is.null(palette)) {
      palette <- wgcna_colors()[-c(1, 4, 8, 9, 11)]
    }
    if (mode == "cell") {
    } else if (mode == "type") {
      groups <- unique(object(x)@meta.data[[ group.by ]])
      patterns <- gs(groups, type_pattern, "")
      palette <- pattern_gradientColor(patterns, groups, palette)
    } else if (mode == "sample") {
      x <- mutate(x, Cell_Sample = paste0(!!rlang::sym(group.by), "_", !!rlang::sym(orig.ident)))
      groups <- unique(object(x)@meta.data[[ group.by ]])
      patterns <- paste0("^", gs(make.names(groups), "[()]", "."), "_")
      Cell_Samples <- unique(object(x)@meta.data[[ "Cell_Sample" ]])
      palette <- pattern_gradientColor(
        patterns, Cell_Samples, palette
      )
      group.by <- "Cell_Sample"
    }
    p <- wrap(as_grob(e(Seurat::DimPlot(
            object(x), reduction = reduction, label = FALSE, pt.size = pt.size,
            group.by = group.by, cols = palette, ...
            ))), 7, 4)
    p <- setLegend(p, "为 {group.by} 的 {reduction} 聚类图。")
    .set_lab(p, sig(x), "The", gs(group.by, "_", "-"))
  })

setMethod("focus", signature = c(x = "job_seurat"),
  function(x, features, group.by = x@params$group.by, 
    sp = FALSE, name = "genes", return_type = c("job", "list"), ...)
  {
    if (is(features, "feature")) {
      features <- resolve_feature(features)
    }
    if (sp) {
      return(focus(.job_seuratSp(object = object(x)), features, ...))
    }
    layout <- calculate_layout(length(features))
    ncol <- layout[[ "cols" ]]
    p.vln <- e(Seurat::VlnPlot(
        object(x), features = features, group.by = group.by,
        pt.size = 0, alpha = .3, cols = color_set(), combine = FALSE
        ))
    p.vln <- lapply(p.vln, function(x) x + theme(legend.position = "none"))
    p.vln <- wrap_layout(patchwork::wrap_plots(p.vln, ncol = ncol), layout)
    p.vln <- .set_lab(
      p.vln, sig(x), paste("violing plot of expression level", name)
    )
    p.vln <- setLegend(p.vln, "基因 {less(features)} 表达水平的小提琴图。")
    p.dim <- wrap_layout(patchwork::wrap_plots(e(Seurat::FeaturePlot(
            object(x), features = features, combine = FALSE
            )), ncol = ncol), layout)
    p.dim <- .set_lab(
      p.dim, sig(x), paste("dimension plot of expression level", name)
    )
    p.dim <- setLegend(p.dim, "基因 {less(features)} 表达水平的 Dimension reduction plot.")
    return_type <- match.arg(return_type)
    if (return_type == "job") {
      x[[ paste0("focus_", name) ]] <- namel(p.vln, p.dim)
      return(x)
    } else {
      return(namel(p.vln, p.dim))
    }
  })

setMethod("map", signature = c(x = "job_seurat", ref = "character"),
  function(x, ref, slot = "scale.data", ...){
    p.heatmap <- wrap(as_grob(e(Seurat::DoHeatmap(
          object(x), features = ref, raster = TRUE, size = 3, slot = slot, ...
          ))))
    p.heatmap <- .set_lab(p.heatmap, sig(x), "heatmap show the reference genes")
    p.heatmap
  })

setMethod("map", signature = c(x = "job_seurat", ref = "sets_feature"),
  function(x, ref, datasets = NULL, names = NULL,
    group.by = x$group.by, sample = "orig.ident",
    pattern_suffix = "_[^_]+$", cell.props.filter = .1,
    feature.nlimit = 3, gene.props.filter = .1, gene.threshold = 0, 
    gene.which = NULL, excludes = "Monocyte", cut.cor = .3, cut.p = .05)
  {
    message("Correlation analysis for two feature sets within each cell group.")
    if (length(ref) != 2) {
      stop('length(ref) != 2.')
    }
    metadata <- object(x)@meta.data
    # filter by cells proportion above sample.
    if (!is.null(cell.props.filter)) {
      freq <- prop.table(table(metadata[[ group.by ]], metadata[[ sample ]]), 1)
      groups <- names(which(apply(freq, 1, function(x) entropy_score(freq = x)) > cell.props.filter))
      filter_out <- groups[ !groups %in% ids(x, group.by) ]
      if (length(filter_out)) {
        message(glue::glue("Filter out cell: \n{crayon::red(showStrings(filter_out))}"))
        x <- methodAdd(x, "\n\n以标准化熵筛选于样本中相对均衡分布的细胞 ({bind(groups)}) (cutoff: {cell.props.filter}) \n\n ({.entropy_score_method()})。\n\n", step = class(ref))
      } else {
        message(glue::glue("{crayon::yellow('No filter out cell.')}"))
      }
    } else {
      groups <- ids(x, group.by)
    }
    if (any(isUnknown <- grpl(groups, "unknown", TRUE))) {
      x <- methodAdd(x, "\n\n去除未知细胞 (unknown)。\n\n", step = class(ref))
      groups <- groups[ !isUnknown ]
    }
    if (!is.null(excludes)) {
      x <- methodAdd(x, "去除 {bind(excludes)} 细胞。\n\n", step = class(ref))
      groups <- groups[ !groups %in% excludes ]
    }
    fun_formal_name <- function(x) {
      gs(make.names(x), "[.]+", "_")
    }
    ## check validity
    lapply(ref, 
      function(fea) {
        if (is(fea, "feature_list")) {
          versus_cells <- .get_versus_cell(names(fea))
          group_names <- fun_formal_name(groups)
          if (all(!versus_cells %in% group_names)) {
            stop('all(!versus_cells %in% group_names), the name of the feature list illegal?')
          }
        }
      })
    ## set function of get genes for cells.
    fun_get_features <- function(fea, group) {
      if (is(fea, "feature_list")) {
        cells <- .get_versus_cell(names(fea), pattern_suffix)
        which <- which(cells == group)
        if (length(which) != 1) {
          NULL
        } else {
          unlist(fea[[ which ]])
        }
      } else if (is(fea, "feature_char")) {
        fea@.Data
      }
    }
    x <- methodAdd(
      x, "\n\n在上述的细胞群中，分析两组 features (即，{snap(ref[[1]])}，与 {snap(ref[[2]])})。\n\n",
      step = class(ref)
    )
    # set up features for each group (cell)
    sets <- sapply(groups, simplify = FALSE,
      function(group) {
        formalGroup <- fun_formal_name(group)
        from <- fun_get_features(ref[[1]], formalGroup)
        to <- fun_get_features(ref[[2]], formalGroup)
        if (!length(from) || !length(to)) {
          NULL
        } else {
          namel(from, to)
        }
      })
    sets <- lst_clear0(sets)
    ## filter the genes (if in `ref`) by expression levels (pct)
    if (!is.null(gene.props.filter)) {
      fun_get_stats <- function(feature) {
        feature <- feature@.Data
        if (all(!feature %in% rownames(object(x)))) {
          stop('all(!feature %in% rownames(object(x))), the feature is not "genes"?')
        }
        stats <- get_high_expressed(
          x, feature, threshold = gene.threshold,
          cutoff = gene.props.filter, group.by = group.by
        )
        split(stats$features, stats$group)
      }
      if (is.null(gene.which)) {
        seqs <- which(vapply(ref, is, logical(1), "feature_char"))
      } else {
        seqs <- gene.which
      }
      if (length(seqs)) {
        x <- methodAdd(x, "\n\n对于基因集，在各组细胞中，以阈值穿透率去除低表达的基因 (例如，去除总体表达为 0 的基因) (阈值：{gene.threshold}，穿透率 cutoff：{gene.props.filter}) ({.fast_penetrate_rate_method()})。\n\n", step = class(ref))
      }
      for (i in seqs) {
        lst_stats <- fun_get_stats(ref[[i]])
        n <- 0L
        sets <- lapply(sets, 
          function(from_and_to) {
            n <<- n + 1L
            isThat <- from_and_to[[i]] %in% lst_stats[[ names(sets)[n] ]]
            from_and_to[[i]] <- from_and_to[[i]][isThat]
            from_and_to
          })
      }
    }
    if (!is.null(feature.nlimit)) {
      sets <- lapply(sets, 
        function(lst) {
          if (length(lst$from) < feature.nlimit || length(lst$to) < feature.nlimit) {
            message('length(lst$from) < feature.nlimit || length(lst$to) < feature.nlimit, skip ...')
          } else {
            lst
          }
        })
      sets <- lst_clear0(sets)
      x <- methodAdd(x, "\n\n如果细胞群中，两组 features 均满足数量大于 {feature.nlimit}，则保留该细胞群，用于后续分析。\n\n", step = class(ref))
    }
    # set function of getting dataset (job_limma for 'cal_corp')
    fun_get_datasets <- function(feature, dataset, group) {
      if (is(dataset, "job_seurat")) {
        object <- asjob_limma(dataset, feature, cell_groups = group, group.by = group.by)
      } else if (is(dataset, "job_limma")) {
        if (is.null(.get_meta(dataset, group.by))) {
          stop('is.null(.get_meta(dataset, group.by)).')
        }
        object <- filter(
          dataset, type = "metadata", !!rlang::sym(group.by) == !!group, 
          add_snap = FALSE, force = TRUE
        )
      } else {
        stop("`datasets` should be a list of 'job_seurat' or 'job_limma'")
      }
      return(object)
    }
    # calculate correlation ...
    x <- methodAdd(
      x, "\n\n分别对各细胞群的两组 features 关联分析。\n\n", step = class(ref)
    )
    res_correlation <- sapply(names(sets), simplify = FALSE,
      function(group) {
        cli::cli_h2(glue::glue("Calculating: {group}"))
        from <- sets[[ group ]]$from
        to <- sets[[ group ]]$to
        lm_from <- fun_get_datasets(from, datasets[[1]], group)
        lm_to <- fun_get_datasets(to, datasets[[2]], group)
        message("Got datasets. Check cell number...")
        lapply(list(lm_from, lm_to), 
          function(x) {
            if (!nrow(x@params$normed_data$targets)[1]) {
              stop('!nrow(x@params$normed_data$targets)[1].')
            }
          })
        cp <- cal_corp(
          lm_from, lm_to, from, to, names = names, gname = FALSE,
          cut.cor = cut.cor, cut.p = cut.p, theme = paste(x@sig, fun_formal_name(group))
        )
        cp$res$hp <- wrap_scale(
          cp$res$hp, length(to), length(from), 4.5, 4 + max(nchar(to)) * .05
        )
        cp$res
      })
    s.com <- vapply(names(res_correlation),
      function(name) {
        n <- nrow(res_correlation[[name]]$sig.corp) %||% 0L
        glue::glue("{name} ({n})")
      }, character(1))
    x <- snapAdd(x, "设定 P 阈值 ({cut.p}) 与关联系数 ({cut.cor}) 阈值，获取各细胞群 (细胞的筛选算法见方法章节) 中的显著关联组。统计为: {bind(s.com)}。", step = class(ref))
    names(res_correlation) <- fun_formal_name(names(res_correlation))
    features <- sapply(
      names(res_correlation), simplify = FALSE, 
        function(name) {
          data <- res_correlation[[name]]$sig.corp
          list(from = unique(data[[1]]), to = unique(data[[2]]))
        }
    )
    analysis <- "细胞群 features 关联分析"
    x$.feature_cfrom <- as_feature(
      lapply(features, function(x) x$from),
      x, analysis = analysis
    )
    x$.feature_cto <- as_feature(
      lapply(features, function(x) x$to), 
      x, analysis = analysis
    )
    x$res_correlation <- res_correlation
    x$.map_snap <- "sets_feature"
    x$.map_heading <- analysis
    return(x)
  })


.get_versus_cell <- function(x, pattern_suffix = "_[^_]+$", 
  split = " - | vs |_vs_", unlist = TRUE)
{
  x <- strsplit(x, split)
  if (!is.null(pattern_suffix)) {
    x <- lapply(x, function(x) unique(s(x, pattern_suffix, "")))
    if (any(lengths(x) != 1)) {
      stop('any(lengths(x) != 1)')
    }
  }
  if (unlist) {
    unlist(x)
  } else {
    x
  }
}

entropy_score <- function(x, freq = NULL) {
  if (is.null(freq)) {
    freq <- prop.table(table(x))
  }
  -sum(freq * log(freq)) / log(length(freq))
}

.entropy_score_method <- function() {
  "给定离散随机变量 $X$，其取值为 ${x_1, x_2,...,x_K}$，对应概率分布为 $P(X = x_i) = p_i$，则 **归一化香农熵** 定义为：$H_{\\text{norm}}(X) = \\frac{ -\\sum_{i=1}^K p_i \\log p_i }{ \\log K }$，取值范围 $0 \\leq H_{\\text{norm}}(X) \\leq 1$"
}

setMethod("map", signature = c(x = "job_seurat", ref = "job_seurat"),
  function(x, ref, use.x, use.ref, name = "map_cell", asIdents = TRUE)
  {
    matched <- match(rownames(object(x)@meta.data), rownames(object(ref)@meta.data))
    object(x)@meta.data[[name]] <- as.character(object(ref)@meta.data[[use.ref]])[matched]
    object(x)@meta.data[[name]] <- ifelse(
      is.na(object(x)@meta.data[[name]]),
      as.character(object(x)@meta.data[[use.x]]),
      object(x)@meta.data[[name]]
    )
    object(x)@meta.data[[name]] <- factor(object(x)@meta.data[[name]],
      levels = sort(unique(object(x)@meta.data[[name]]))
    )
    if (asIdents) {
      if (identical(names(object(x)@active.ident), rownames(object(x)@meta.data))) {
        SeuratObject::Idents(object(x)) <- nl(rownames(object(x)@meta.data), object(x)@meta.data[[name]], FALSE)
      } else {
        stop("identical(names(object(x)@active.ident), rownames(object(x)@meta.data)) == FALSE")
      }
    }
    x$group.by <- name
    return(x)
  })

prepare_10x <- function(target, pattern, single = FALSE, col.gene = 1, check = TRUE) {
  if (!single) {
    files <- list.files(target, "\\.gz$", full.names = TRUE)
    files <- files[ grepl(pattern, files) ]
    dir <- paste0(target, "/", get_realname(files)[1])
    dir.create(dir, FALSE)
    lapply(files,
      function(x)
        file.rename(x, paste0(dir, "/", gs(basename(x), ".*_", "")))
    )
    expected <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    alls <- list.files(dir, full.names = TRUE)
    allnames <- basename(alls)
    if (!all(expected %in% allnames)) {
      message("Got all files: ", paste0(allnames, collapse = ", "))
      notGet <- expected[ !expected %in% allnames ]
      Get <- expected[ expected %in% allnames ]
      message("Got: ", paste0(Get, collapse = ", "), "\n",
        "Not got: ", crayon::yellow(paste0(notGet, collapse = ", ")))
      isUnExpected <- !allnames %in% expected
      if (length(which(isUnExpected)) == 1 & length(notGet) == 1) {
        file.rename(alls[ isUnExpected ], paste0(dir, "/", notGet))
        message(crayon::red("Rename: "), alls[ isUnExpected ], " as ", paste0(dir, "/", notGet))
      }
    }
    dir
  } else {
    if (!missing(pattern)) {
      stop("`pattern` not supported while single == TRUE")
    }
    path <- dirname(target)
    name <- get_realname(target)
    dir <- paste0(path, "/", name)
    if (file.exists(dir)) {
      unlink(dir, TRUE, TRUE)
    }
    data <- ftibble(target)
    if (length(col.gene) > 1) {
      message("Use first element of `col.gene`, and removed the others.")
      data <- dplyr::relocate(data, dplyr::all_of(col.gene))
      data <- dplyr::select(data, -dplyr::all_of(seq_along(col.gene)[-1]))
    }
    print(data)
    if (check) {
      isCon <- usethis::ui_yeah("Continue?")
      if (!isCon) {
        stop("...")
      }  
    }
    data[[ 1 ]] <- gs(data[[ 1 ]], "\\.[0-9]*$", "")
    data <- dplyr::distinct(data, !!rlang::sym(colnames(data)[1]), .keep_all = TRUE)
    ## as Matrix
    features <- data[[ 1 ]]
    data <- as.matrix(data[, -1])
    rownames(data) <- features
    require(Matrix)
    data <- as(data, "sparseMatrix")
    DropletUtils::write10xCounts(dir, data)
    return(dir)
  }
}

setMethod("skel", signature = c(x = "job_seurat"),
  function(x, sig = "sr"){
    code <- c('',
      'sr <- job_seurat("")',
      'sr <- step1(sr)',
      'sr@plots$step1$p.qc',
      'sr <- step2(sr, 0, 7500, 35)',
      'sr@plots$step2$p.pca_rank',
      'sr <- step3(sr, 1:15, 1.2)',
      'sr@step <- 4L',
      'sr@plots$step3$p.umap',
      'sr <- step5(sr, 5)',
      'sr <- step6(sr, "Kidney")',
      'sr@plots$step6$p.map_scsa',
      "clear(sr)"
    )
    code <- gs(code, "\\bsr\\b", sig)
    writeLines(code)
  })

setGeneric("asjob_seurat", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_seurat"))

setMethod("mutate", signature = c(x = "job_seurat"),
  function(x, ...){
    object(x)@meta.data <- dplyr::mutate(object(x)@meta.data, ...)
    return(x)
  })

setMethod("cal_corp", signature = c(x = "job_seurat", y = "NULL"),
  function(x, y, from, to, names = NULL)
  {
    data <- object(x)[[ object(x)@active.assay ]]@data
    data <- data[ rownames(data) %in% unique(c(from, to)), ]
    anno <- data.frame(symbol = rownames(data))
    data <- data.frame(Matrix::as.matrix(data))
    .cal_corp.elist(data, anno, use = "symbol", from, to, names)
  })

setMethod("meta", signature = c(x = "job_seurat"),
  function(x){
    as_tibble(object(x)@meta.data)
  })

setMethod("feature", signature = c(x = "job_seurat"),
  function(x, mode = .all_features(x), ...){
    if (missing(mode)) {
      mode <- ".feature"
    } else {
      if (!grpl(mode, "^\\.feature")) {
        mode <- paste0(".feature_", mode)
      }
      mode <- match.arg(mode)
    }
    feas <- x[[ mode ]]
    if (!is(feas, "feature")) {
      feas <- as_feature(feas, x, ...)
    }
    feas
  })

setMethod("set_remote", signature = c(x = "job_seurat"),
  function(x, wd = glue::glue("~/seurat_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

.all_features <- function(x) {
  if (!is(x, "job")) {
    stop('!is(x, "job").')
  }
  alls <- names(x@params)
  alls[ grpl(alls, "^\\.feature") ]
}

applyGptcelltype <- function(input, tissuename, model = c("gpt-3.5-turbo", "gpt-4"),
  topgenenumber = 10)
{
  # https://platform.openai.com/api-keys
  OPENAI_API_KEY <- getOption("OPENAI_API_KEY")
  if (is.null(OPENAI_API_KEY)) {
    API.flag <- 0L
  } else {
    API.flag <- 1L
  }
  model <- match.arg(model)
  if (is(input, "list")) {
    input <- sapply(input, paste, collapse = ',')
  } else {
    if (!is.data.frame(input)) {
      stop("`input` must be either data.frame or list.")
    }
    if (!all(c("cluster", "gene") %in% colnames(input))) {
      stop("The data.frame format of `input` must contain columns: cluster, gene.")
    }
    input <- tapply(input$gene, list(input$cluster),
      function(i) paste0(i[1:topgenenumber], collapse = ','))
  }
  if (!API.flag){
   stop('Identify cell types of ', tissuename,
     ' cells using the following markers separately for each\n row.',
     ' Only provide the cell type name. Do not show numbers before the name.\n',
     ' Some can be a mixture of multiple cell types.\n',
     paste0(names(input), ':', unlist(input), collapse = "\n"))
  } else {
    cutnum <- ceiling(length(input) / 30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(seq_along(input), cutnum))
    } else {
      cid <- rep(1, length(input))
    }
    allres <- sapply(1:cutnum, simplify = FALSE,
      function(i) {
        id <- which(cid == i)
        content <- paste0(
          'Identify cell types of ', tissuename, ', ',
          ' using the following markers separately for each row.',
          ' Only provide the cell type name.',
          ' Do not show numbers before the name.\n',
          ' Some can be a mixture of multiple cell types.\n\n',
          paste(input[id], collapse = '\n')
        )
        flag <- 0L
        while (!flag) {
          k <- openai::create_chat_completion(model = model,
            messages = list(list(role = "user", content = "How are you?")),
            openai_api_key = OPENAI_API_KEY
          )
          stop_debug(k)
          res <- strsplit(k$choices[, 'message.content'], '\n')[[1]]
          if (length(res) == length(id))
            flag <- 1
        }
        names(res) <- names(input)[id]
        res
      })
    gsub(',$', '', unlist(allres))
  }
}

identify.mouseMacroPhe <- function(x, use = "scsa_cell",
  cell.name = "Macrophage",
  markers = x@tables$step5$all_markers_no_filter, top.ref = 10, ...)
{
  message("This function only support for organism of 'mouse'.\n",
    "Make sure the dataset has been subset, only retain the Macrophage cells.")
  ref <- get_data.nmt2015()
  Macrophage_M0 <- dplyr::filter(ref@object$sig.m0, FC.M1_vs_M0 < 1)$gene
  Macrophage_M1 <- dplyr::filter(ref@object$sig.m1, FC.M1_vs_M0 > 1)$gene
  Macrophage_M2 <- dplyr::filter(ref@object$sig.m2, FC.M2_vs_M0 > 1)$gene
  ref.markers <- namel(Macrophage_M0, Macrophage_M1, Macrophage_M2)
  ref.markers <- lapply(ref.markers, function(i) i[ i %in% rownames(object(x)@assays$SCT$scale.data) ])
  if (!is.null(top.ref)) {
    ref.markers <- lapply(ref.markers, head, n = top.ref)
  }
  ref.markers <- as_df.lst(ref.markers, "cell", "markers")
  lst <- scsa_annotation(x, "All", ref.markers, org = "M", res.col = "macrophage_phenotypes",
    onlyUseRefMarkers = TRUE, ...)
  x <- lst$x
  ref.markers <- .set_lab(as_tibble(ref.markers), sig(x), "the markers for Macrophage phenotypes annotation")
  x$p.macrophage <- lst$p.map_scsa
  x$t.macrophage_ref <- ref.markers
  x$p.macrophage_hp_ref <- map(x, ref.markers$markers)
  ref@object <- NULL
  .add_internal_job(ref)
  return(x)
}

scsa_annotation <- function(
  x, tissue, ref.markers = NULL, onlyUseRefMarkers = !is.null(ref.markers),
  filter.p = 0.01, filter.fc = .5, reset = NULL,
  org = c("Human", "Mouse"),
  cmd = pg("scsa"), db = pg("scsa_db"), res.col = "scsa_annotation",
  cache = "scsa", ...)
{
  if (grpl(tissue, "(?<!\\\\)\\s", perl = TRUE)) {
    tissue <- gs(tissue, "\\s", "\\\\ ")
  }
  org <- match.arg(org)
  prop <- as.list(prop.table(table(grpl(rownames(object(x)), "^[A-Z][a-z]"))))
  if (!is.null(prop[[ "TRUE" ]])) {
    if (prop[[ "TRUE" ]] > .5 && org == "Human") {
      isThat <- usethis::ui_yeah("Match org of 'Human', but most of genes is Mouse like symbol, continue?")
      if (!isThat) {
        stop("...")
      }
    }  
  }
  ## cache
  hash <- e(
    digest::digest(
      list(
        x@sig, x@tables$step5$all_markers_no_filter,
        tissue, ref.markers, onlyUseRefMarkers, filter.p, filter.fc, org
      ), "md5"
    )
  )
  message(glue::glue("Got hash: {hash}"))
  cache <- paste0(cache, "_", hash)
  dir.create(cache, FALSE)
  mkfile <- function(name, fileext = NULL) {
    file.path(cache, paste0(name, fileext))
  }
  ## results file
  marker_file <- mkfile("marker_file", fileext = ".csv")
  message(glue::glue("Marker file: {marker_file}"))
  result_file <- mkfile("result_file")
  message(glue::glue("Results file: {result_file}"))
  if (!file.exists(result_file)) {
    if (!is.null(ref.markers)) {
      .check_columns(ref.markers, c("cell", "markers"))
      ref.markers <- dplyr::relocate(ref.markers, cell, markers)
      ref.markers_file <- mkfile("ref.markers_file", fileext = ".tsv")
      message(glue::glue("Reference markers: {ref.markers_file}"))
      write_tsv(ref.markers, ref.markers_file, col.names = FALSE)
      x@params$ref.markers_file <- ref.markers_file
      ref.markers.cmd <- paste0(" -M ", ref.markers_file)
    } else {
      ref.markers.cmd <- ""
    }
    all_markers <- dplyr::rename(x@tables$step5$all_markers_no_filter, avg_logFC = avg_log2FC)
    all_markers <- dplyr::select(all_markers, -rownames)
    all_markers <- dplyr::mutate(all_markers, gene = gs(gene, "\\.[0-9]*", ""))
    data.table::fwrite(all_markers, marker_file)
    cli::cli_alert_info(cmd)
    cdRun(cmd,
      " -s seurat", " -i ", marker_file,
      " -k ", tissue, " -d ", db, " ", ref.markers.cmd,
      " -p ", filter.p, " -f ", filter.fc,
      " -E -g ", org, " -m txt",
      if (onlyUseRefMarkers) " -N " else NULL,
      " -o ", result_file, " > /tmp/log_scsa.txt"
    )
  }
  scsa_res <- dplyr::rename_all(ftibble(result_file), make.names)
  scsa_res <- scsa_res_all <- dplyr::arrange(scsa_res, Cluster, dplyr::desc(Z.score))
  ## add into seurat object
  scsa_res <- dplyr::distinct(scsa_res, Cluster, .keep_all = TRUE)
  if (!is.null(reset)) {
    if (is.null(names(reset))) {
      stop('is.null(names(reset)).')
    }
    if (is.character(reset)) {
      reset <- as.list(reset)
    }
    scsa_res <- dplyr::mutate(
      scsa_res, Cell.Type = dplyr::recode(
        Cell.Type, !!!reset, .default = Cell.Type
      )
    )
  }
  clusters <- object(x)@meta.data$seurat_clusters
  cell_types <- scsa_res$Cell.Type[match(clusters, scsa_res$Cluster)]
  cell_types <- ifelse(is.na(cell_types), "Unknown", cell_types)
  if (!is.null(reset)) {
    cell_types <- dplyr::recode(cell_types, !!!reset, .default = cell_types)
  }
  object(x)@meta.data[[ res.col ]] <- factor(cell_types)
  ## plot
  p.map_scsa <- e(Seurat::DimPlot(
      object(x), reduction = "umap", label = FALSE, pt.size = .7,
      group.by = res.col, cols = color_set()
      ))
  p.map_scsa <- wrap(as_grob(p.map_scsa), 7, 4)
  p.map_scsa <- .set_lab(p.map_scsa, sig(x), "SCSA", "Cell type annotation")
  p.map_scsa <- setLegend(p.map_scsa, "为 SCSA 细胞注释结果的 UMAP 图。")
  .add_internal_job(.job(method = "`SCSA` (python) used for cell type annotation",
      cite = "[@ScsaACellTyCaoY2020]"))
  namel(x, p.map_scsa, res.col, scsa_res_all, scsa_res)
}

matchCellMarkers <- function(lst.markers, ref, least = 2) {
  lst <- vapply(lst.markers, FUN.VALUE = logical(1),
    function(x) {
      x <- as.list(table(ref %in% x))
      x[[ "TRUE" ]] >= least
    }
  )
  lst
}

parse_GPTfeedback <- function(feedback) {
  res <- gs(feedback, "^[0-9. ]+|\\.$|\\*|\\(.*\\)\\s*", "")
  res <- strsplit(res, "\\. ")
  cells <- vapply(res, function(x) strsplit(x[1], "; ")[[1]][1], character(1))
  cells <- gs(cells, "^\\s*|\\s*$", "")
  markers <- lapply(res, function(x) strsplit(x[2], "; "))
  data <- tibble::tibble(cells = cells, markers = markers)
  markers <- markers[ order(cells) ]
  markers <- unique(gs(unlist(markers), "\\s", ""))
  namel(markers, cells, data)
}

prepare_GPTmessage_for_celltypes <- function(tissue, 
  marker_list, n = 10, toClipboard = TRUE, save = !toClipboard,
  save_name = "query_message_for_gpt.txt")
{
  message <- glue::glue("Identify cell types of {tissue} cells using the following markers separately for each row. Only provide the cell type name (for each row). Show numbers before the name. Some can be a mixture of multiple cell types (separated by '; ', end with '. '). Then, provide at least 1 and at most 5 markers  that distinguish the cell type from other cells (separated by '; ').  (e.g., 1. X Cell; Y Cell. Marker1; Marker2; Marker3)")
  message("`marker_list` should be table output from seurat (column: cluster, gene)")
  marker_list <- split(marker_list$gene, marker_list$cluster)
  ncluster <- length(marker_list)
  if (!is.null(n)) {
    marker_list <- lapply(marker_list, head, n = n)
  }
  markers <- vapply(marker_list, paste0, character(1), collapse = ",")
  markers <- paste0(seq_along(markers), ". ", markers)
  markers <- paste(markers, collapse = "\n")
  query <- paste0(message, "\n", markers)
  if (toClipboard) {
    gett(query)
  } else {
    writeLines(query, save_name)
  }
  return(list(query = query, message = message, ncluster = ncluster))
}

.modify_gpt_marker_and_cell <- function(x, all_markers, 
  query, res, keep_markers)
{
  ## is High discriminative genes?
  all_markers$cells <- res$cells[match(all_markers$cluster, 0:(query$ncluster - 1))]
  all_markers <- dplyr::filter(all_markers, gene %in% res$markers)
  numInMarkers <- vapply(res$markers,
    function(marker) {
      cells <- unique(dplyr::filter(all_markers, gene == !!marker)$cells)
      length(cells)
    }, integer(1))
  isPrimMarkers <- numInMarkers <= 2L
  notDiscMarkers <- res$markers[!isPrimMarkers]
  message(
    glue::glue(
      "Total not discriminative markers {length(notDiscMarkers)} / {length(res$markers)}:
      \t{bind(notDiscMarkers)}"
    )
  )
  cell_markers <- reframe_col(res$data, "markers", unlist)
  cell_markers$num <- numInMarkers[ match(cell_markers$markers, res$markers) ]
  cell_markers <- dplyr::arrange(cell_markers, num)
  ## keep 2-3 markers, for each cells.
  cell_markers <- split_lapply_rbind(
    cell_markers, ~ cells, head, n = keep_markers
  )
  x$hp_cell_markers <- cell_markers
  excluOverMarkers <- unique(cell_markers$markers)
  message(
    glue::glue(
      "Markers excluded overlaps: {length(excluOverMarkers)}
      (Not discriminative {length(which(excluOverMarkers %in% notDiscMarkers))})"
    )
  )
  gptCells <- object(x)@meta.data[["ChatGPT_cell"]]
  allCells <- unique(cell_markers$cells)
  cells_filter <- c()
  for (i in allCells) {
    dat <- dplyr::filter(cell_markers, cells == i)
    if (all(dat$num > 3)) {
      gptCells <- ifelse(gptCells == i, "Unknown", gptCells)
      cells_filter <- c(cells_filter, i)
      next
    } else if (all(dat$num > 1)) {
      dat_overlap <- dplyr::filter(cell_markers, markers %in% dat$markers)
      maybeParent <- dat_overlap$cells[ which.min(nchar(dat_overlap$cells)) ]
      if (maybeParent != i && grpl(i, maybeParent, TRUE)) {
        message(glue::glue("Re assign cells: {i} -> {maybeParent}"))
        gptCells <- ifelse(gptCells == i, maybeParent, gptCells)
        cells_filter <- c(cells_filter, i)
      }
    }
  }
  if (length(cells_filter)) {
    cell_markers <- dplyr::filter(cell_markers, !cells %in% cells_filter)
  }
  cell_markers <- dplyr::arrange(cell_markers, cells)
  excluOverMarkers <- unique(cell_markers$markers)
  message(glue::glue("Post Markers excluded overlaps: {length(excluOverMarkers)}"))
  object(x)@meta.data[["ChatGPT_cell"]] <- gptCells
  x$gpt_excluOverMarkers <- excluOverMarkers
  return(x)
}

.plot_marker_heatmap <- function(x, markers, group.by,
  show = NULL, extra.after = NULL,
  order.by = c("smart", "raw"), 
  max = 2, soft = TRUE, notShow = NULL, scale = FALSE)
{
  if ((is(markers, "list") || is(markers, "feature")) && !is.null(max)) {
    lst <- markers
    markers <- unlist(markers)
  } else {
    max <- NULL
  }
  if (!is.null(show)) {
    if (!is(show, "list")) {
      show <- list(Extra = show)
    }
    if (!is.null(extra.after)) {
      lst <- append(lst, show, after = extra.after)
      markers <- unlist(lst)
    } else {
      markers <- c(markers, show[[1]])
      lst <- c(lst, show)
    }
  }
  if (!is.null(notShow)) {
    markers <- markers[ !markers %in% notShow ]
    lst <- lapply(lst, 
      function(x) {
        x[ !x %in% notShow ]
      })
  }
  avgExpr <- Seurat::AverageExpression(
    object(x), features = markers,
    assays = object(x)@active.assay,
    group.by = group.by, slot = "data"
  )
  avgExpr <- data.frame(avgExpr[[ 1 ]], check.names = FALSE)
  colnames(avgExpr) <- sort(unique(object(x)@meta.data[[ group.by ]]))
  avgExpr <- dplyr::rename(as_tibble(avgExpr), Gene = rownames)
  avgExpr <- tidyr::pivot_longer(avgExpr, -Gene, names_to = "Cell", values_to = "Expression")
  avgExpr <- dplyr::mutate(
    avgExpr, Expression = log2(Expression + 1)
  )
  if (!is.null(max)) {
    exprCutoff <- fivenum(avgExpr$Expression)[4]
    marker_keep <- lapply(lst,
      function(gene) {
        if (soft && length(gene) <= max) {
          return(gene)
        }
        data <- dplyr::filter(avgExpr, Gene %in% !!gene)
        data <- dplyr::group_by(data, Gene)
        data <- dplyr::summarise(
          data, score = .calculate_deviation_score(Expression), 
          maxExpr = max(Expression)
        )
        data <- dplyr::arrange(data, dplyr::desc(score))
        if (soft) {
          testData <- dplyr::filter(data, maxExpr >= exprCutoff)
          if (nrow(testData) > max) {
            data <- testData
          } else {
            ## if too less, re-search...
            nkeep <- max
            for (i in seq(max, nrow(data))) {
              if (any(data$maxExpr[seq_len(i)] >= exprCutoff)) {
                nkeep <- i
                break
              }
            }
            max <- nkeep
          }
        }
        head(data$Gene, n = max)
      })
    avgExpr <- dplyr::filter(avgExpr, Gene %in% unlist(marker_keep))
  }
  order.by <- match.arg(order.by)
  if (order.by == "raw") {
    avgExpr <- dplyr::mutate(
      avgExpr, Gene = factor(Gene, levels = unique(markers))
    )
  } else if (order.by == "smart") {
    allMarkers <- unique(avgExpr$Gene)
    cutoff <- mean(range(avgExpr$Expression))
    dat <- lapply(split(avgExpr, ~ Cell),
      function(x) {
        x <- dplyr::arrange(x, dplyr::desc(Expression))
        dplyr::mutate(x,
          isOutlier = Expression >= find_outliers(Expression),
          isHigh = Expression >= cutoff
        )
      })
    dat <- dplyr::arrange(
      rbind_list(dat), dplyr::desc(isOutlier)
    )
    glevels <- unique(dat$Gene)
    avgExpr <- dplyr::mutate(avgExpr, Gene = factor(Gene, levels = glevels))
  }
  if (scale) {
    avgExpr <- dplyr::mutate(avgExpr, Expression = scale(Expression))
    fun_palette <- fun_color(
      values = avgExpr$Expression, category = "div", rev = TRUE
    )
  } else {
    fun_palette <- fun_color(
      values = avgExpr$Expression, category = "seq", rev = FALSE
    )
  }
  p.markers <- tidyHeatmap::heatmap(
    avgExpr, Gene, Cell, Expression, cluster_columns = FALSE, 
    palette_value = fun_palette,
    cluster_rows = FALSE
  )
  p.markers@arguments <- list()
  p.markers
}

find_outliers <- function(x, try_gap = TRUE) {
  iqr <- find_outliers_iqr(x)
  if (length(iqr) || !try_gap) {
    min(iqr)
  } else {
    gap <- find_outliers_gap(x)
    if (length(gap)) {
      min(gap)
    } else {
      Inf
    }
  }
}

find_outliers_gap <- function(x) {
  sorted_x <- sort(x)
  diffs <- diff(sorted_x)
  max_gap_index <- which.max(diffs)
  if (max_gap_index < length(sorted_x)) {
    outliers <- sorted_x[(max_gap_index + 1):length(sorted_x)]
  } else {
    outliers <- numeric(0)
  }
  outliers
}

find_outliers_iqr <- function(x, coef = 1.5) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  upper_bound <- q[2] + coef * iqr
  outliers <- x[x > upper_bound]
  outliers
}

.calculate_deviation_score <- function(x) {
  max <- max(x)
  skewness <- e1071::skewness(x)
  dispersion <- max / median(x)
  gap <- (max - sort(x, decreasing = TRUE)[2]) / max
  skewness + dispersion + gap
}

plot_cells_proportion <- function(metadata, sample = "orig.ident", cell = "ChatGPT_cell") {
  ntypes <- length(unique(metadata[[ cell ]]))
  data <- data.frame(
    prop.table(table(droplevels(metadata[[ cell ]]), metadata[[ sample ]]), 1)
  )
  p.props <- ggplot(data, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(stat = "identity", width = .7, size = .5) +
    coord_flip() +
    labs(x = "Cells", y = "Ratio", fill = "Sample") +
    scale_fill_manual(values = color_set()) +
    rstyle("theme")
  wrap(p.props, 7, .4 * ntypes)
}

get_high_expressed <- function(x, features, threshold = 0, cutoff = .3,
  group.by = x$group.by, data = object(x)$RNA$data)
{
  if (!is(x, "job_seurat")) {
    stop('!is(x, "job_seurat").')
  }
  if (!is(data, "dgCMatrix")) {
    stop('!is(data, "dgCMatrix").')
  }
  data <- data[ rownames(data) %in% features, ]
  if (!nrow(data)) {
    stop('!nrow(data).')
  }
  groups <- object(x)@meta.data[[ group.by ]]
  if (is.null(groups)) {
    stop('is.null(groups): object(x)@meta.data[[ group.by ]]')
  }
  stats <- fast_penetrate_rate(data, groups, threshold)
  stats$features <- rownames(data)[stats$row]
  stats <- tibble::as_tibble(stats)
  if (!is.null(cutoff)) {
    stats <- dplyr::filter(stats, ratio > !!cutoff)
  }
  stats
}

fast_penetrate_rate <- function(sp_mat, groups, threshold = 0) {
  if (!is(sp_mat, "dgCMatrix")) {
    stop('!is(sp_mat, "dgCMatrix").')
  }
  require(data.table)
  # Convert groups into numerical indexes
  groups <- factor(groups)
  groupLevels <- levels(groups)
  group_idx <- as.integer(groups)
  n_groups <- max(group_idx)
  
  # Extract non-zero coordinates and values
  smry <- Matrix::summary(sp_mat)
  dt <- data.table(
    row = smry$i,
    col = smry$j,
    val = smry$x,
    group = group_idx[smry$j]
  )
  
  # Filter expression values that exceed the threshold
  dt_pass <- dt[val > threshold, .(count = .N), by = .(row, group)]
  
  # Calculate the total number of cells in each group
  group_size <- data.table(group = group_idx)[, .(total = .N), by = group]
  
  # Merge calculation of penetration rate
  result <- dt_pass[
    group_size, 
    on = "group", 
    .(row, group, ratio = count / total)
  ]
  result <- result[!is.na(ratio), ]
  result$group <- groupLevels[ result$group ]
  result
}

.fast_penetrate_rate_method <- function(feature = "基因", 
  type = "细胞", level = "表达")
{
  glue::glue(
    "设某{{{feature}}} $g$ 在{{{type}}}群体 $C$ 中的{{{level}}}值集合为 ${e_c | c \\in C}$，给定阈值 $\\tau$，则 **阈值穿透率** 定义为：$\\text{Penetration}(g, C, \\tau) = \\frac{ \\sum_{c \\in C} \\mathbf{1}_{\\{e_c > \\tau\\}} }{ |C| } \\times 100\\%$ ($\\mathbf{1}_{{e_c > \\tau}}$ 是指示函数，当 $e_c > \\tau$ 时为 1，否则为 0)",
    .open = "{{{", .close = "}}}"
  )
}

pgc <- pattern_gradientColor <- function(pattern, names,
  colors = ggsci::pal_rickandmorty()(10), ...)
{
  if (length(pattern) > length(colors)) {
    message("`colors` provided not enough.")
    colors <- color_set()
  }
  names <- as.character(unique(names))
  colors <- colors[ seq_along(pattern) ]
  n <- 0L
  palette <- lapply(pattern,
    function(pt) {
      n <<- n + 1L
      x <- sort(names[ grpl(names, pt, ...) ])
      colors <- colorRampPalette(c("white", colors[n]))(length(x) + 2)
      colors <- colors[-c(1, length(x) - 1)]
      names(colors) <- x
      return(colors)
    })
  unlist(palette)
}
