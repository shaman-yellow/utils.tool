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
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019; @ScsaACellTyCaoY2020]",
    method = "Seurat used for scRNA-seq processing; SCSA used for cell type annotation"
    ))

job_seurat <- function(dir, project = get_filename(sub("/$", "", dir)),
  min.cells = 3, min.features = 200, ...)
{
  data <- e(Seurat::Read10X(dir))
  object <- e(Seurat::CreateSeuratObject(counts = data, project = project,
      min.cells = min.cells, min.features = min.features, ...))
  .job_seurat(object = object)
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
    p.qc <- plot_qc.seurat(x)
    p.qc <- .set_lab(p.qc, sig(x), "Quality", "Control")
    x@plots[[ 1 ]] <- list(p.qc = p.qc)
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000,
    use = "nFeature_RNA")
  {
    step_message("This contains several execution: Subset the data,
      Normalization, Feature selection, Scale the data, Linear dimensional
      reduction, Select dimensionality.
      red{{`min.features`}} and red{{`max.features`}} were needed for subset.
      Then `object(x)` were performed with:
      `Seurat::SCTransform`;
      `Seurat::RunPCA`.
      All plots were in `x@plots[[ 2 ]]`
      NOTE: inspect the plots and red{{Determined dims}} for downstream analysis. "
    )
    if (missing(min.features) | missing(max.features))
      stop("missing(min.features) | missing(max.features)")
    object(x) <- e(SeuratObject:::subset.Seurat(
        object(x), subset = !!rlang::sym(use) > min.features &
          !!rlang::sym(use) < max.features & percent.mt < max.percent.mt
        ))
    object(x) <- e(Seurat::SCTransform(
        object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T,
        assay = SeuratObject::DefaultAssay(object(x))
        ))
    object(x) <- e(Seurat::RunPCA(
        object(x), features = SeuratObject::VariableFeatures(object(x))
        ))
    x@plots[[ 2 ]] <- plot_pca.seurat(object(x))
    x@plots[[ 2 ]]$p.pca_rank <- .set_lab(
      x@plots[[ 2 ]]$p.pca_rank, sig(x), "Ranking of principle components")
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat"),
  function(x, dims, resolution, force = F){
    step_message("This contains several execution: Cluster the cells, UMAP
      reduction, Cluster biomarker (Differential analysis). Object were
      performed with:
      `Seurat::FindNeighbors`;
      `Seurat::FindClusters` grey{{(clusters can be accessed by `SeuratObject::Idents`)}};
      yellow{{The parameter `resolution` between 0.4-1.2 typically returns good results for
      single-cell datasets of around 3000 cells (ncol). Optimal resolution
      often increases for larger datasets.}}
      grey{{Consider use Seurat::FindMarkers Seurat::FindAllMarkers(..., test.use = 'roc') for
      customize differential analysis.}}
      grey{{Consider use `Seurat::VlnPlot` or `Seurat::FeaturePlot`}} for
      customize mapping of features."
    )
    if (missing(dims) | missing(resolution))
      stop("missing(dims) | missing(resolution)")
    if (!force) {
      if (ncol(object(x)) > 3000 & resolution < 1.2) {
        stop("ncol(object(x)) > 3000 but resolution < 1.2. ",
          "Please consider higher `resolution`. ")
      }
    }
    # if (isNamespaceLoaded("future")) {
      # e(pkgload::unload("future"))
    # }
    object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims))
    object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution))
    object(x) <- e(Seurat::RunUMAP(object(x), dims = dims))
    p.umap <- e(Seurat::DimPlot(object(x), reduction = "umap", cols = color_set()))
    p.umap <- wrap(p.umap, 6, 5)
    p.umap <- .set_lab(p.umap, sig(x), "UMAP", "Clustering")
    x@plots[[ 3 ]] <- namel(p.umap)
    return(x)
  })

setMethod("step4", signature = c(x = "job_seurat"),
  function(x, use = "SingleR", use.level = c("label.main", "label.fine"),
    ref = celldex::HumanPrimaryCellAtlasData())
  {
    use <- "SingleR"
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
          object(x), reduction = "umap", label = F, pt.size = .7,
          group.by = "SingleR_cell", cols = color_set()
          ))
      p.map_SingleR <- p.map_SingleR
      p.map_SingleR <- wrap(p.map_SingleR, 10, 7)
      x@tables[[ 4 ]] <- namel(anno_SingleR)
      x@plots[[ 4 ]] <- namel(p.score_SingleR, p.map_SingleR)
      x@params$group.by <- "SingleR_cell"
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_seurat"),
  function(x, workers = 3){
    step_message("Find all Marders for Cell Cluster.")
    fun <- function(x) {
      markers <- as_tibble(
        e(Seurat::FindAllMarkers(object(x), min.pct = 0.25,
            logfc.threshold = 0.25, only.pos = T))
      )
    }
    if (!is.null(workers)) {
      markers <- parallel(x, fun, workers)
    } else {
      markers <- fun(x)
    }
    all_markers_no_filter <- markers
    markers <- dplyr::filter(markers, p_val_adj < .05)
    tops <- slice_max(group_by(markers, cluster), avg_log2FC, n = 10)
    if (F) {
      p.toph <- e(Seurat::DoHeatmap(object(x), features = tops$gene, raster = T))
      p.toph <- wrap(p.toph, 14, 12)
      x@plots[[ 5 ]] <- namel(p.toph)
    }
    x@tables[[ 5 ]] <- list(all_markers = markers, all_markers_no_filter = all_markers_no_filter)
    return(x)
  })

setMethod("step6", signature = c(x = "job_seurat"),
  function(x, tissue, ref.markers = NULL, filter.p = 0.01, filter.fc = .5,
    cmd = "python3 ~/SCSA/SCSA.py", db = "~/SCSA/whole_v2.db")
  {
    step_message("Use SCSA to annotate cell types (<https://github.com/bioinfo-ibms-pumc/SCSA>).")
    if (!is.null(ref.markers)) {
      .check_columns(ref.markers, c("cell", "markers"))
      ref.markers <- dplyr::relocate(ref.markers, cell, markers)
      ref.markers_file <- tempfile("ref.markers_file", fileext = ".tsv")
      write_tsv(ref.markers, ref.markers_file, col.names = F)
      x@params$ref.markers_file <- ref.markers_file
      ref.markers.cmd <- paste0(" -M ", ref.markers_file)
    } else {
      ref.markers.cmd <- ""
    }
    marker_file <- tempfile("marker_file", fileext = ".csv")
    result_file <- tempfile("result_file")
    all_markers <- dplyr::rename(x@tables$step5$all_markers_no_filter, avg_logFC = avg_log2FC)
    all_markers <- dplyr::select(all_markers, -rownames)
    all_markers <- dplyr::mutate(all_markers, gene = gs(gene, "\\.[0-9]*", ""))
    data.table::fwrite(all_markers, marker_file)
    cli::cli_alert_info(cmd)
    cdRun(cmd,
      " -s seurat", " -i ", marker_file,
      " -k ", tissue, " -d ", db, " ", ref.markers.cmd,
      " -p ", filter.p, " -f ", filter.fc,
      " -E -g Human -m txt",
      " -o ", result_file, " > /tmp/log_scsa.txt"
    )
    x@params$marker_file <- marker_file
    x@params$scsa_res_file <- result_file
    x@params$scsa_log <- readLines("/tmp/log_scsa.txt")
    scsa_res <- dplyr::rename_all(ftibble(result_file), make.names)
    scsa_res <- dplyr::arrange(scsa_res, Cluster, dplyr::desc(Z.score))
    x@tables[[ 6 ]] <- namel(scsa_res)
    ## add into seurat object
    scsa_res <- dplyr::distinct(scsa_res, Cluster, .keep_all = T)
    clusters <- object(x)@meta.data$seurat_clusters
    cell_types <- scsa_res$Cell.Type[match(clusters, scsa_res$Cluster)]
    cell_types <- ifelse(is.na(cell_types), "Unknown", cell_types)
    object(x)@meta.data$scsa_cell <- factor(cell_types,
      levels = c(unique(scsa_res$Cell.Type), "Unknown"))
    x@params$group.by <- "scsa_cell"
    ## plot
    p.map_scsa <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = F, pt.size = .7,
        group.by = "scsa_cell", cols = color_set()
        ))
    p.map_scsa <- wrap(as_grob(p.map_scsa), 7, 4)
    p.map_scsa <- .set_lab(p.map_scsa, sig(x), "SCSA", "Cell type annotation")
    x@plots[[ 6 ]] <- list(p.map_scsa = p.map_scsa)
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

setMethod("diff", signature = c(x = "job_seurat"),
  function(x, group.by, contrasts){
    if (is.data.frame(contrasts)) {
      contrasts <- apply(contrasts, 1, c, simplify = F)
    }
    if (is.null(x@params$contrasts)) {
      res <- e(lapply(contrasts,
          function(con) {
            data <- Seurat::FindMarkers(object(x),
              ident.1 = con[1], ident.2 = con[2],
              group.by = group.by,
              test.use = "negbinom"
            )
            data$gene <- rownames(data)
            data
          }))
      names(res) <- vapply(contrasts, function(x) paste0(x[1], "_vs_", x[2]), character(1))
      res <- dplyr::as_tibble(data.table::rbindlist(res, idcol = T))
      res <- dplyr::rename(res, contrast = .id)
      res <- dplyr::filter(res, p_val_adj < .05)
      x@params$contrasts <- res
    } else {
      res <- x@params$contrasts
    }
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
    tops <- unlist(tops, recursive = F)
    x$diff_sets_intersection <- tops
    p.sets_intersection <- new_upset(lst = tops, trunc = NULL)
    p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "contrasts-DEGs-intersection")
    x$p.diff_sets_intersection <- p.sets_intersection
    return(x)
  })

diff_group <- function(x, group.by, patterns) {
  group <- combn(as.character(ids(x, group.by)), 2, simplify = F)
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
  function(x, ref, p.adjust = .05, log2fc = 1, use_all = T, prop = .6, print = T){
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
        alls <- data.table::rbindlist(alls, fill = T, idcol = T)
        rename(alls, unique_level = .id)
      })
    map_prop <- data.table::rbindlist(stats, fill = T, idcol = T)
    map_prop <- rename(map_prop, cluster = .id)
    x@tables$map_prop <- map_prop
    map_prop <- mutate(map_prop, unique_level = as.integer(unique_level))
    map_prop <- arrange(map_prop, dplyr::desc(prop))
    if (use_all) {
      mapu <- distinct(map_prop, cell_type, .keep_all = T)
      mapu <- rbind(mapu, map_prop)
      map_res <- distinct(mapu, cluster, .keep_all = T)
    } else {
      map_res <- filter(map_prop, prop > !!prop)
      map_res <- distinct(map_res, cluster, .keep_all = T)
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
  function(x, id = x@params$group.by, unique = T){
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
      col = ggsci::pal_npg()(2)[2], combine = F
      ))
  p.pca_pcComponents <- lapply(p.pca_pcComponents,
    function(p) {
      p$layers[[1]]$aes_params$size <- 4
      p$layers[[1]]$aes_params$alpha <- .5
      p + theme_minimal()
    })
  p.pca_pcComponents <- patchwork::wrap_plots(p.pca_pcComponents, ncol = 2)
  p.pca_1v2 <- e(Seurat::DimPlot(x, reduction = "pca", cols = ggsci::pal_npg()(10)))
  p.pca_1v2 <- wrap(p.pca_1v2, 6, 5)
  p.pca_heatmap <- e(Seurat::DimHeatmap(
      x, dims = 1:10, cells = 500,
      balanced = TRUE, fast = F, combine = F
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
  p.pca_rank <- e(Seurat::ElbowPlot(x))
  p.pca_rank$layers[[1]]$aes_params$size <- 5
  p.pca_rank$layers[[1]]$aes_params$alpha <- .5
  p.pca_rank$layers[[1]]$aes_params$colour <- "brown"
  p.pca_rank <- wrap(p.pca_rank, 4, 4)
  namel(
    p.pca_pcComponents,
    p.pca_1v2, p.pca_heatmap, p.pca_rank
  )
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
  function(x, group.by = x@params$group.by, pt.size = .7, palette = x$palette){
    if (is.null(palette)) {
      palette <- color_set()
    }
    p <- wrap(e(Seurat::DimPlot(
          object(x), reduction = "umap", label = F, pt.size = pt.size,
          group.by = group.by, cols = palette
          )), 7, 4)
    .set_lab(p, sig(x), "The", gs(group.by, "_", "-"))
  })

setMethod("focus", signature = c(x = "job_seurat"),
  function(x, features, group.by = x@params$group.by){
    p.vln <- wrap(e(Seurat::VlnPlot(
        object(x), features = features, group.by = group.by,
        pt.size = 0, alpha = .3, cols = color_set()
        )))
    p.vln <- .set_lab(p.vln, sig(x), "violing plot of expression level of the genes")
    p.dim <- wrap(e(Seurat::FeaturePlot(
        object(x), features = features
        )))
    p.dim <- .set_lab(p.dim, sig(x), "dimension plot of expression level of the genes")
    namel(p.vln, p.dim)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_seurat"),
  function(x, ref, use.x, use.ref, name = "cell_mapped", asIdents = T){
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
        SeuratObject::Idents(object(x)) <- nl(rownames(object(x)@meta.data), object(x)@meta.data[[name]], F)
      } else {
        stop("identical(names(object(x)@active.ident), rownames(object(x)@meta.data)) == F")
      }
    }
    return(x)
  })

prepare_10x <- function(target, pattern, single = F) {
  if (!single) {
    files <- list.files(target, "\\.gz$", full.names = T)
    files <- files[ grepl(pattern, files) ]
    dir <- paste0(target, "/", get_realname(files)[1])
    dir.create(dir, F)
    lapply(files,
      function(x)
        file.rename(x, paste0(dir, "/", gs(get_filename(x), ".*_", "")))
    )
    dir
  } else {
    path <- get_path(target)
    name <- get_realname(target)
    dir <- paste0(path, "/", name)
    if (file.exists(dir)) {
      unlink(dir, T, T)
    }
    data <- ftibble(target)
    data[[ 1 ]] <- gs(data[[ 1 ]], "\\.[0-9]*$", "")
    data <- dplyr::distinct(data, !!rlang::sym(colnames(data)[1]), .keep_all = T)
    ## as Matrix
    features <- data[[ 1 ]]
    data <- as.matrix(data[, -1])
    rownames(data) <- features
    require(Matrix)
    mtx <- as(data, "Matrix")
    DropletUtils::write10xCounts(dir, mtx)
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

setGeneric("asjob_seurat", 
  function(x, ...) standardGeneric("asjob_seurat"))

setMethod("mutate", signature = c(DF_object = "job_seurat"),
  function(DF_object, ...){
    object(DF_object)@meta.data <- dplyr::mutate(object(DF_object)@meta.data, ...)
    return(DF_object)
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
