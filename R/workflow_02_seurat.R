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
    info = c("Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html")
    ))

job_seurat <- function(dir, project = make.names(date()),
  min.cells = 3, min.features = 200, ...)
{
  data <- Seurat::Read10X(dir)
  object <- Seurat::CreateSeuratObject(counts = data, project = project,
    min.cells = min.cells, min.features = min.features, ...)
  .job_seurat(object = object)
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    callNextMethod()
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
    x <- callNextMethod()
    step_message("QC and selecting cells for further analysis.",
      "This do:",
      "Generate column in `object(x)@meta.data`; ",
      "Plots in `x@plots[[ 1 ]]`.",
      "Check the plots then decided the ",
      crayon::red("min nFeature_RNA"), " and ",
      crayon::red("max nFeature_RNA"), ", ",
      "as well as ", crayon::red("percent.mt"), "."
    )
    object(x)[[ "percent.mt" ]] <- Seurat::PercentageFeatureSet(
      object(x), pattern = "^MT-"
    )
    require(patchwork)
    p.qc <- plot_qc.seurat(object(x))
    x@plots[[ 1 ]] <- list(p.qc = p.qc)
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000){
    x <- callNextMethod(x)
    step_message("This contains several execution: Subset the data,
      Normalization, Feature selection, Scale the data, Linear dimensional
      reduction, Select dimensionality.
      red{{`min.features`}} and red{{`max.features`}} were needed for subset.
      Then `object(x)` were performed with:
      `Seurat::NormalizeData`;
      `Seurat::FindVariableFeatures`;
      `Seurat::SCTransform`;
      `Seurat::RunPCA`.
      All plots were in `x@plots[[ 2 ]]`
      NOTE: inspect the plots and red{{Determined dims}} for downstream analysis. "
    )
    if (missing(min.features) | missing(max.features))
      stop("missing(min.features) | missing(max.features)")
    message("\n`Run SeuratObject::subset.Seurat`")
    object(x) <- SeuratObject:::subset.Seurat(
      object(x), subset = nFeature_RNA > min.features &
        nFeature_RNA < max.features & percent.mt < max.percent.mt
    )
    message("\nRun `Seurat::NormalizeData`")
    object(x) <- Seurat::NormalizeData(
      object(x), normalization.method = "LogNormalize", scale.factor = 10000
    )
    message("\nRun: `Seurat::FindVariableFeatures`")
    object(x) <- Seurat::FindVariableFeatures(
      object(x), selection.method = "vst", nfeatures = 2000
    )
    p.var2000 <- plot_var2000(object(x))
    message("\nRun: `Seurat::SCTransform`")
    object(x) <- Seurat::SCTransform(
      object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T,
    )
    message("\nRun: `Seurat::RunPCA`")
    object(x) <- Seurat::RunPCA(
      object(x), features = SeuratObject::VariableFeatures(object(x))
    )
    message("\nPlot PCA ...")
    x@plots[[ 2 ]] <- c(namel(p.var2000), plot_pca.seurat(object(x)))
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat"),
  function(x, dims, resolution){
    x <- callNextMethod(x)
    step_message("This contains several execution: Cluster the cells, UMAP
      reduction, Cluster biomarker (Differential analysis). Object were
      performed with:
      `Seurat::FindNeighbors`;
      `Seurat::FindClusters` grey{{(clusters can be accessed by `SeuratObject::Idents`)}};
      `Seurat::RunUMAP`;
      `Seurat::FindAllMarkers`.
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
    if (ncol(object(x)) > 3000 & resolution < 1.2) {
      stop("ncol(object(x)) > 3000 but resolution < 1.2. ",
        "Please consider higher `resolution`. ")
    }
    object(x) <- Seurat::FindNeighbors(object(x), dims = dims)
    object(x) <- Seurat::FindClusters(object(x), resolution = resolution)
    object(x) <- Seurat::RunUMAP(object(x), dims = dims)
    p.umap <- Seurat::DimPlot(object(x), reduction = "umap", cols = color_set())
    markers <- as_tibble(
      Seurat::FindAllMarkers(object(x), min.pct = 0.25, logfc.threshold = 0.25)
    )
    tops <- dplyr::filter(markers, p_val_adj < .05)
    tops <- slice_max(group_by(markers, cluster), avg_log2FC, n = 10)
    p.toph <- Seurat::DoHeatmap(object(x), features = tops$gene, raster = F)
    p.toph <- wrap(p.toph, 14, 12)
    x@tables[[ 3 ]] <- list(all_markers = markers)
    x@plots[[ 3 ]] <- namel(p.umap, p.toph)
    return(x)
  })

setMethod("step4", signature = c(x = "job_seurat"),
  function(x, ref = celldex::HumanPrimaryCellAtlasData()){
    x <- callNextMethod(x)
    step_message("Use `SingleR` and `celldex` to annotate cell types.
      By default, red{{`celldex::HumanPrimaryCellAtlasData`}} was used
      as red{{`ref`}} dataset. This annotation would generate red{{'SingleR_cell'}}
      column in `object(x)@meta.data`. Plots were generated in `x@plots[[ 4 ]]`;
      tables in `x@tables[[ 4 ]]`.
      ")
      ref <- celldex::HumanPrimaryCellAtlasData()
      clusters <- object(x)@meta.data$seurat_clusters
      anno_SingleR <- SingleR::SingleR(object(x)@assays$SCT@scale.data,
        ref = ref, labels = ref$label.fine, clusters = clusters
      )
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
      p.map_SingleR <- Seurat::DimPlot(
        object(x), reduction = "umap", label = TRUE, pt.size = 0.5,
        group.by = "SingleR_cell", cols = color_set()
      )
      p.map_SingleR <- p.map_SingleR + theme(legend.position = "none")
      p.map_SingleR <- wrap(p.map_SingleR, 10, 7)
      x@tables[[ 4 ]] <- namel(anno_SingleR)
      x@plots[[ 4 ]] <- namel(p.score_SingleR, p.map_SingleR)
      return(x)
  })

plot_var2000 <- function(x) {
  top20 <- head(SeuratObject::VariableFeatures(x), 20)
  p.var2000 <- Seurat::VariableFeaturePlot(x)
  p.var2000 <- Seurat::LabelPoints(p.var2000, points = top20, repel = TRUE)
  p.var2000 <- wrap(p.var2000, 12, 8)
  p.var2000
}

plot_pca.seurat <- function(x) {
  p.pca_pcComponents <- Seurat::VizDimLoadings(
    x, dims = 1:2, reduction = "pca",
    col = ggsci::pal_npg()(2)[2], combine = F
  )
  p.pca_pcComponents <- lapply(p.pca_pcComponents,
    function(p) {
      p$layers[[1]]$aes_params$size <- 4
      p$layers[[1]]$aes_params$alpha <- .5
      p + theme_minimal()
    })
  p.pca_pcComponents <- patchwork::wrap_plots(p.pca_pcComponents, ncol = 2)
  p.pca_1v2 <- Seurat::DimPlot(x, reduction = "pca", cols = ggsci::pal_npg()(10))
  p.pca_heatmap <- Seurat::DimHeatmap(
    x, dims = 1:10, cells = 500,
    balanced = TRUE, fast = F, combine = F
  )
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
  p.pca_rank <- Seurat::ElbowPlot(x)
  p.pca_rank$layers[[1]]$aes_params$size <- 5
  p.pca_rank$layers[[1]]$aes_params$alpha <- .5
  p.pca_rank$layers[[1]]$aes_params$colour <- "brown"
  namel(
    p.pca_pcComponents,
    p.pca_1v2, p.pca_heatmap, p.pca_rank
  )
}

plot_qc.seurat <- function(x) {
  p.feature_count_mt <- Seurat::VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA",
      "percent.mt"), ncol = 3, pt.size = 0, alpha = .3, cols = c("lightyellow"))
  p.qcv1 <- Seurat::FeatureScatter(
    x, feature1 = "nCount_RNA", feature2 = "percent.mt",
    cols = c("lightblue")
  )
  p.qcv2 <- Seurat::FeatureScatter(
    x, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
    cols = c("lightblue")
  )
  p.qc <- p.feature_count_mt / (p.qcv1 + p.qcv2)
  p.qc <- wrap(p.qc, 14, 12)
  p.qc
}
