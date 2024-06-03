# ==========================================================================
# workflow of seuratSp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seuratSp <- setClass("job_seuratSp", 
  contains = c("job_seurat"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("Tutorial: https://satijalab.org/seurat/articles/spatial_vignette.html",
      "\nAdditional: https://ludvigla.github.io/STUtility_web_site/index.html"),
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019]",
    method = "R package `Seurat` used for spatial scRNA-seq analysis",
    tag = "scrna:strna"
    ))

job_seuratSp <- function(dir)
{
  obj <- e(Seurat::Load10X_Spatial(dir))
  .job_seuratSp(object = obj)
}

setMethod("step0", signature = c(x = "job_seuratSp"),
  function(x){
    step_message("Prepare your data with function `job_seuratSp`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_seuratSp"),
  function(x){
    step_message("Quality control (QC).")
    x <- callNextMethod(x)
    p.nCount_spatial <- e(Seurat::SpatialFeaturePlot(object(x), features = "nCount_Spatial")) +
      theme(legend.position = "right")
    p.nCount_spatial <- wrap(as_grob(p.nCount_spatial), 5.5, 4)
    x@plots[[ 1 ]] <- c(x@plots[[ 1 ]], namel(p.nCount_spatial))
    return(x)
  })

setMethod("step2", signature = c(x = "job_seuratSp"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000,
    use = "nFeature_Spatial")
  {
    x <- callNextMethod(x, min.features, max.features, max.percent.mt, nfeatures, use)
    return(x)
  })

setMethod("step3", signature = c(x = "job_seuratSp"),
  function(x, dims, resolution){
    x <- callNextMethod(x, dims, resolution)
    palette <- .setPaletteForSpatialJob(x, "seurat_clusters")
    p.spatial_cluster <- e(Seurat::SpatialDimPlot(object(x), label = TRUE, label.size = 3, cols = palette)) +
      guides(fill = guide_legend(override.aes = list(size = 3)))
    p.spatial_cluster <- wrap(x@plots[[ 3 ]]$p.umap@data + p.spatial_cluster, 13, 5.5)
    x@plots[[ 3 ]] <- namel(p.spatial_cluster)
    return(x)
  })

setMethod("vis", signature = c(x = "job_seuratSp"),
  function(x, group.by = x@params$group.by, pt.size = .7){
    p.umap <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = F, pt.size = pt.size,
        group.by = group.by, cols = color_set()
        ))
    palette <- .setPaletteForSpatialJob(x, group.by)
    p.spatial <- e(Seurat::SpatialDimPlot(object(x),
        group.by = group.by, label = F, label.size = 3, cols = palette)) +
      guides(fill = guide_legend(override.aes = list(size = 3)))
    p <- wrap(as_grob(p.umap + p.spatial), 13, 5.5)
    .set_lab(p, sig(x), "The", gs(group.by, "_", "-"))
  })

setMethod("step4", signature = c(x = "job_seuratSp"),
  function(x, ref = celldex::HumanPrimaryCellAtlasData()){
    x <- callNextMethod(x, ref)
    return(x)
  })

setMethod("step5", signature = c(x = "job_seuratSp"),
  function(x, workers = 3, spatial = F){
    x <- callNextMethod(x, workers)
    if (spatial) {
      object(x) <- e(Seurat::FindSpatiallyVariableFeatures(object(x),
          assay = "SCT", selection.method = "moransi",
          features = Seurat::VariableFeatures(object(x)), verbose = T))
      moI.top_features <- SpatiallyVariableFeatures_workaround(object(x))
      x@params$moI.top_features <- moI.top_features
    }
    return(x)
  })

.setPaletteForSpatialJob <- function(x, group.by = x@params$group.by) {
  groups <- x@object@meta.data[[ group.by ]]
  if (is.factor(groups)) {
    levels <- levels(groups)
  } else {
    levels <- sort(unique(groups))
  }
  nl(levels, color_set()[1:length(levels)], F)
}

SpatiallyVariableFeatures_workaround <- function(object, assay = "SCT", selection.method = "moransi") {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[data[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  # Return row names of the sorted data frame
  rownames(sorted_data)
}

