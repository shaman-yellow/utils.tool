# ==========================================================================
# workflow of seurat_sub
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat_sub <- setClass("job_seurat_sub", 
  contains = c("job_seurat"),
  prototype = prototype(
    analysis = "Seurat 细胞亚群分析"
    ))

setGeneric("asjob_seurat_sub", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_seurat_sub"))

setMethod("asjob_seurat_sub", signature = c(x = "job_seurat"),
  function(x, ref, group, group.by = "scsa_cell", cellName = "cell",
    get_after = "seurat_clusters")
  {
    if (missing(ref)) {
      ref <- as_tibble(object(x)@meta.data, idcol = cellName)
    }
    if (is(ref, "job_seurat")) {
      if (is.null(object(ref))) {
        message(glue::glue("is.null(object(ref)), try get metadata from `ref$final_metadata`"))
        metadata <- ref$final_metadata
      } else {
        metadata <- as_tibble(object(ref)@meta.data, idcol = cellName)
      }
    } else {
      metadata <- ref
    }
    if (!is(metadata, "data.frame")) {
      stop('!is(metadata, "data.frame").')
    }
    if (x@step != 1) {
      stop('x@step != 1.')
    }
    .check_columns(metadata, c(cellName, group.by), "metadata")
    if (!is.null(get_after)) {
      metadata <- dplyr::select(
        metadata, !!rlang::sym(cellName), 
        !!rlang::sym(group.by), !!rlang::sym(get_after):last_col()
      )
    } else {
      metadata <- dplyr::select(
        metadata, !!rlang::sym(cellName), !!rlang::sym(group.by)
      )
    }
    if (!all(colnames(object(x)) %in% metadata[[cellName]])) {
      stop('!all(colnames(object(x)) %in% metadata[[cellName]]).')
    }
    if (!any(metadata[[ group.by ]] %in% group)) {
      stop('!any(metadata[[ group.by ]] %in% group).')
    }
    orderSub <- match(colnames(object(x)), metadata[[cellName]])
    object(x)@meta.data <- cbind(
      object(x)@meta.data, metadata[ orderSub,  ]
    )
    snapAdd_onExit("x", "提取 {bind(group)} 对其二次降维聚类以分析其亚群。")
    x <- getsub(x, !!rlang::sym(group.by) %in% !!group)
    x <- .job_seurat_sub(object = object(x))
    x$group.by <- group.by
    return(x)
  })

setMethod("step0", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Prepare your data with function `job_seurat_sub`.")
  })

setMethod("step1", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Render as job_seurat5n.")
    x <- .job_seurat5n(x)
    return(x)
  })

# setMethod("step2", signature = c(x = "job_seurat_sub"),
#   function(x){
#     if (x$sub_from == "job_seurat5n") {
#       step_message("Render as `job_seurat5n`.")
#       x <- .job_seurat5n(x)
#       x$JoinLayers <- TRUE
#     } else {
#       step_message("Do nothing.")
#     }
#     return(x)
#   })

# setMethod("step2", signature = c(x = "job_seurat_sub"),
#   function(x, reset = FALSE, sct = FALSE) {
#     step_message("Reset variable features.")
#     if (reset || sct) {
#       if (sct) {
#         object(x) <- e(Seurat::SCTransform(
#             object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE,
#             assay = SeuratObject::DefaultAssay(object(x))
#             ))
#       } else {
#         object(x) <- e(Seurat::NormalizeData(object(x)))
#         object(x) <- e(Seurat::FindVariableFeatures(object(x)))
#         object(x) <- e(Seurat::ScaleData(object(x)))
#       }
#       object(x) <- e(Seurat::RunPCA(
#           object(x), 
#           features = SeuratObject::VariableFeatures(object(x)), reduction.name = "pca"
#           ))
#       p.pca_rank <- e(Seurat::ElbowPlot(object(x), 30))
#       p.pca_rank <- wrap(pretty_elbowplot(p.pca_rank), 4, 4)
#       p.pca_rank <- .set_lab(p.pca_rank, sig(x), "Standard deviations of PCs")
#       p.pca_rank <- setLegend(p.pca_rank, "为主成分 (PC) 的 Standard deviations。")
#       x@plots[[ 2 ]] <- namel(p.pca_rank)
#     }
#     return(x)
#   })

# setMethod("step3", signature = c(x = "job_seurat_sub"),
#   function(x, dims = 1:15, resolution = 1.2, force = TRUE, ...)
#   {
#     x <- callNextMethod(x, dims, resolution, force = force, ...)
#     return(x)
#   })

as_markers <- function(cell_markers, snap = NULL, df = NULL) {
  if (!is.null(df)) {
    colnames(df) <- c("cell", "markers")
    cell_markers <- df
  } else {
    cell_markers <- as_df.lst(cell_markers, "cell", "markers")
  }
  if (!is.null(snap)) {
    snap(cell_markers) <- snap
  }
  cell_markers
}

