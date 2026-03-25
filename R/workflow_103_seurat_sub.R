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
    if (!all(colnames(object(x)) %in% metadata[[cellName]])) {
      stop('!all(colnames(object(x)) %in% metadata[[cellName]]).')
    }
    if (!any(metadata[[ group.by ]] %in% group)) {
      stop('!any(metadata[[ group.by ]] %in% group).')
    }
    exclude <- NULL
    if (TRUE) {
      metaEx <- dplyr::filter(metadata, !!rlang::sym(group.by) %in% !!group)
      freq <- table(metaEx$orig.ident)
      if (any(freq < 2)) {
        exclude <- names(freq)[ freq < 2 ]
        methodAdd_onExit("x", " (去除了对应细胞 ({bind(group)}) 数量过少 (&lt; 2) 的样本，防止 Seurat 运行错误，这些样本为: {bind(exclude)}) ")
        message(glue::glue("Pre exclude sample: {bind(exclude)}"))
        message("Prevent only one target cell from remaining in the sample, causing dgcMatrix drop to become numeric and resulting in errors.")
        metadata <- dplyr::filter(metadata, !orig.ident %in% exclude)
        object(x) <- object(x)[, !object(x)@meta.data$orig.ident %in% exclude ]
      }
    }
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
    orderSub <- match(colnames(object(x)), metadata[[cellName]])
    object(x)@meta.data <- cbind(
      object(x)@meta.data, metadata[ orderSub,  ]
    )
    x <- getsub(x, !!rlang::sym(group.by) %in% !!group)
    isDropped <- vapply(
      object(x)@assays$RNA@layers, FUN.VALUE = logical(1),
      function(x) {
        is.numeric(x)
      }
    )
    if (any(isDropped)) {
      sampleNames <- s(names(object(x)@assays$RNA@layers), "^counts\\.", "")
      stop(glue::glue("Please remove the Dropped assay layers: {bind(sampleNames[isDropped])}"))
      # object(x)@assays$RNA@layers <- object(x)@assays$RNA@layers[!isDropped]
      # object(x)@assays$RNA@default <- length(object(x)@assays$RNA@layers)
      # lgMap <- object(x)@assays$RNA@cells@.Data
      # mata <- object(x)@meta.data <- dplyr::filter(object(x)@meta.data, !orig.ident %in% !!sampleNames)
      # object(x)@assays$RNA@cells@.Data <- lgMap[
      #   rownames(lgMap) %in% rownames(meta), colnames(lgMap) %in% meta$orig.ident
      #   ]
      # object(x)@active.ident <- as.factor(object(x)@meta.data$orig.ident)
    }
    validObject(object(x))
    x <- .job_seurat_sub(object = object(x))
    x <- methodAdd(x, "提取 {bind(group)} 对其二次降维聚类 (重新对其按照完整工作流整合不同样本) 以分析其亚群。")
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

as_markers <- function(cell_markers, snap = NULL, df = NULL, ref = "pmid")
{
  if (!is.null(df)) {
    colnames(df) <- c("cell", "markers")
    cell_markers <- df
  } else {
    type <- vapply(cell_markers, class, character(1))
    if (all(type == "character")) {
      cell_markers <- as_df.lst(cell_markers, "cell", "markers")
    } else if (all(type == "list")) {
      lst <- lapply(cell_markers, 
        function(x) {
          setNames(
            tibble::tibble(markers = x[["markers"]], ref = bind(x[[ref]])),
            c("markers", ref)
          )
        })
      cell_markers <- dplyr::bind_rows(lst, .id = "cell")
      refs <- bind(unique(unlist(strsplit(cell_markers[[ref]], ", "))))
      snap(cell_markers) <- glue::glue("{toupper(ref)}: {refs}")
    }
  }
  if (!is.null(snap)) {
    snap(cell_markers) <- snap
  }
  cell_markers
}

