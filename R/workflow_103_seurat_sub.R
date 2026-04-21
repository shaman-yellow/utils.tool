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
    get_after = "seurat_clusters", refine = NULL, cut.refine = 0)
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
      object(x)@meta.data, metadata[ orderSub, !colnames(metadata) %in% colnames(object(x)@meta.data)]
    )
    x <- getsub(x, !!rlang::sym(group.by) %in% !!group)
    methodAdd_onExit("x", "提取 {bind(group)} 对其二次降维聚类 (重新对其按照完整工作流整合不同样本) 以分析其亚群。")
    isDropped <- vapply(
      object(x)@assays$RNA@layers, FUN.VALUE = logical(1),
      function(x) {
        is.numeric(x)
      }
    )
    if (any(isDropped)) {
      sampleNames <- s(names(object(x)@assays$RNA@layers), "^counts\\.", "")
      stop(glue::glue("Please remove the Dropped assay layers: {bind(sampleNames[isDropped])}"))
    }
    if (!is.null(refine)) {
      if (!is(refine, "refine")) {
        stop('!is(refine, "refine"), the refine is not valid?')
      }
      if (length(unlist(refine$cells)) != ncol(object(x))) {
        stop(
          'length(unlist(refine$cells)) != ncol(object(x)), not match cells number?'
        )
      }
      message(
        glue::glue("Before refine cells, dim: {bind(dim(object(x)))}")
      )
      ncell_pre <- ncol(object(x))
      cells_keep <- unlist(
        dplyr::filter(tibble::as_tibble(refine), purity > cut.refine)$cells
      )
      object(x) <- e(SeuratObject:::subset.Seurat(object(x), cells = cells_keep))
      message(
        glue::glue("After refine cells, dim: {bind(dim(object(x)))}")
      )
      ncell_aft <- ncol(object(x))
      methodAdd_onExit("x", "{snap(refine)}在对细胞提纯之前，子集细胞总数量为 {ncell_pre}，在提纯之后，细胞总数量为 {ncell_aft}。")
    }
    validObject(object(x))
    x <- .job_seurat_sub(object = object(x))
    x$group.by <- group.by
    return(x)
  })

# object(x)@assays$RNA@layers <- object(x)@assays$RNA@layers[!isDropped]
# object(x)@assays$RNA@default <- length(object(x)@assays$RNA@layers)
# lgMap <- object(x)@assays$RNA@cells@.Data
# mata <- object(x)@meta.data <- dplyr::filter(object(x)@meta.data, !orig.ident %in% !!sampleNames)
# object(x)@assays$RNA@cells@.Data <- lgMap[
#   rownames(lgMap) %in% rownames(meta), colnames(lgMap) %in% meta$orig.ident
#   ]
# object(x)@active.ident <- as.factor(object(x)@meta.data$orig.ident)

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

setMethod("refine", signature = c(x = "job_seurat"),
  function(x, pos, neg, group.by = "seurat_clusters")
  {
    if (x@step < 3L) {
      stop('x@step < 3L.')
    }
    object(x) <- Seurat::AddModuleScore(
      object(x),
      features = list(pos),
      name = "tmp_pos"
    )
    object(x) <- Seurat::AddModuleScore(
      object(x),
      features = list(neg),
      name = "tmp_neg"
    )
    raw <- meta <- object(x)@meta.data
    meta$purity <- meta$tmp_pos1 - meta$tmp_neg1
    meta <- data.table::as.data.table(meta)
    stat <- meta[ , .(
      tmp_pos1 = mean(tmp_pos1, na.rm = TRUE),
      tmp_neg1 = mean(tmp_neg1, na.rm = TRUE),
      purity = mean(purity, na.rm = TRUE),
      n_cells = .N), by = group.by
      ]
    stat <- tibble::as_tibble(stat)
    cells <- split(rownames(raw), raw[[ group.by ]])
    stat$cells <- lapply(stat[[ group.by ]],
      function(group) {
        cells[[ as.character(group) ]]
      })
    stat$n_cells <- lengths(stat$cells)
    stat <- structure(stat, class = c(class(stat), "refine"))
    snap(stat) <- glue::glue("{.note_refine_cluster}\n计算目标谱系得分所使用的 marker 为 {bind(pos)}，计算非目标谱系得分的 marker 为 {bind(neg)}。")
    stat
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

.note_refine_cluster <- "鉴于全局注释阶段主要关注细胞大类识别，当提取子集 (尤其是子集细胞数相对较少时) 目标细胞子集中仍可能包含转录特征相近的非目标细胞或边缘状态细胞，因此在子集分析前，需要对目标细胞纯度进行再次评估与优化。具体而言，基于目标细胞的经典标志基因集合构建目标谱系得分（positive lineage score），同时基于潜在混杂细胞类型的特异性标志基因集合构建非目标谱系得分（negative lineage score）。该目标谱系得分通过 `Seurat::AddModuleScore` 计算。随后以二者差值定义细胞纯度评分（purity score），用于衡量单细胞转录特征与目标谱系的一致性。进一步按照子集重聚类结果，在簇水平汇总各细胞群的平均纯度评分及细胞数量，识别纯度较低或富集非目标谱系特征的细胞簇，并将其排除。"
