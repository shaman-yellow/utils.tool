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
  function(x, ..., groups = NULL, group.by = x$group.by, sample = 1L)
  {
    sub_from <- class(x)
    if (!missing(sample)) {
      seed <- x$seed
      x <- getsub(x, sample = sample, ...)
      x <- .job_seurat_sub(object = object(x))
      x$seed <- seed
      x <- snapAdd(x, "根据样本和细胞类型分组，将细胞随机抽样 (各组比例为：{sample}) (细胞数量较多，通过随机抽样的方式减少计算负担) (随机种子：{x$seed})。")
    } else {
      object <- object(x)
      x <- .job_seurat_sub()
      if (!is.null(groups)) {
        cells <- which(object@meta.data[[ group.by ]] %in% groups)
        x <- snapAdd(x, "提取 {bind(groups)}，分析亚群。")
      } else {
        data <- trace_filter(object@meta.data, ...)
        cells <- rownames(data)
        cells <- which(rownames(object@meta.data) %in% cells)
        x <- snapAdd(x, "{snap(data)}分析其亚群。")
      }
      object(x) <- e(SeuratObject:::subset.Seurat(object, cells = cells))
    }
    x$group.by <- group.by
    x$sub_from <- sub_from
    return(x)
  })

setMethod("step0", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Prepare your data with function `job_seurat_sub`.")
  })

setMethod("step1", signature = c(x = "job_seurat_sub"),
  function(x, assay = "RNA", by = "orig.ident"){
    if (x$sub_from == "job_seurat5n") {
      step_message("Split datasets, and render as `job_seurat5n`.")
      object(x)[[assay]] <- split(
        object(x)[[assay]], 
        do.call(`$`, list(object(x), by))
      )
      x <- .job_seurat5n(x)
    } else {
      step_message("Do nothing.")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat_sub"),
  function(x, reset = FALSE, sct = FALSE) {
    step_message("Reset variable features.")
    if (reset || sct) {
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
      p.pca_rank <- e(Seurat::ElbowPlot(object(x), 30))
      p.pca_rank <- wrap(pretty_elbowplot(p.pca_rank), 4, 4)
      p.pca_rank <- .set_lab(p.pca_rank, sig(x), "Standard deviations of PCs")
      p.pca_rank <- setLegend(p.pca_rank, "为主成分 (PC) 的 Standard deviations。")
      x@plots[[ 2 ]] <- namel(p.pca_rank)
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat_sub"),
  function(x, dims = 1:15, resolution = 1.2, force = TRUE, ...)
  {
    x <- callNextMethod(x, dims, resolution, force = force, ...)
    return(x)
  })

as_markers <- function(cell_markers) {
  as_df.lst(cell_markers, "cell", "markers")
}

