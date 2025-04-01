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
  function(x, ..., groups = NULL, group.by = x$group.by){
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
    x$group.by <- group.by
    return(x)
  })

setMethod("step0", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Prepare your data with function `job_seurat_sub`.")
  })

setMethod("step1", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Do nothing.")
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat_sub"),
  function(x){
    step_message("Reset variable features.")
    object(x) <- e(Seurat::FindVariableFeatures(object(x)))
    object(x) <- e(Seurat::ScaleData(object(x)))
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

