# ==========================================================================
# workflow of aucell
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_aucell <- setClass("job_aucell", 
  contains = c("job"),
  prototype = prototype(
    pg = "aucell",
    info = c("https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html"),
    cite = "[@SCENIC_single_Aibar_2017]",
    method = "",
    tag = "aucell",
    analysis = "AUCell 识别细胞的基因集活性"
    ))

setGeneric("asjob_aucell", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_aucell"))

setMethod("asjob_aucell", signature = c(x = "job_seurat"),
  function(x, sets, assay = "RNA", ...){
    mtx <- object(x)@assays[[assay]]$counts
    if (is.null(mtx) || is.null(rownames(mtx))) {
      stop('is.null(mtx) || is.null(rownames(mtx)).')
    }
    snapAdd_onExit("x", "对单细胞数据集 {x@sig} 以 AUCell 识别细胞的基因集活性 (使用 {snap(sets)}) ")
    x <- job_aucell(mtx, sets, ...)
    return(x)
  })

job_aucell <- function(mtx, sets)
{
  if (!is(mtx, "dgCMatrix")) {
    stop('!is(mtx, "dgCMatrix").')
  }
  if (!is(sets, "GeneSetCollection")) {
    sets <- as_collection(sets)
  }
  gids <- unique(unlist(GSEABase::geneIds(sets)))
  isIns <- gids %in% rownames(mtx)
  message(glue::glue("Has genes: {try_snap(isIns)}"))
  if (all(!isIns)) {
    stop('all(!isIns).')
  }
  x <- .job_aucell(object = mtx)
  x <- methodAdd(x, "以 R 包 `AUCell` ({packageVersion('AUCell')}) {cite_show('SCENIC_single_Aibar_2017')} 识别单细胞数据集的基因集调控活性。")
  x$sets <- sets
  return(x)
}

setMethod("step0", signature = c(x = "job_aucell"),
  function(x){
    step_message("Prepare your data with function `job_aucell`.")
  })

setMethod("step1", signature = c(x = "job_aucell"),
  function(x, workers = NULL){
    step_message("Running...")
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L,
        {
          x <- step1(x, workers = "{workers}")
        }
      )
    } else {
      if (!is.null(workers)) {
        workers <- e(BiocParallel::MulticoreParam(workers))
      }
      x$res <- e(AUCell::AUCell_run(object(x), x$sets, BPPARAM = workers))
      object(x) <- NULL
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_aucell"),
  function(x, ref, type = "AUC", scale = FALSE){
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    if (type == "AUC") {
      res <- t(e(AUCell::getAUC(ref$res)))
      colnames(res) <- paste0("AUC_", colnames(res))
    }
    res <- res[match(rownames(object(x)@meta.data), rownames(res)), , drop = FALSE]
    if (scale) {
      res <- scale(res)
    }
    object(x)@meta.data <- object(x)@meta.data[, !colnames(object(x)@meta.data) %in% colnames(res)]
    object(x)@meta.data <- cbind(object(x)@meta.data, res)
    if (ncol(res) <= 20) {
      x <- focus(x, colnames(res), name = "AUCell", cols = c("skyblue", "blue", "black"))
    }
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_aucell"),
  function(x, wd = glue::glue("~/aucell_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

