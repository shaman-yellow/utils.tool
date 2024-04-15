# ==========================================================================
# workflow of seuratn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seuratn <- setClass("job_seuratn", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/seurat/wiki"),
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019]",
    method = "R package `Seurat` used for multiple dataset integration"
    ))

job_seuratn <- function(dirs, names = NULL)
{
  n <- 0L
  object <- pbapply::pblapply(dirs,
    function(dir) {
      n <<- n + 1L
      project <- names[n]
      suppressMessages(job_seurat(dir, project = project))@object
    })
  .job_seuratn(object = object)
}

setMethod("step0", signature = c(x = "job_seuratn"),
  function(x){
    step_message("Prepare your data with function `job_seuratn`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_seuratn"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000){
    step_message("
      QC and Preparing integration on datasets normalized with SCTransform.
      "
    )
    if (missing(min.features) | missing(max.features))
      stop("missing(min.features) | missing(max.features)")
    object(x) <- e(lapply(object(x),
        function(obj) {
          obj[[ "percent.mt" ]] <- Seurat::PercentageFeatureSet(
            obj, pattern = "^MT-"
          )
          if (!is.null(min.features)) {
            obj <- SeuratObject:::subset.Seurat(
                obj, subset = nFeature_RNA > min.features &
                  nFeature_RNA < max.features & percent.mt < max.percent.mt
                )
          }
          return(obj)
        }))
    object(x) <- e(lapply(object(x), Seurat::SCTransform,
        method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T,
        ))
    x@params$anchor.features <- e(Seurat::SelectIntegrationFeatures(object(x), nfeatures = 3000))
    object(x) <- e(Seurat::PrepSCTIntegration(object(x), anchor.features = x@params$anchor.features))
    return(x)
  })

setMethod("step2", signature = c(x = "job_seuratn"),
  function(x, workers = 3){
    step_message("Perform integration (time-consumed).
      Prarameter red{{workers}} would call `future::plan`.
      If it meet error (Seurat::IntegrateData), the object would be returned.
      Transform the job as 'job_seurat'.  And conduct some visualization.
      "
    )
    fun <- function(x) {
      object(x) <- e(Seurat::FindIntegrationAnchors(object(x), normalization.method = "SCT",
          anchor.features = x@params$anchor.features))
      return(x)
    }
    if (!is.null(workers)) {
      x <- parallel(x, fun, workers)
    } else {
      x <- fun(x)
    }
    res <- try(e(Seurat::IntegrateData(object(x), normalization.method = "SCT")))
    if (inherits(res, "try-error")) {
      return(x)
    } else {
      object(x) <- res
    }
    object(x) <- e(Seurat::RunPCA(object(x), verbose = FALSE))
    x <- .job_seurat(object = object(x), step = 2L)
    x@plots[[ 2 ]] <- plot_pca.seurat(object(x))
    x
  })
