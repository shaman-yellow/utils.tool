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
    info = c("Tutorial: https://github.com/satijalab/seurat/wiki")
    ))

job_seuratn <- function(dirs)
{
  object <- pbapply::pblapply(dirs,
    function(dir) {
      suppressMessages(job_seurat(dir))@object
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
  function(x){
    step_message(".
      Performing integration on datasets normalized with SCTransform.
      "
    )
    object(x) <- e(lapply(object(x),
        function(obj) {
          obj[[ "percent.mt" ]] <- Seurat::PercentageFeatureSet(
            obj, pattern = "^MT-"
          )
          return(obj)
        }))
    object(x) <- e(lapply(object(x), Seurat::SCTransform,
        method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T,
        ))
    features <- e(Seurat::SelectIntegrationFeatures(object(x), nfeatures = 3000))
    object(x) <- e(Seurat::PrepSCTIntegration(object(x), anchor.features = features))
    object(x) <- e(Seurat::FindIntegrationAnchors(object(x), normalization.method = "SCT",
        anchor.features = features))
    object(x) <- e(Seurat::IntegrateData(object(x), normalization.method = "SCT"))
    object(x) <- e(Seurat::RunPCA(object(x), verbose = FALSE))
    return(x)
  })

setMethod("step2", signature = c(x = "job_seuratn"),
  function(x){
    step_message("Transform the job as 'job_seurat' in step 2.
      And conduct some visualization.
      "
    )
    x <- .job_seurat(object = object(x), step = 2L)
    x@plots[[ 2 ]] <- plot_pca.seurat(object(x))
    x
  })
