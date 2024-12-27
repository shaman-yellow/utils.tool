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
    info = c("https://satijalab.org/seurat/articles"),
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019]",
    method = "R package `Seurat` used for multiple dataset integration",
    tag = "scrna:anno",
    analysis = "Seurat 集成单细胞数据分析"
    ))

job_seuratn <- function(dirs, names = NULL, mode = c("sc", "st"), st.filename = "filtered_feature_bc_matrix.h5")
{
  mode <- match.arg(mode)
  n <- 0L
  object <- pbapply::pblapply(dirs,
    function(dir) {
      n <<- n + 1L
      project <- names[n]
      if (mode == "sc") {
        suppressMessages(job_seurat(dir, project = project))@object
      } else {
        suppressMessages(job_seuratSp(dir, filename = st.filename))@object
      }
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
    fakeRNAassayFromSpatial <- FALSE
    object(x) <- e(lapply(object(x),
        function(x) {
          if (any(names(x@assays) == "Spatial")) {
            # https://github.com/satijalab/seurat/issues/8216
            message("Due to issue, copy 'Spatial' to 'RNA'.")
            x[["RNA"]] <- x[["Spatial"]]
            fakeRNAassayFromSpatial <<- TRUE
          }
          if (!any(names(x@assays) == "SCT")) {
            Seurat::SCTransform(x,
              method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE,
              assay = x@active.assay
            )
          } else x
        }))
    x@params$anchor.features <- e(Seurat::SelectIntegrationFeatures(object(x), nfeatures = nfeatures))
    object(x) <- e(Seurat::PrepSCTIntegration(object(x), anchor.features = x@params$anchor.features))
    if (fakeRNAassayFromSpatial) {
      object(x) <- lapply(object(x),
        function(x) {
          x[[ "RNA" ]] <- NULL
          x
        })
    }
    meth(x)$step1 <- glue::glue("使用 Seurat R 包 ({packageVersion('Seurat')}) 进行数据质量控制 (QC) 和下游分析。一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。使用 Seurat::SCTransform 函数对数据标准化。随后，以 Seurat::SelectIntegrationFeatures (返还{nfeatures}个基因) 和 Seurat::PrepSCTIntegration 选择用于整合多重数据集的基因。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_seuratn"),
  function(x, workers = NULL){
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
    x <- .job_seurat(object = object(x), step = 2L, meth = meth(x), sig = x@sig)
    x@plots[[ 2 ]] <- plot_pca.seurat(object(x))
    meth(x)$step2 <- glue::glue("以 Seurat::FindIntegrationAnchors，和 Seurat::IntegrateData 整合多个数据集，并 PCA 降维。")
    x
  })
