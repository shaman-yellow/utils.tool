# ==========================================================================
# workflow of seurat5n
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat5n <- setClass("job_seurat5n", 
  contains = c("job_seurat"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "seurat5n",
    info = c("https://satijalab.org/seurat/articles/integration_introduction"),
    cite = "[@DictionaryLearHaoY2024]",
    method = "",
    tag = "seurat5n",
    analysis = "Seurat 集成单细胞数据分析"
    ))

job_seurat5n <- function(dirs, names = NULL, mode = c("sc", "st"), st.filename = "filtered_feature_bc_matrix.h5")
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
  x <- .job_seurat5n(object = object)
  object(x) <- e(SeuratObject:::merge.Seurat(object(x)[[1]], object(x)[-1]))
  object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(object(x), pattern = "^MT-"))
  p.qc_pre <- plot_qc.seurat(x)
  p.qc_pre <- .set_lab(p.qc_pre, sig(x), "Pre-Quality control")
  x@params$p.qc_pre <- p.qc_pre
  meth(x)$step0 <- glue::glue("使用 Seurat R 包 ({packageVersion('Seurat')}) 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。")
  return(x)
}

setMethod("step0", signature = c(x = "job_seurat5n"),
  function(x){
    step_message("Prepare your data with function `job_seurat5n`.")
  })

setMethod("step1", signature = c(x = "job_seurat5n"),
  function(x, min.features, max.features, max.percent.mt = 5)
  {
    step_message("Quality control (QC).")
    if (!is.null(min.features)) {
      object(x) <- e(SeuratObject:::subset.Seurat(
          object(x), subset = nFeature_RNA > min.features &
            nFeature_RNA < max.features & percent.mt < max.percent.mt
          ))
      p.qc_aft <- plot_qc.seurat(x)
      p.qc_aft <- .set_lab(p.qc_aft, sig(x), "After Quality control")
      x@params$p.qc_aft <- p.qc_aft
      meth(x)$step1 <- glue::glue("一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat5n"),
  function(x, ndims = 20){
    step_message("Run standard anlaysis workflow")
    object(x) <- e(Seurat::NormalizeData(object(x)))
    object(x) <- e(Seurat::FindVariableFeatures(object(x)))
    object(x) <- e(Seurat::ScaleData(object(x)))
    object(x) <- e(Seurat::RunPCA(object(x)))
    p.pca_rank <- e(Seurat::ElbowPlot(object(x), ndims))
    p.pca_rank <- wrap(pretty_elbowplot(p.pca_rank), 4, 4)
    p.pca_rank <- .set_lab(p.pca_rank, sig(x), "Standard deviations of PCs")
    x@plots[[ 2 ]] <- namel(p.pca_rank)
    meth(x)$step <- glue::glue("执行标准 Seurat 分析工作流 (`NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`)。以 `ElbowPlot` 判断后续分析的 PC 维度。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat5n"),
  function(x, dims = 1:15, resolution = 2, use = "HarmonyIntegration")
  {
    step_message("Identify clusters of cells")
    object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = "pca"))
    object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution,
        cluster.name = "unintegrated_clusters"))
    object(x) <- e(Seurat::RunUMAP(object(x), dims = dims,
        reduction = "pca", reduction.name = "umap_unintegrated"))
    p.umapUint <-  e(Seurat::DimPlot(object(x), reduction = "umap_unintegrated",
        group.by = c("orig.ident", "seurat_clusters"), cols = color_set(TRUE)))
    p.umapUint <- .set_lab(wrap(p.umapUint, 10, 5), sig(x), "UMAP Unintegrated")
    ## integrated
    methods <- list(CCAIntegration = Seurat::CCAIntegration,
      HarmonyIntegration = Seurat::HarmonyIntegration)
    use <- match.arg(use, names(methods))
    object(x) <- e(Seurat::IntegrateLayers(object = object(x),
        method = methods[[ use ]], orig.reduction = "pca",
        new.reduction = use, verbose = FALSE))
    object(x)[["RNA"]] <- e(SeuratObject::JoinLayers(object(x)[["RNA"]]))
    ## SeuratObject::DefaultDimReduc, search in case of UMAP
    object(x)@reductions$umap_unintegrated <- NULL
    ## method of job_seurat
    x <- callNextMethod(x, dims, resolution, reduction = use)
    p.umapInt <-  e(Seurat::DimPlot(object(x),
        group.by = c("orig.ident", "seurat_clusters"), cols = color_set(TRUE)))
    p.umapInt <- .set_lab(wrap(p.umapInt, 10, 5), sig(x), "UMAP Integrated")
    plots <- namel(p.umapUint, p.umapInt)
    x@plots[[ 3 ]] <- c(x@plots[[ 3 ]], plots)
    meth(x)$step3 <- glue::glue("以 `Seurat::IntegrateLayers` 集成数据，去除批次效应 (使用 {use} 方法)。在 1-{max(dims)} PC 维度下，以 `Seurat::FindNeighbors` 构建 Nearest-neighbor Graph。随后在 {resolution} 分辨率下，以 `Seurat::FindClusters` 函数识别细胞群并以 `Seurat::RunUMAP` 进行 UMAP 聚类。")
    return(x)
  })
