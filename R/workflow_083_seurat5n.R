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
  # https://satijalab.org/seurat/articles/parsebio_sketch_integration
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
  # if (!is.null(names)) {
  #   x <- snapAdd(x, "读取 {bind(names)} 样本的数据集。")
  # }
  object(x) <- e(SeuratObject:::merge.Seurat(object(x)[[1]], object(x)[-1]))
  object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(object(x), pattern = "^MT-"))
  p.qc_pre <- plot_qc.seurat(x)
  x$p.qc_pre <- p.qc_pre
  x <- methodAdd(x, "使用 Seurat R 包 ({packageVersion('Seurat')}) 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。")
  return(x)
}

setMethod("step0", signature = c(x = "job_seurat5n"),
  function(x){
    step_message("Prepare your data with function `job_seurat5n`.")
  })

setMethod("step1", signature = c(x = "job_seurat5n"),
  function(x, min.features, max.features, max.count, max.percent.mt = 5)
  {
    step_message("Quality control (QC).")
    if (!is.null(min.features)) {
      object(x) <- e(SeuratObject:::subset.Seurat(
          object(x), subset = nFeature_RNA > min.features &
            nFeature_RNA < max.features & percent.mt < max.percent.mt &
            nCount_RNA < max.count
          ))
      p.qc_aft <- plot_qc.seurat(x)
      p.qc_aft <- set_lab_legend(
        p.qc_aft,
        glue::glue("{x@sig} After Quality control"),
        glue::glue("数据过滤后的 QC 图|||{.seurat_qc_note}") #__REVISE__ set_lab_legend 2026-03-23_21:59:08
      )
      p.qc_pre <- set_lab_legend(
        x$p.qc_pre,
        glue::glue("{x@sig} before Quality control"),
        glue::glue("质量控制 (QC) 图 (数据过滤前) |||{.seurat_qc_note}") #__REVISE__ set_lab_legend 2026-03-23_21:51:34
      )
      x <- plotsAdd(x, p.qc_pre = p.qc_pre, p.qc_aft = p.qc_aft)
      x <- methodAdd(
        x, "前期质量控制{aref(p.qc_pre)}，一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。保留总基因表达量小于 {max.count} 细胞。过滤后{aref(p.qc_aft)}，所有样本共包含{ncol(object(x))}个细胞用于后续分析。" #__REVISE__ methodAdd 2026-03-23_22:06:48
      )
      # x <- methodAdd(x, "一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat5n"),
  function(x, ndims = 20, sct = FALSE, jk = FALSE, workers = NULL){
    step_message("Run standard anlaysis workflow or `SCTransform`.")
    if (is.remote(x)) {
      if (is.null(workers)) {
        stop('is.null(workers).')
      }
      x <- run_job_remote(x, wait = 1,
        {
          x <- step2(x, ndims = "{ndims}", sct = "{sct}", workers = "{workers}")
        }
      )
      return(x)
    }
    if (sct) {
      if (!is.null(workers)) {
        future::plan(future::multicore, workers = workers)
      }
      object(x) <- e(Seurat::SCTransform(
          object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE,
          assay = SeuratObject::DefaultAssay(object(x))
          ))
      message(glue::glue("Shift assays to {SeuratObject::DefaultAssay(object(x))}"))
      x <- methodAdd(
        x, "使用 `Seurat::SCTransform` (默认参数) 对数据集归一化 (<https://satijalab.org/seurat/articles/sctransform_vignette>) 。" #__REVISE__ methodAdd 2026-03-23_22:41:41
      )
    } else {
      object(x) <- e(Seurat::NormalizeData(object(x)))
      object(x) <- e(Seurat::FindVariableFeatures(object(x)))
      object(x) <- e(Seurat::ScaleData(object(x)))
      x <- methodAdd(
        x, "执行标准 Seurat 分析工作流 (`NormalizeData`, `FindVariableFeatures`, `ScaleData`)。"
      )
      p.varfeature <- e(Seurat::VariableFeaturePlot(object(x)))
      p.varfeature <- set_lab_legend(
        wrap(p.varfeature),
        glue::glue("{x@sig} Variable Feature Plot"), #__REVISE__ set_lab_legend 2026-03-23_22:13:35
        glue::glue("高变基因图|||红色代表高变基因，横坐标为基因在所有细胞中的表达水平（log10对数值），纵坐标为基因在所有细胞中的表达水平的标准差，数值越大，表示该基因在细胞中的表达水平越不稳定。生物学差异（如细胞类型、状态等差异）通常会导致某些基因在不同细胞之间表现出较大变异，因此更有可能提供关于生物学现象的信息。") #__REVISE__ set_lab_legend 2026-03-23_22:15:27
      )
      x <- plotsAdd(x, p.varfeature)
    }
    object(x) <- e(Seurat::RunPCA(object(x)))
    x <- methodAdd(x, "随后 PCA 聚类 (`RunPCA`)。")
    if (jk && !sct) {
      object(x) <- e(Seurat::JackStraw(object(x), dims = ndims))
      object(x) <- e(Seurat::ScoreJackStraw(object(x)))
      p.jackPlot <- Seurat::JackStrawPlot(object(x))
      p.jackPlot <- set_lab_legend(
        wrap(p.jackPlot),
        glue::glue("{x@sig} Jack Straw plot"),
        glue::glue("Jackstraw 置换检验 |||通过对原始数据进行多次置换，构建一个零假设分布，然后将实际观测到的主成分得分与该零假设分布进行比较。每个点表示基因在某个主成分上的投影得分与随机背景的比较，大于或等于实际观测主成分得分的比例就是 p 值。p &lt; 0.05 通常认为在该显著性水平下，实际观测到的主成分得分显著高于随机情况下的得分，说明该主成分具有统计学意义，不是由随机因素导致的。通过量化主成分的显著性强度，与均匀分布（虚线）比较，判断哪些主成分更具有统计学意义，富含低p值基因较多的主成分更有统计学意义。")
      )
      x <- plotsAdd(x, p.jackPlot)
      x <- methodAdd(x, "通过 Jackstraw 函数置换检验重新聚类以检验 PC 的选择结果{aref(p.jackPlot)}（P &lt; 0.05）。")
    }
    p.pca_rank <- e(Seurat::ElbowPlot(object(x), ndims))
    # add Seurat::PCAPlot
    p.pca_rank <- set_lab_legend(
      wrap(pretty_elbowplot(p.pca_rank), 4, 4),
      glue::glue("{x@sig} Standard deviations of PCs"),
      glue::glue("主成分 (PC) 的标准化方差 (Standard deviations)|||横坐标为主成分数目，纵坐标代表基于每个主成分对方差解释率的排名（每个主成分的解释方差是其特征值（eigenvalue），表示它解释了总变异的比例），图中每个点表示一个主成分的方差解释比例。") #__REVISE__ set_lab_legend 2026-03-23_22:01:31
    )
    x <- methodAdd(
      x, "使用 ElbowPlot 函数绘制肘图{aref(p.pca_rank)}，帮助确定用于下游分析的主成分以进行后续分析。" #__REVISE__ methodAdd 2026-03-23_22:27:28
    )
    x <- plotsAdd(x, p.pca_rank)
    # x <- snapAdd(x, "数据归一化，PCA 聚类 (Seurat 标准工作流，见方法章节) 后。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat5n"),
  function(x, dims = 1:15, resolution = 1.2,
    use = c("HarmonyIntegration", "CCAIntegration", "RPCAIntegration"), ...)
  {
    step_message("Identify clusters of cells")
    use <- match.arg(use)
    if (!is.null(x$JoinLayers) && x$JoinLayers) {
      message("Job is 'job_seurat5n', but 'JoinLayers' has been performed.")
      object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = use))
      object(x) <- e(
        Seurat::FindClusters(object(x), resolution = resolution, ...)
      )
      object(x) <- e(
        Seurat::RunUMAP(object(x), dims = dims, reduction = use, ...)
      )
    } else {
      if (is.null(x$.before_IntegrateLayers)) {
        object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = "pca"))
        object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution,
            cluster.name = "unintegrated_clusters"))
        object(x) <- e(Seurat::RunUMAP(object(x), dims = dims,
            reduction = "pca", reduction.name = "umap_unintegrated"))
        x$.before_IntegrateLayers <- TRUE
      }
      p.umapUint <-  e(Seurat::DimPlot(object(x), reduction = "umap_unintegrated",
          group.by = c("orig.ident", "seurat_clusters"), cols = color_set(TRUE)))
      p.umapUint <- set_lab_legend(
        wrap(p.umapUint, 10, 5),
        glue::glue("{x@sig} UMAP Unintegrated"),
        glue::glue("去除批次效应之前的 UMAP 聚类图|||不同颜色代表不同cluster。横纵坐标是 UMAP 降维的两个维度。UMAP能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。") #__REVISE__ set_lab_legend 2026-03-23_22:44:03
      )
      x <- plotsAdd(x, p.umapUint)
      ## integrated
      methods <- list(CCAIntegration = Seurat::CCAIntegration,
        HarmonyIntegration = Seurat::HarmonyIntegration,
        RPCAIntegration = Seurat::RPCAIntegration
      )
      use <- match.arg(use, names(methods))
      object <- object(x)
      res <- try(e(Seurat::IntegrateLayers(object = object,
            method = methods[[ use ]], orig.reduction = "pca",
            new.reduction = use, verbose = FALSE,
            normalization.method = if (object@active.assay == "SCT") "SCT" else "LogNormalize")))
      if (!inherits(res, "try-error")) {
        object(x) <- res
      } else {
        warning("Got error while perform `Seurat::IntegrateLayers`, return the job.")
        return(x)
      }
      object(x)[["RNA"]] <- e(SeuratObject::JoinLayers(object(x)[["RNA"]]))
      ## SeuratObject::DefaultDimReduc, search in case of UMAP
      object(x)@reductions$umap_unintegrated <- NULL
      ## method of job_seurat
      x <- callNextMethod(
        x, dims, resolution, reduction = use, ...
      )
      x <- methodAdd(
        x, "以 `Seurat::IntegrateLayers` 集成数据，去除批次效应 (使用 {use} 方法)。", add = FALSE
      )
    }
    p.umapInt <-  e(Seurat::DimPlot(object(x),
        group.by = c("orig.ident", "seurat_clusters"), cols = color_set(TRUE)))
    p.umapInt <- set_lab_legend(
      wrap(p.umapInt, 10, 5),
      glue::glue("{x@sig} UMAP Integrated"), #__REVISE__ set_lab_legend 2026-03-23_22:45:48
      glue::glue("去除批次效应之后的 UMAP 聚类图|||不同颜色代表不同cluster。横纵坐标是 UMAP 降维的两个维度。UMAP能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
    )
    p.umapLabel <-  e(Seurat::DimPlot(object(x),
        group.by = c("seurat_clusters"), 
        cols = color_set(TRUE), label = TRUE))
    x <- plotsAdd(x, p.umapInt, p.umapLabel)
    x <- methodAdd(x, "在 1-{max(dims)} PC 维度下，以 `Seurat::FindNeighbors` 构建 Nearest-neighbor Graph。随后在 {resolution} 分辨率下，以 `Seurat::FindClusters` 函数识别细胞群并以 `Seurat::RunUMAP` 进行 UMAP 聚类。")
    x <- methodAdd(x, "在去除批次效应前，UMAP 图{aref(p.umapUint)}中各样本保持离散。去除批次效应后{aref(p.umapInt)}，各样本相互均匀混合，即批次效应已被良好地处理。")
    x$JoinLayers <- TRUE
    return(x)
  })

setMethod("asjob_limma", signature = c(x = "job_seurat"),
  function(x, features, cells, cell_groups, group.by = x$group.by,
    slot = "data", fun_norm = function(x) log2(x + 1), gname = TRUE)
  {
    metadata <- dplyr::mutate(
      as_tibble(object(x)@meta.data), 
      sample = rownames, group = !!rlang::sym(group.by), .before = 1
    )
    if (missing(cells) && !missing(cell_groups)) {
      cells <- metadata[[ group.by ]] %in% cell_groups
    }
    metadata <- metadata[ cells, ]
    assay <- object(x)@assays[[ object(x)@active.assay ]]
    data <- SeuratObject::LayerData(object = assay, layer = slot)
    genes <- data.frame(gene = rownames(data))
    if (gname) {
      genes <- dplyr::mutate(genes, gene = gname(gene))
    }
    if (is.character(features)) {
      features <- genes$gene %in% features
    }
    genes <- genes[ features, , drop = FALSE]
    data <- data[ features, cells ]
    data <- data.frame(data, check.names = FALSE)
    data <- fun_norm(data)
    object <- new_from_package(
      "EList", "limma", list(E = data, targets = metadata, genes = genes)
    )
    validObject(object)
    x <- .job_limma()
    x$normed_data <- object
    x$from_seurat <- TRUE
    return(x)
  })

.seurat_qc_note <- "**上方小提琴图**：每个样本对应一个‘小提琴’，小提琴的宽度代表相应数据的密度，宽度越大表示在该区域内的数据点越密集，更多数据点集中于此区域；宽度越小则表示密度越小，即数据相对较少；过滤标准：nCount 和 nFeature 过高可能是双细胞，过低可能是细胞碎片；percent.mt（线粒体基因表达比例，是细胞内线粒体基因表达量占所有基因表达量的比例）表明细胞状态，值过高可能是细胞正在经历压力或死亡。**下方点图**：每一个点代表一个细胞，不同颜色代表不同样本；左图横坐标为总基因表达数，纵坐标为线粒体基因比例；右图横坐标为 nCount（总基因表达数），纵坐标为 Feature（总基因数）；正常情况下，nCount 越多那么 nFeature 就越高，呈现出正相关关系，因此检测到的基因表达数应与检测到的基因数目在细胞间高度相关，而线粒体基因比例则不相关（若呈正相关，横坐标越大纵坐标也越大；若呈负相关，横坐标越大纵坐标应越小）。"

