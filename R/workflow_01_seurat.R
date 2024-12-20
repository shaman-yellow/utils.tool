# ==========================================================================
# workflow of seurat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat <- setClass("job_seurat", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html"),
    cite = "[@IntegratedAnalHaoY2021; @ComprehensiveIStuart2019]",
    method = "The R package `Seurat` used for scRNA-seq processing",
    tag = "scrna:anno",
    analysis = "Seurat 单细胞数据分析"
    ))

job_seurat <- function(dir = NULL, project = basename(sub("/$", "", dir)),
  min.cells = 3, min.features = 200, file_h5 = NULL, ...)
{
  if (!is.null(file_h5)) {
    data <- e(Seurat::Read10X_h5(file_h5))
  } else {
    data <- e(Seurat::Read10X(dir))
  }
  object <- e(Seurat::CreateSeuratObject(counts = data, project = project,
      min.cells = min.cells, min.features = min.features, ...))
  x <- .job_seurat(object = object)
  meth(x)$step0 <- glue::glue("")
  x
}

setMethod("step0", signature = c(x = "job_seurat"),
  function(x){
    step_message("Prepare your data with function `job_seurat`. ",
      crayon::red("Required: Seurat v5. "),
      "The `dir` to read would passed to `Seurat::Read10X` to ", 
      "load single-cell data. The directory must contains three ",
      "files: ", crayon::red("barcodes.tsv"), "; ",
      crayon::red("genes.tsv"), "; ", crayon::red("matrix.mtx"), "."
    )
  })

setMethod("step1", signature = c(x = "job_seurat"),
  function(x){
    step_message("QC and selecting cells for further analysis.",
      "This do:",
      "Generate column in `object(x)@meta.data`; ",
      "Plots in `x@plots[[ 1 ]]`.",
      "Check the plots then decided the ",
      crayon::red("min nFeature_RNA"), " and ",
      crayon::red("max nFeature_RNA"), ", ",
      "as well as ", crayon::red("percent.mt"), "."
    )
    object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(
      object(x), pattern = "^MT-"
    ))
    p.qc <- plot_qc.seurat(x)
    p.qc <- .set_lab(p.qc, sig(x), "Quality", "Control")
    x@plots[[ 1 ]] <- list(p.qc = p.qc)
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat"),
  function(x, min.features, max.features, max.percent.mt = 5, nfeatures = 2000,
    use = "nFeature_RNA")
  {
    step_message("This contains several execution: Subset the data,
      Normalization, Feature selection, Scale the data, Linear dimensional
      reduction, Select dimensionality.
      red{{`min.features`}} and red{{`max.features`}} were needed for subset.
      Then `object(x)` were performed with:
      `Seurat::SCTransform`;
      `Seurat::RunPCA`.
      All plots were in `x@plots[[ 2 ]]`
      NOTE: inspect the plots and red{{Determined dims}} for downstream analysis. "
    )
    if (missing(min.features) | missing(max.features))
      stop("missing(min.features) | missing(max.features)")
    object(x) <- e(SeuratObject:::subset.Seurat(
        object(x), subset = !!rlang::sym(use) > min.features &
          !!rlang::sym(use) < max.features & percent.mt < max.percent.mt
        ))
    object(x) <- e(Seurat::SCTransform(
        object(x), method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T,
        assay = SeuratObject::DefaultAssay(object(x))
        ))
    object(x) <- e(Seurat::RunPCA(
        object(x), features = SeuratObject::VariableFeatures(object(x))
        ))
    x@plots[[ 2 ]] <- plot_pca.seurat(object(x))
    x@plots[[ 2 ]]$p.pca_rank <- .set_lab(
      x@plots[[ 2 ]]$p.pca_rank, sig(x), "Ranking of principle components")
    meth(x)$step1 <- glue::glue("使用 Seurat R 包 ({packageVersion('Seurat')}) 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。建议线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。使用 Seurat::SCTransform 函数对数据标准化。再PCA降维。")
    message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat"),
  function(x, dims = 1:15, resolution = 1.2, reduction = "pca", force = F){
    step_message("This contains several execution: Cluster the cells, UMAP
      reduction, Cluster biomarker (Differential analysis). Object were
      performed with:
      `Seurat::FindNeighbors`;
      `Seurat::FindClusters` grey{{(clusters can be accessed by `SeuratObject::Idents`)}};
      yellow{{The parameter `resolution` between 0.4-1.2 typically returns good results for
      single-cell datasets of around 3000 cells (ncol). Optimal resolution
      often increases for larger datasets.}}
      grey{{Consider use Seurat::FindMarkers Seurat::FindAllMarkers(..., test.use = 'roc') for
      customize differential analysis.}}
      grey{{Consider use `Seurat::VlnPlot` or `Seurat::FeaturePlot`}} for
      customize mapping of features."
    )
    if (missing(dims) | missing(resolution))
      stop("missing(dims) | missing(resolution)")
    if (!force) {
      if (ncol(object(x)) > 3000 & resolution < 1.2) {
        stop("ncol(object(x)) > 3000 but resolution < 1.2. ",
          "Please consider higher `resolution`. ")
      }
    }
    # if (isNamespaceLoaded("future")) {
      # e(pkgload::unload("future"))
    # }
    object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = reduction))
    object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution))
    object(x) <- e(Seurat::RunUMAP(object(x), dims = dims, reduction = reduction))
    p.umap <- e(Seurat::DimPlot(object(x), cols = color_set(T)))
    p.umap <- wrap(p.umap, 6, 5)
    p.umap <- .set_lab(p.umap, sig(x), "UMAP", "Clustering")
    x@plots[[ 3 ]] <- namel(p.umap)
    meth(x)$step3 <- glue::glue("在 1-{max(dims)} PC 维度下，以 Seurat::FindNeighbors 构建 Nearest-neighbor Graph。随后在 {resolution} 分辨率下，以 Seurat::FindClusters 函数识别细胞群并以 Seurat::RunUMAP 进行 UMAP 聚类。")
    return(x)
  })

setMethod("step4", signature = c(x = "job_seurat"),
  function(x, use = "SingleR", use.level = c("label.main", "label.fine"),
    ref = celldex::HumanPrimaryCellAtlasData())
  {
    if (use == "scAnno") {
      ## PMID: 37183449
      ## https://github.com/liuhong-jia/scAnno
      stop("Deprecated. Too many bugs in that package.")
    } else if (use == "SingleR") {
      step_message("Use `SingleR` and `celldex` to annotate cell types.
        By default, red{{`celldex::HumanPrimaryCellAtlasData`}} was used
        as red{{`ref`}} dataset. This annotation would generate red{{'SingleR_cell'}}
        column in `object(x)@meta.data`. Plots were generated in `x@plots[[ 4 ]]`;
        tables in `x@tables[[ 4 ]]`.
        "
      )
      ref <- e(celldex::HumanPrimaryCellAtlasData())
      clusters <- object(x)@meta.data$seurat_clusters
      use.level <- match.arg(use.level)
      anno_SingleR <- e(SingleR::SingleR(object(x)@assays$SCT@scale.data,
          ref = ref, labels = ref[[ use.level ]], clusters = clusters
          ))
      score <- as.matrix(anno_SingleR$scores)
      rownames(score) <- rownames(anno_SingleR)
      p.score_SingleR <- callheatmap(
        new_heatdata(as_data_long(score, 
            row_var = "Cluster", col_var = "Cell_type"))
      )
      p.score_SingleR <- wrap(p.score_SingleR, 30, 8)
      anno_SingleR <- tibble::tibble(
        seurat_clusters = rownames(anno_SingleR),
        SingleR_cell = anno_SingleR$labels
      )
      object(x)@meta.data$SingleR_cell <-
        anno_SingleR$SingleR_cell[match(clusters,
          anno_SingleR$seurat_clusters)]
      p.map_SingleR <- e(Seurat::DimPlot(
          object(x), reduction = "umap", label = F, pt.size = .7,
          group.by = "SingleR_cell", cols = color_set()
          ))
      p.map_SingleR <- p.map_SingleR
      p.map_SingleR <- wrap(p.map_SingleR, 10, 7)
      x@tables[[ 4 ]] <- namel(anno_SingleR)
      x@plots[[ 4 ]] <- namel(p.score_SingleR, p.map_SingleR)
      x@params$group.by <- "SingleR_cell"
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_seurat"),
  function(x, workers = NULL, min.pct = .25, logfc.threshold = .25)
  {
    step_message("Find all Marders for Cell Cluster.")
    fun <- function(x) {
      markers <- as_tibble(
        e(Seurat::FindAllMarkers(object(x), min.pct = min.pct,
            logfc.threshold = logfc.threshold, only.pos = T))
      )
    }
    if (!is.null(workers)) {
      markers <- parallel(x, fun, workers)
    } else {
      markers <- fun(x)
    }
    all_markers_no_filter <- markers
    markers <- dplyr::filter(markers, p_val_adj < .05)
    tops <- dplyr::slice_max(dplyr::group_by(markers, cluster), avg_log2FC, n = 10)
    if (F) {
      p.toph <- e(Seurat::DoHeatmap(object(x), features = tops$gene, raster = T))
      p.toph <- wrap(p.toph, 14, 12)
      x@plots[[ 5 ]] <- namel(p.toph)
    }
    x@tables[[ 5 ]] <- list(all_markers = markers, all_markers_no_filter = all_markers_no_filter)
    meth(x)$step5 <- glue::glue("以 `Seurat::FindAllMarkers` (LogFC 阈值 {logfc.threshold}; 最小检出率 {min.pct}) 为所有细胞群寻找 Markers。")
    return(x)
  })

setMethod("step6", signature = c(x = "job_seurat"),
  function(x, tissue, ref.markers = NULL, filter.p = 0.01, filter.fc = 1.5, filter.pct = .7,
    org = c("Human", "Mouse"),
    cmd = pg("scsa"), db = pg("scsa_db"), res.col = "scsa_cell",
    method = c("gpt", "scsa"), n = 30)
  {
    method <- match.arg(method)
    if (method == "gpt") {
      step_message("Use ChatGPT (ChatGPT-4o) to annotate cell types.")
      all_markers <- x@tables$step5$all_markers
      if (!is.null(filter.pct)) {
        all_markers <- dplyr::filter(all_markers, pct.1 >= filter.pct)
      }
      query <- prepare_GPTmessage_for_celltypes(tissue, all_markers, n = n)
      if (usethis::ui_yeah("Get results from clipboard?")) {
        feedback <- get_clipboard()
        if (length(feedback) != query$ncluster) {
          stop(glue::glue("The res number ({nrow(res)}) not match cluster number ({query$ncluster})."))
        }
        res <- parse_GPTfeedback(feedback)
        seqs <- match(object(x)@meta.data[[ "seurat_clusters"]], 0:(query$ncluster - 1))
        object(x)@meta.data[[ "ChatGPT_cell"]] <- res$cells[ seqs ]
        ## heatmap
        p.markers <- e(Seurat::DoHeatmap(object(x), features = res$markers,
            group.by = "ChatGPT_cell", raster = T, group.colors = color_set(), label = F))
        p.markers <- .set_lab(p.markers, sig(x), "Markers in cell types")
        attr(p.markers, "lich") <- new_lich(namel(ChatGPT_Query = query$query, feedback), sep = "\n")
        ## dim plot
        p.map_gpt <- e(Seurat::DimPlot(
            object(x), reduction = "umap", label = F, pt.size = .7,
            group.by = "ChatGPT_cell", cols = color_set()
            ))
        p.map_gpt <- wrap(as_grob(p.map_gpt), 7, 4)
        p.map_gpt <- .set_lab(p.map_gpt, sig(x), "ChatGPT", "Cell type annotation")
        x@plots[[ 6 ]] <- namel(p.map_gpt, p.markers)
        meth(x)$step6 <- glue::glue("以 ChatGPT-4o 注释细胞类型 {cite_show('Assessing_GPT_4_Hou_W_2024')}。将每一个细胞群的 Top {n} 基因提供给 ChatGPT，使其注释细胞类型。询问信息为：\n\n{text_roundrect(query$message)}")
      } else {
        stop("Terminated.")
      }
    } else {
      step_message("Use SCSA to annotate cell types (<https://github.com/bioinfo-ibms-pumc/SCSA>).")
      lst <- do.call(scsa_annotation, as.list(environment()))
      x <- lst$x
      x@tables[[ 6 ]] <- list(scsa_res_all = lst$scsa_res_all)
      x@params$group.by <- lst$res.col
      x@plots[[ 6 ]] <- list(p.map_scsa = lst$p.map_scsa)
    }
    return(x)
  })

setMethod("step7", signature = c(x = "job_seurat"),
  function(x, classifier, db = org.Hs.eg.db::org.Hs.eg.db){
    step_message("
      Use `garnett::classify_cells` to anntate cells.
      Prarameter red{{`classifier`}} specify the pre-difined classifier.
      "
    )
    object <- e(SeuratWrappers::as.cell_data_set(object(x),
        group.by = x@params$group.by, ...))
    object <- e(garnett::classify_cells(object, classifier, db = db, 
      cds_gene_id_type = "SYMBOL"))
    object(x) <- e(SeuratWrappers:::as.Seurat.cell_data_set(object))
    x@params$group.by <- "cell_type"
    return(x)
  })

setMethod("diff", signature = c(x = "job_seurat"),
  function(x, group.by, contrasts, name = "contrasts", force = T)
  {
    if (is.data.frame(contrasts)) {
      contrasts <- apply(contrasts, 1, c, simplify = F)
    }
    if (is.null(x@params[[ name ]]) || force) {
      res <- e(lapply(contrasts,
          function(con) {
            data <- Seurat::FindMarkers(object(x),
              ident.1 = con[1], ident.2 = con[2],
              group.by = group.by
            )
            data$gene <- rownames(data)
            data
          }))
      names(res) <- vapply(contrasts, function(x) paste0(x[1], "_vs_", x[2]), character(1))
      res <- dplyr::as_tibble(data.table::rbindlist(res, idcol = T))
      res <- dplyr::rename(res, contrast = .id)
      res <- dplyr::filter(res, p_val_adj < .05)
      res <- .set_lab(res, sig(x), "DEGs of the contrasts")
      x@params[[ name ]] <- res
    } else {
      res <- x@params[[ name ]]
    }
    ## contrast intersection
    tops <- split(res, ~ contrast)
    use.gene <- "gene"
    fun_filter <- function(x) rm.no(x)
    tops <- lapply(tops,
      function(data){
        up <- dplyr::filter(data, avg_log2FC > 0)[[ use.gene ]]
        down <- dplyr::filter(data, avg_log2FC < 0)[[ use.gene ]]
        lst <- list(up = up, down = down)
        lapply(lst, fun_filter)
      })
    tops <- unlist(tops, recursive = F)
    x[[ paste0(name, "_intersection") ]] <- tops
    p.sets_intersection <- new_upset(lst = tops, trunc = NULL)
    p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "contrasts-DEGs-intersection")
    x[[ paste0("p.", name, "_intersection") ]] <- p.sets_intersection
    return(x)
  })

diff_group <- function(x, group.by, patterns) {
  group <- combn(as.character(ids(x, group.by)), 2, simplify = F)
  try_contrast(group, patterns)
}

try_contrast <- function(combn, patterns) {
  lapply(combn,
    function(x) {
      res <- lapply(patterns,
        function(pattern) {
          matched <- grepl(pattern, x)
          if (any(matched))
            c(x[matched], x[!matched])
        })
      res <- lst_clear0(res)
      res[[1]]
    })
}

setMethod("map", signature = c(x = "job_seurat", ref = "marker_list"),
  function(x, ref, p.adjust = .05, log2fc = 1, use_all = T, prop = .6, print = T){
    if (x@step < 5L)
      stop("x@step < 5L")
    all_markers <- tables(x, 5)$all_markers
    all_markers <- filter(all_markers, p_val_adj < p.adjust, avg_log2FC > log2fc)
    all_markers <- select(all_markers, cluster, gene)
    all_markers <- split(all_markers$gene, all_markers$cluster)
    refs <- lapply(ref@level, function(n) ref@marker[[ n ]])
    stats <- lapply(all_markers,
      function(markers) {
        alls <- lapply(refs,
          function(ref) {
            all <- vapply(ref, FUN.VALUE = numeric(1),
              function(ref.) {
                res <- prop.table(table(ref. %in% markers))
                res <- unname(res[ names(res) == "TRUE" ])
                if (!length(res)) 0
                else res
              })
            tibble::tibble(cell_type = names(all), prop = unname(all))
          })
        names(alls) <- ref@level
        alls <- data.table::rbindlist(alls, fill = T, idcol = T)
        rename(alls, unique_level = .id)
      })
    map_prop <- data.table::rbindlist(stats, fill = T, idcol = T)
    map_prop <- rename(map_prop, cluster = .id)
    x@tables$map_prop <- map_prop
    map_prop <- mutate(map_prop, unique_level = as.integer(unique_level))
    map_prop <- arrange(map_prop, dplyr::desc(prop))
    if (use_all) {
      mapu <- distinct(map_prop, cell_type, .keep_all = T)
      mapu <- rbind(mapu, map_prop)
      map_res <- distinct(mapu, cluster, .keep_all = T)
    } else {
      map_res <- filter(map_prop, prop > !!prop)
      map_res <- distinct(map_res, cluster, .keep_all = T)
    }
    if (print)
      print(map_res, n = 100)
    object(x)@meta.data$map_cell_type <- do.call(
      dplyr::recode,
      c(list(as.character(object(x)@meta.data$seurat_clusters)),
        nl(map_res$cluster, map_res$cell_type))
    )
    x@params$group.by <- "map_cell_type"
    p.map_marker_list <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = TRUE, pt.size = 0.5,
        group.by = "map_cell_type", cols = color_set()
        ))
    x@plots[[ 4 ]]$p.map_marker_list <- wrap(as_grob(p.map_marker_list))
    message(crayon::yellow(paste0("Use `plots(",
          as.character(substitute(x, env = parent.frame(1))),
          ", 4)$p.map_marker_list` to show result.")))
    return(x)
  })

setMethod("ids", signature = c(x = "job_seurat"),
  function(x, id = x@params$group.by, unique = T){
    ids <- object(x)@meta.data[[ id ]]
    if (unique)
      ids <- unique(ids)
    ids
  })

plot_var2000 <- function(x) {
  top20 <- head(SeuratObject::VariableFeatures(x), 20)
  p.var2000 <- e(Seurat::VariableFeaturePlot(x))
  p.var2000 <- e(Seurat::LabelPoints(p.var2000, points = top20, repel = TRUE))
  p.var2000 <- wrap(p.var2000, 12, 8)
  p.var2000
}

plot_pca.seurat <- function(x) {
  p.pca_pcComponents <- e(Seurat::VizDimLoadings(
      x, dims = 1:2, reduction = "pca",
      col = ggsci::pal_npg()(2)[2], combine = F
      ))
  p.pca_pcComponents <- lapply(p.pca_pcComponents,
    function(p) {
      p$layers[[1]]$aes_params$size <- 4
      p$layers[[1]]$aes_params$alpha <- .5
      p + theme_minimal()
    })
  p.pca_pcComponents <- patchwork::wrap_plots(p.pca_pcComponents, ncol = 2)
  p.pca_1v2 <- e(Seurat::DimPlot(x, reduction = "pca", cols = color_set()))
  p.pca_1v2 <- wrap(p.pca_1v2, 6, 5)
  p.pca_heatmap <- e(Seurat::DimHeatmap(
      x, dims = 1:10, cells = 500,
      balanced = TRUE, fast = F, combine = F
      ))
  n <- 0
  p.pca_heatmap <- lapply(p.pca_heatmap,
    function(p) {
      n <<- n + 1
      # p <- suppressMessages(p + scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF"))
      p <- p + ggtitle(paste0("PC_", n))
      p <- p + theme(plot.title = element_text(hjust = .5, face = "bold"))
      if (n < length(p.pca_heatmap))
        p <- p + theme(legend.position = "none")
      p
    })
  p.pca_heatmap <- patchwork::wrap_plots(p.pca_heatmap)
  p.pca_heatmap <- wrap(p.pca_heatmap, 12, 6)
  p.pca_rank <- pretty_elbowplot(e(Seurat::ElbowPlot(x)))
  p.pca_rank <- wrap(p.pca_rank, 4, 4)
  namel(
    p.pca_pcComponents,
    p.pca_1v2, p.pca_heatmap, p.pca_rank
  )
}

pretty_elbowplot <- function(plot) {
  plot$layers[[1]]$aes_params$size <- 5
  plot$layers[[1]]$aes_params$alpha <- .5
  plot$layers[[1]]$aes_params$colour <- "brown"
  plot
}

setMethod("getsub", signature = c(x = "job_seurat"),
  function(x, ...){
    object(x) <- e(SeuratObject:::subset.Seurat(object(x), ...))
    return(x)
  })

setMethod("active", signature = c(x = "job_seurat"),
  function(x, assay = "RNA"){
    object(x)@active.assay <- assay
    x
  })

plot_qc.seurat <- function(x) {
  require(patchwork)
  if (is(x, "job_seuratSp")) {
    nFeature <- "nFeature_Spatial"
    nCount <- "nCount_Spatial"
  } else if (is(x, "job_seurat")) {
    nFeature <- "nFeature_RNA"
    nCount <- "nCount_RNA"
  }
  x <- object(x)
  p.feature_count_mt <- e(Seurat::VlnPlot(x, features = c(nFeature, nCount,
        "percent.mt"), ncol = 3, pt.size = 0, alpha = .3, cols = color_set()))
  p.qcv1 <- e(Seurat::FeatureScatter(
      x, feature1 = nCount, feature2 = "percent.mt",
      cols = color_set()))
  p.qcv2 <- e(Seurat::FeatureScatter(
      x, feature1 = nCount, feature2 = nFeature,
      cols = color_set()))
  p.qc <- p.feature_count_mt / (p.qcv1 + p.qcv2)
  p.qc <- wrap(p.qc, 13, 10)
  p.qc
}

setMethod("vis", signature = c(x = "job_seurat"),
  function(x, group.by = x@params$group.by, pt.size = .7, palette = x$palette){
    if (is.null(palette)) {
      palette <- color_set()
    }
    p <- wrap(as_grob(e(Seurat::DimPlot(
          object(x), reduction = "umap", label = F, pt.size = pt.size,
          group.by = group.by, cols = palette
          ))), 7, 4)
    .set_lab(p, sig(x), "The", gs(group.by, "_", "-"))
  })

setMethod("focus", signature = c(x = "job_seurat"),
  function(x, features, group.by = x@params$group.by, sp = F, ...){
    if (sp) {
      return(focus(.job_seuratSp(object = object(x)), features, ...))
    }
    p.vln <- wrap(e(Seurat::VlnPlot(
        object(x), features = features, group.by = group.by,
        pt.size = 0, alpha = .3, cols = color_set()
        )))
    p.vln <- .set_lab(p.vln, sig(x), "violing plot of expression level of the genes")
    p.dim <- wrap(e(Seurat::FeaturePlot(
        object(x), features = features
        )))
    p.dim <- .set_lab(p.dim, sig(x), "dimension plot of expression level of the genes")
    namel(p.vln, p.dim)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "character"),
  function(x, ref, slot = "scale.data", ...){
    p.heatmap <- wrap(as_grob(e(Seurat::DoHeatmap(
          object(x), features = ref, raster = T, size = 3, slot = slot, ...
          ))))
    p.heatmap <- .set_lab(p.heatmap, sig(x), "heatmap show the reference genes")
    p.heatmap
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_seurat"),
  function(x, ref, use.x, use.ref, name = "cell_mapped", asIdents = T){
    matched <- match(rownames(object(x)@meta.data), rownames(object(ref)@meta.data))
    object(x)@meta.data[[name]] <- as.character(object(ref)@meta.data[[use.ref]])[matched]
    object(x)@meta.data[[name]] <- ifelse(
      is.na(object(x)@meta.data[[name]]),
      as.character(object(x)@meta.data[[use.x]]),
      object(x)@meta.data[[name]]
    )
    object(x)@meta.data[[name]] <- factor(object(x)@meta.data[[name]],
      levels = sort(unique(object(x)@meta.data[[name]]))
    )
    if (asIdents) {
      if (identical(names(object(x)@active.ident), rownames(object(x)@meta.data))) {
        SeuratObject::Idents(object(x)) <- nl(rownames(object(x)@meta.data), object(x)@meta.data[[name]], F)
      } else {
        stop("identical(names(object(x)@active.ident), rownames(object(x)@meta.data)) == F")
      }
    }
    return(x)
  })

prepare_10x <- function(target, pattern, single = F, col.gene = 1, check = T) {
  if (!single) {
    files <- list.files(target, "\\.gz$", full.names = T)
    files <- files[ grepl(pattern, files) ]
    dir <- paste0(target, "/", get_realname(files)[1])
    dir.create(dir, F)
    lapply(files,
      function(x)
        file.rename(x, paste0(dir, "/", gs(basename(x), ".*_", "")))
    )
    expected <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    alls <- list.files(dir, full.names = T)
    allnames <- basename(alls)
    if (!all(expected %in% allnames)) {
      message("Got all files: ", paste0(allnames, collapse = ", "))
      notGet <- expected[ !expected %in% allnames ]
      Get <- expected[ expected %in% allnames ]
      message("Got: ", paste0(Get, collapse = ", "), "\n",
        "Not got: ", crayon::yellow(paste0(notGet, collapse = ", ")))
      isUnExpected <- !allnames %in% expected
      if (length(which(isUnExpected)) == 1 & length(notGet) == 1) {
        file.rename(alls[ isUnExpected ], paste0(dir, "/", notGet))
        message(crayon::red("Rename: "), alls[ isUnExpected ], " as ", paste0(dir, "/", notGet))
      }
    }
    dir
  } else {
    if (!missing(pattern)) {
      stop("`pattern` not supported while single == T")
    }
    path <- dirname(target)
    name <- get_realname(target)
    dir <- paste0(path, "/", name)
    if (file.exists(dir)) {
      unlink(dir, T, T)
    }
    data <- ftibble(target)
    if (length(col.gene) > 1) {
      message("Use first element of `col.gene`, and removed the others.")
      data <- dplyr::relocate(data, dplyr::all_of(col.gene))
      data <- dplyr::select(data, -dplyr::all_of(seq_along(col.gene)[-1]))
    }
    print(data)
    if (check) {
      isCon <- usethis::ui_yeah("Continue?")
      if (!isCon) {
        stop("...")
      }  
    }
    data[[ 1 ]] <- gs(data[[ 1 ]], "\\.[0-9]*$", "")
    data <- dplyr::distinct(data, !!rlang::sym(colnames(data)[1]), .keep_all = T)
    ## as Matrix
    features <- data[[ 1 ]]
    data <- as.matrix(data[, -1])
    rownames(data) <- features
    require(Matrix)
    mtx <- as(data, "Matrix")
    DropletUtils::write10xCounts(dir, mtx)
    return(dir)
  }
}

setMethod("skel", signature = c(x = "job_seurat"),
  function(x, sig = "sr"){
    code <- c('',
      'sr <- job_seurat("")',
      'sr <- step1(sr)',
      'sr@plots$step1$p.qc',
      'sr <- step2(sr, 0, 7500, 35)',
      'sr@plots$step2$p.pca_rank',
      'sr <- step3(sr, 1:15, 1.2)',
      'sr@step <- 4L',
      'sr@plots$step3$p.umap',
      'sr <- step5(sr, 5)',
      'sr <- step6(sr, "Kidney")',
      'sr@plots$step6$p.map_scsa',
      "clear(sr)"
    )
    code <- gs(code, "\\bsr\\b", sig)
    writeLines(code)
  })

setGeneric("asjob_seurat", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_seurat"))

setMethod("mutate", signature = c(x = "job_seurat"),
  function(x, ...){
    object(x)@meta.data <- dplyr::mutate(object(x)@meta.data, ...)
    return(x)
  })

setMethod("cal_corp", signature = c(x = "job_seurat", y = "NULL"),
  function(x, y, from, to, names = NULL)
  {
    data <- object(x)[[ object(x)@active.assay ]]@data
    data <- data[ rownames(data) %in% unique(c(from, to)), ]
    anno <- data.frame(symbol = rownames(data))
    data <- data.frame(Matrix::as.matrix(data))
    .cal_corp.elist(data, anno, use = "symbol", from, to, names)
  })

setMethod("meta", signature = c(x = "job_seurat"),
  function(x){
    as_tibble(object(x)@meta.data)
  })

applyGptcelltype <- function(input, tissuename, model = c("gpt-3.5-turbo", "gpt-4"),
  topgenenumber = 10)
{
  # https://platform.openai.com/api-keys
  OPENAI_API_KEY <- getOption("OPENAI_API_KEY")
  if (is.null(OPENAI_API_KEY)) {
    API.flag <- 0L
  } else {
    API.flag <- 1L
  }
  model <- match.arg(model)
  if (is(input, "list")) {
    input <- sapply(input, paste, collapse = ',')
  } else {
    if (!is.data.frame(input)) {
      stop("`input` must be either data.frame or list.")
    }
    if (!all(c("cluster", "gene") %in% colnames(input))) {
      stop("The data.frame format of `input` must contain columns: cluster, gene.")
    }
    input <- tapply(input$gene, list(input$cluster),
      function(i) paste0(i[1:topgenenumber], collapse = ','))
  }
  if (!API.flag){
   stop('Identify cell types of ', tissuename,
     ' cells using the following markers separately for each\n row.',
     ' Only provide the cell type name. Do not show numbers before the name.\n',
     ' Some can be a mixture of multiple cell types.\n',
     paste0(names(input), ':', unlist(input), collapse = "\n"))
  } else {
    cutnum <- ceiling(length(input) / 30)
    if (cutnum > 1) {
      cid <- as.numeric(cut(1:length(input), cutnum))
    } else {
      cid <- rep(1, length(input))
    }
    allres <- sapply(1:cutnum, simplify = F,
      function(i) {
        id <- which(cid == i)
        content <- paste0(
          'Identify cell types of ', tissuename, ', ',
          ' using the following markers separately for each row.',
          ' Only provide the cell type name.',
          ' Do not show numbers before the name.\n',
          ' Some can be a mixture of multiple cell types.\n\n',
          paste(input[id], collapse = '\n')
        )
        flag <- 0L
        while (!flag) {
          k <- openai::create_chat_completion(model = model,
            messages = list(list(role = "user", content = "How are you?")),
            openai_api_key = OPENAI_API_KEY
          )
          stop_debug(k)
          res <- strsplit(k$choices[, 'message.content'], '\n')[[1]]
          if (length(res) == length(id))
            flag <- 1
        }
        names(res) <- names(input)[id]
        res
      })
    gsub(',$', '', unlist(allres))
  }
}

identify.mouseMacroPhe <- function(x, use = "scsa_cell",
  cell.name = "Macrophage",
  markers = x@tables$step5$all_markers_no_filter, top.ref = 10, ...)
{
  message("This function only support for organism of 'mouse'.\n",
    "Make sure the dataset has been subset, only retain the Macrophage cells.")
  ref <- get_data.nmt2015()
  Macrophage_M0 <- dplyr::filter(ref@object$sig.m0, FC.M1_vs_M0 < 1)$gene
  Macrophage_M1 <- dplyr::filter(ref@object$sig.m1, FC.M1_vs_M0 > 1)$gene
  Macrophage_M2 <- dplyr::filter(ref@object$sig.m2, FC.M2_vs_M0 > 1)$gene
  ref.markers <- namel(Macrophage_M0, Macrophage_M1, Macrophage_M2)
  ref.markers <- lapply(ref.markers, function(i) i[ i %in% rownames(object(x)@assays$SCT$scale.data) ])
  if (!is.null(top.ref)) {
    ref.markers <- lapply(ref.markers, head, n = top.ref)
  }
  ref.markers <- as_df.lst(ref.markers, "cell", "markers")
  lst <- scsa_annotation(x, "All", ref.markers, org = "M", res.col = "macrophage_phenotypes",
    onlyUseRefMarkers = T, ...)
  x <- lst$x
  ref.markers <- .set_lab(as_tibble(ref.markers), sig(x), "the markers for Macrophage phenotypes annotation")
  x$p.macrophage <- lst$p.map_scsa
  x$t.macrophage_ref <- ref.markers
  x$p.macrophage_hp_ref <- map(x, ref.markers$markers)
  ref@object <- NULL
  .add_internal_job(ref)
  return(x)
}

scsa_annotation <- function(
  x, tissue, ref.markers = NULL, filter.p = 0.01, filter.fc = .5,
  org = c("Human", "Mouse"),
  cmd = pg("scsa"), db = pg("scsa_db"), res.col = "scsa_annotation",
  onlyUseRefMarkers = F)
{
  if (grpl(tissue, "(?<!\\\\)\\s", perl = T)) {
    tissue <- gs(tissue, "\\s", "\\\\ ")
  }
  org <- match.arg(org)
  prop <- as.list(prop.table(table(grpl(rownames(object(x)), "^[A-Z][a-z]"))))
  if (!is.null(prop[[ "TRUE" ]])) {
    if (prop[[ "TRUE" ]] > .5 && org == "Human") {
      isThat <- usethis::ui_yeah("Match org of 'Human', but most of genes is Mouse like symbol, continue?")
      if (!isThat) {
        stop("...")
      }
    }  
  }
  if (!is.null(ref.markers)) {
    .check_columns(ref.markers, c("cell", "markers"))
    ref.markers <- dplyr::relocate(ref.markers, cell, markers)
    ref.markers_file <- tempfile("ref.markers_file", fileext = ".tsv")
    write_tsv(ref.markers, ref.markers_file, col.names = F)
    x@params$ref.markers_file <- ref.markers_file
    ref.markers.cmd <- paste0(" -M ", ref.markers_file)
  } else {
    ref.markers.cmd <- ""
  }
  marker_file <- tempfile("marker_file", fileext = ".csv")
  result_file <- tempfile("result_file")
  all_markers <- dplyr::rename(x@tables$step5$all_markers_no_filter, avg_logFC = avg_log2FC)
  all_markers <- dplyr::select(all_markers, -rownames)
  all_markers <- dplyr::mutate(all_markers, gene = gs(gene, "\\.[0-9]*", ""))
  data.table::fwrite(all_markers, marker_file)
  cli::cli_alert_info(cmd)
  cdRun(cmd,
    " -s seurat", " -i ", marker_file,
    " -k ", tissue, " -d ", db, " ", ref.markers.cmd,
    " -p ", filter.p, " -f ", filter.fc,
    " -E -g ", org, " -m txt",
    if (onlyUseRefMarkers) " -N " else NULL,
    " -o ", result_file, " > /tmp/log_scsa.txt"
  )
  scsa_res <- dplyr::rename_all(ftibble(result_file), make.names)
  scsa_res <- scsa_res_all <- dplyr::arrange(scsa_res, Cluster, dplyr::desc(Z.score))
  ## add into seurat object
  scsa_res <- dplyr::distinct(scsa_res, Cluster, .keep_all = T)
  clusters <- object(x)@meta.data$seurat_clusters
  cell_types <- scsa_res$Cell.Type[match(clusters, scsa_res$Cluster)]
  cell_types <- ifelse(is.na(cell_types), "Unknown", cell_types)
  object(x)@meta.data[[ res.col ]] <- factor(cell_types,
    levels = c(unique(scsa_res$Cell.Type), "Unknown"))
  ## plot
  p.map_scsa <- e(Seurat::DimPlot(
      object(x), reduction = "umap", label = F, pt.size = .7,
      group.by = res.col, cols = color_set()
      ))
  p.map_scsa <- wrap(as_grob(p.map_scsa), 7, 4)
  p.map_scsa <- .set_lab(p.map_scsa, sig(x), "SCSA", "Cell type annotation")
  .add_internal_job(.job(method = "`SCSA` (python) used for cell type annotation",
      cite = "[@ScsaACellTyCaoY2020]"))
  namel(x, p.map_scsa, res.col, scsa_res_all)
}

matchCellMarkers <- function(lst.markers, ref, least = 2) {
  lst <- vapply(lst.markers, FUN.VALUE = logical(1),
    function(x) {
      x <- as.list(table(ref %in% x))
      x[[ "TRUE" ]] >= least
    }
  )
  lst
}

parse_GPTfeedback <- function(feedback) {
  res <- gs(feedback, "^[0-9]+\\. |\\.$", "")
  res <- strsplit(res, "\\. ")
  cells <- vapply(res, function(x) strsplit(x[1], "; ")[[1]][1], character(1))
  markers <- lapply(res, function(x) strsplit(x[2], "; "))
  markers <- markers[ order(cells) ]
  markers <- unique(gs(unlist(markers), "\\s", ""))
  namel(markers, cells)
}

prepare_GPTmessage_for_celltypes <- function(tissue, marker_list, n = 10, toClipboard = T)
{
  message <- glue::glue("Identify cell types of {tissue} cells using the following markers separately for each row. Only provide the cell type name (for each row). Show numbers before the name. Some can be a mixture of multiple cell types (separated by '; ', end with '. '). Then, provide at least 1 and at most 3 classical markers that distinguish the cell type from other cells (separated by '; ').  (e.g., 1. X Cell; Y Cell. Marker1; Marker2; Marker3)")
  message("`marker_list` should be table output from seurat (column: cluster, gene)")
  marker_list <- split(marker_list$gene, marker_list$cluster)
  ncluster <- length(marker_list)
  if (!is.null(n)) {
    marker_list <- lapply(marker_list, head, n = n)
  }
  markers <- vapply(marker_list, paste0, character(1), collapse = ",")
  markers <- paste0(seq_along(markers), ". ", markers)
  markers <- paste(markers, collapse = "\n")
  query <- paste0(message, "\n", markers)
  if (toClipboard) {
    gett(query)
  }
  return(list(query = query, message = message, ncluster = ncluster))
}
