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
  function(x, sets = "H", pattern = NULL, name = pattern, gather = TRUE,
    assay = SeuratObject::DefaultAssay(object(x)), ...)
  {
    mtx <- Seurat::GetAssayData(object(x), assay = assay, layer = "data")
    if (is.null(mtx) || is.null(rownames(mtx))) {
      stop('is.null(mtx) || is.null(rownames(mtx)).')
    }
    if (!is(sets, "feature") && is.character(sets) && length(sets) == 1) {
      j <- .set_msig_db(.job(), sets)
      mode <- sets
      sets <- as_feature(
        lapply(split(j$db_anno$gene_symbol, j$db_anno$gs_name), unique),
        glue::glue("MSigDB {mode} 基因集")
      )
      methodAdd_onExit("x", meth(j)$step0)
      if (!is.null(pattern)) {
        sets <- sets[ grp(names(sets), pattern, TRUE) ]
        methodAdd_onExit("x", "在 {mode} 基因集中获取与 {pattern} 相关的基因。")
        if (gather) {
          genes <- unique(resolve_feature(sets))
          sets <- as_feature(setNames(list(genes), name), glue::glue("{name} 相关基因集"))
          methodAdd_onExit("x", "对其去重、合并，作为 {name} 基因集用于 AUCell 评分分析。{name} 相关基因集共包含 {length(genes)} 个基因。")
        } else {
          methodAdd_onExit("x", "将 {name} 基因集用于 AUCell 评分分析。该基因集共包含 {length(sets)} 个子集，各子集包含基因统计为：{try_snap(sets)}。")
          snap(sets) <- glue::glue("{name} 相关基因集")
        }
      }
    } else {
      methodAdd_onExit("x", "使用{snap(sets)}作为 AUCell 输入。")
    }
    if (!is(sets, "feature")) {
      stop('!is(sets, "feature").')
    }
    pr <- params(x)
    x <- job_aucell(mtx, sets, ...)
    x@params <- append(x@params, pr)
    x$.feature_genesSets <- sets
    x$name <- name
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
  x <- methodAdd(x, "AUCell 是一种基于单细胞转录组数据评估基因集活性的分析方法，其主要目的是在单细胞分辨率下量化预定义基因集（如信号通路、转录因子靶基因集或细胞状态特征基因）的活跃程度。该方法通过对每个细胞内基因表达进行排序，并计算目标基因集在高表达基因中的富集面积（AUC score），从而评估该基因集在不同细胞中的相对活性。")
  x <- methodAdd(x, "以 R 包 `AUCell` ⟦pkgInfo('AUCell')⟧ {cite_show('SCENIC_single_Aibar_2017')} 识别单细胞数据集的基因集调控活性。")
  x$sets <- sets
  x$gids <- gids
  return(x)
}

setMethod("step0", signature = c(x = "job_aucell"),
  function(x){
    step_message("Prepare your data with function `job_aucell`.")
  })

setMethod("step1", signature = c(x = "job_aucell"),
  function(x, workers = NULL, group.by = x$group.by,
    rerun = FALSE, fun_name = function(x) s(x, "^HALLMARK_", ""))
  {
    step_message("Running...")
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L,
        {
          x <- step1(x, workers = "{workers}")
        }
      )
      return(x)
    }
    fun_show <- function(string) stringr::str_wrap(gs(string, "_", " "), 20)
    fun_aucell <- function(...) {
      if (!is.null(workers)) {
        workers <- e(BiocParallel::MulticoreParam(workers))
      }
      e(AUCell::AUCell_run(object(x), x$sets, BPPARAM = workers))
    }
    res_aucell <- expect_local_data(
      "tmp", "AUcell", fun_aucell,
      list(colnames(object(x)), sig(x), names(x$sets), x$gids), 
      rerun = rerun
    )
    x$res_aucell <- e(AUCell::getAUC(res_aucell))
    if (!is.null(fun_name)) {
      rownames(x$res_aucell) <- fun_name(rownames(x$res_aucell))
    }
    x <- methodAdd(x, "采用 AUCell 算法对预定义功能基因集在单细胞水平进行活性评分，并提取各细胞对应的 AUC（Area Under the Curve）值。")
    if (!is.null(group.by)) {
      metadata <- x$metadata
      fun_mean <- function(...) {
        auc <- .get_auc_from_job_aucell(x)
        data <- cbind(
          metadata[, group.by, drop = FALSE],
          as.data.frame(auc)
        )
        require(data.table)
        data <- as.data.table(data)
        data <- data[ ,
          lapply(.SD, mean),
          by = group.by,
          .SDcols = setdiff(names(data), group.by)
          ]
        tibble::as_tibble(data)
      }
      data <- expect_local_data(
        "tmp", "aucell_mean", fun_mean,
        list(rownames(x$res_aucell), metadata$cell, group.by, x$gids),
        rerun = rerun
      )
      x <- methodAdd(x, "基于细胞注释信息对同一细胞群内所有细胞的 AUC 值取平均，以评估不同细胞类型的整体功能状态。")
      x$res_aucell_mean <- data
      if (ncol(data) < 21L) {
        layout <- wrap_layout(NULL, ncol(data) - 1L, 3)
        data <- tidyr::pivot_longer(
          data, -!!rlang::sym(group.by), 
          names_to = "Function", values_to = "Activity"
        )
        if (!is.null(fun_show)) {
          data <- dplyr::mutate(data, Function = fun_show(Function))
        }
        p.aucell_mean <- ggplot(data, aes(x = reorder(!!rlang::sym(group.by), Activity), y = Activity)) +
          geom_col() +
          facet_wrap(~ Function, ncol = layout$ncol) +
          labs(x = "Cell types", y = "Activity") +
          coord_flip() +
          theme_minimal()
        p.aucell_mean <- set_lab_legend(
          add(layout, p.aucell_mean),
          glue::glue("{x@sig} Mean AUCell Activity"),
          glue::glue("各细胞类型评价 AUCell 功能活性|||每个分面代表一个独立的功能通路或生物学过程（Function），横坐标表示不同细胞群体（按平均活性值排序），纵坐标表示该群体的平均 AUCell 活性评分（Activity）。")
        )
        x <- snapAdd(
          x, "对于每个功能基因集，计算对应细胞群体的平均活性分数，并以分面柱状图形式展示{aref(p.aucell_mean)}。"
        )
        x <- plotsAdd(x, p.aucell_mean)
      }
    }
    return(x)
  })

setMethod("quantile", signature = c(x = "job_aucell"),
  function(x, cols = NULL, cut = .75,
    gather = c("merge", "intersect", "respective"), 
    name = x$name, group.by = x$group.by, ...)
  {
    gather <- match.arg(gather)
    if (is.null(cols)) {
      auc <- .get_auc_from_job_aucell(x)
      cols <- colnames(auc)
    }
    data <- x$res_aucell_mean
    if (is.null(data)) {
      stop('is.null(data), no `x$res_aucell_mean` data.')
    }
    celltypes <- quantile(
      data.frame(data), cols, get = group.by, cut = cut, ...
    )
    if (gather == "intersect") {
      if (length(cols) > 1) {
        celltypes <- ins(lst = celltypes)
      } else {
        celltypes <- unlist(celltypes)
      }
      snap <- glue::glue("{bind(cols)} 基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 的细胞类型")
      as_feature(as.character(celltypes), snap, nature = "cell")
    } else if (gather == "respective") {
      snap <- glue::glue("基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 的细胞类型")
      as_feature(
        lapply(celltypes, as.character), snap, nature = "cell"
      )
    } else if (gather == "merge") {
      snap <- glue::glue("{name} 各基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 合并后的细胞类型")
      celltypes <- unique(as.character(unlist(celltypes)))
      as_feature(setNames(list(celltypes), name), snap, nature = "cell")
    }
  })

setMethod("step2", signature = c(x = "job_aucell"),
  function(x, group.by = "seurat_clusters")
  {
    step_message("Annotate for clusters.")
    auc <- .get_auc_from_job_aucell(x)
    require(data.table)
    dt <- as.data.table(auc)
    dt$cluster <- x$metadata[[ group.by ]]
    if (is.null(dt$cluster)) {
      stop('is.null(dt$cluster), no value in of column `group.by` in metadata.')
    }
    cluster_mean <- dt[ ,
      lapply(.SD, mean),
      by = cluster,
      .SDcols = setdiff(names(dt), "cluster")
      ]
    mean_scaled <- t(scale(t(cluster_mean[, -1])))
    x <- methodAdd(x, "依据细胞分群结果 ({group.by}) 对同一亚群内细胞的 AUC 分数取平均，以获得各细胞亚群的整体功能活性特征。为消除不同基因集之间评分尺度差异，对每个亚群的平均 AUC 矩阵按 Cluster 进行 Z-score 标准化处理。")
    mean_scaled <- tibble::as_tibble(
      mean_scaled
    )
    mean_scaled <- dplyr::mutate(
      mean_scaled, cluster = cluster_mean$cluster, 
      .before = 1
    )
    mean_scaled <- tidyr::pivot_longer(
      mean_scaled, -cluster,
      names_to = "Function", values_to = "Activity"
    )
    annotation <- dplyr::group_by(mean_scaled, cluster)
    annotation <- dplyr::summarise(
      annotation, Function = Function[ which.max(Activity) ]
    )
    x <- methodAdd(x, "对于每个细胞亚群，进一步筛选其标准化活性值最高的功能通路，并将该通路作为该亚群的主要功能注释（Annotation）。")
    data <- dplyr::filter(
      mean_scaled, Function %in% !!annotation$Function
    )
    data <- map(
      data, "cluster", annotation, "cluster", "Function", col = "Annotation"
    )
    args <- list(
      .data = data, .row = quote(Function), .column = quote(cluster),
      .value = quote(Activity), group_by = quote(Annotation),
      cluster_columns = TRUE, column_names_rot = 45,
      cluster_rows = TRUE,
      row_names_max_width = grobWidth(textGrob(data$Function, gpar(fontsize = 10, fontface = 1)))
    )
    rm(dt, auc)
    p.hp <- wrap_scale_heatmap(
      funPlot(heatmap_with_group, args),
      data$cluster, data$Function, pre_width = 6
    )
    p.hp <- set_lab_legend(
      p.hp,
      glue::glue("{x@sig} Cluster functional enrichment score heatmap"),
      glue::glue("功能富集得分热图|||热图展示不同细胞亚群与代表性功能状态之间的对应关系，其中横坐标细胞亚群，纵坐标表示表示功能基因集（Function），颜色梯度表示标准化后的相对活性强弱（Activity），暖色代表该亚群中该功能相对激活，冷色代表相对低活性。热图的 Cluster 对应有注释类型。热图行列均基于功能活性模式进行层次聚类。")
    )
    types <- unique(annotation$Function)
    x <- snapAdd(
      x, "如图{aref(p.hp)} (热图仅展示有 cluster 注释的功能的活性)，⟦mark$red('AUCell 亚群功能富集一共注释了 {length(types)} 种类型的功能，分别为：{bind(types)}')⟧。"
    )
    x$metadata <- map(
      x$metadata, group.by, annotation, "cluster", "Function",
      col = "AUCell_Function"
    )
    x <- plotsAdd(x, p.hp)
    return(x)
  })

setMethod("step3", signature = c(x = "job_aucell"),
  function(x, use.trait, data_trait = x$metadata, rerun = FALSE)
  {
    step_message("Correlation with trait data.")
    auc <- .get_auc_from_job_aucell(x)
    types <- unique(x$metadata[[ "AUCell_Function" ]])
    data_trait <- data_trait[, use.trait, drop = FALSE]
    data_aucell <- as.data.frame(auc[, colnames(auc) %in% types])
    fun_cor <- function(...) {
      cli::cli_alert_info("safe_fortify_cor")
      safe_fortify_cor(data_trait, data_aucell)
    }
    x$cor_trait_aucell <- expect_local_data(
      "tmp", "aucell_trait_activity_cor", fun_cor, rerun = rerun,
      list(
        x$metadata$cell, colnames(data_aucell), 
        colnames(data_trait)
      )
    )
    snap_cor <- .stat_ggcor_table_list(
      x$cor_trait_aucell, "Trait", "Function"
    )
    p.cor_trait_aucell <- .ggcor_add_general_style(ggcor::quickcor(x$cor_trait_aucell))
    p.cor_trait_aucell <- set_lab_legend(
      wrap_scale_heatmap(p.cor_trait_aucell, length(use.trait), length(types), raw = FALSE),
      glue::glue("{x@sig} trait correlation with AUCell"),
      glue::glue("{bind(use.trait)} 与 AUCell 活性关联分析热图|||热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。P 值以 * 标注 ({.md_p_significant})。")
    )
    x <- snapAdd(x, "对表型 ({bind(use.trait)}) 与 AUCell Function 活性之间关联分析，如图{aref(p.cor_trait_aucell)}，{snap_cor}")
    x <- plotsAdd(x, p.cor_trait_aucell)
    return(x)
  })

.get_auc_from_job_aucell <- function(x) {
  if (!identical(colnames(x$res_aucell), x$metadata$cell)) {
    stop('!identical(colnames(x$res_aucell), x$metadata$cell).')
  }
  t(x$res_aucell)
}

setMethod("map", signature = c(x = "job_aucell", ref = "job_seurat"),
  function(x, ref, use.trait = NULL, use.function = NULL,
    group.by = ref$group.by, pal = NULL, .name = "seurat")
  {
    fun_show <- function(string) stringr::str_wrap(gs(string, "_", " "), 20)
    fun_rename_title <- function(lst) {
      lst$title <- fun_show(lst$title)
      lst
    }
    if (x@step == 1L) {
      if (is.null(use.function)) {
        auc <- .get_auc_from_job_aucell(x)
        use.function <- colnames(auc)
      }
      meta <- dplyr::select(x$metadata, cell)
      meta <- cbind(meta, auc)
      meta <- data.frame(meta[, -1L, drop = FALSE], row.names = meta$cell)
      object(ref) <- SeuratObject::AddMetaData(object(ref), meta)
      layout <- wrap_layout(NULL, length(use.function))
      ps.map <- e(Seurat::FeaturePlot(object(ref), 
          features = use.function, combine = FALSE
          ))
      if (!is(ps.map, "list")) {
        ps.map <- list(ps.map)
      }
      ps.map <- lapply(ps.map, 
        function(x) {
          x + theme(plot.title = element_text(face = "plain", size = 10))
        })
      if (!is.null(fun_show)) {
        ps.map <- lapply(
          ps.map, .set_ggplot_content,
          fun = fun_rename_title,
          slot = "labels"
        )
      }
      p.map <- add(layout, ps.map, TRUE)
      p.map <- set_lab_legend(p.map,
        glue::glue("{x@sig} AUCell Activity UMAP mapping"),
        glue::glue("AUCell 功能活性 UMAP 图||| {bind(use.function)} 的 AUCell 功能活性 UMAP 图。")
      )
      x[[ glue::glue("map_{.name}") ]] <- namel(p.map, metadata = meta)
    } else if (x@step > 1L) {
      meta <- x$metadata[, !duplicated(colnames(x$metadata))]
      meta <- dplyr::select(meta, cell, AUCell_Function)
      col_r_trait <- NULL
      col_trait <- NULL
      if (!is.null(use.trait)) {
        if (is.null(x$cor_trait_aucell)) {
          stop('is.null(x$cor_trait_aucell), but !is.null(use.trait)')
        }
        data_cor <- dplyr::filter(
          x$cor_trait_aucell, .row.names %in% use.trait
        )
        lst_cor <- split(data_cor, ~ .row.names)
        col_r_trait <- glue::glue("r_{use.trait}")
        for (i in seq_along(lst_cor)) {
          meta <- map(
            meta, "AUCell_Function", lst_cor[[i]], ".col.names", 
            "r", col = col_r_trait[ i ]
          )
        }
        allAvai <- c(colnames(meta(ref)), colnames(meta))
        col_trait <- use.trait[ use.trait %in% allAvai ]
        if (any(isNot <- !use.trait %in% col_trait)) {
          warning("Can not got trait from `ref`: ", bind(use.trait[ isNot ]))
        }
      }
      meta <- data.frame(meta[, -1L, drop = FALSE], row.names = meta$cell)
      object(ref) <- SeuratObject::AddMetaData(object(ref), meta)
      group <- c("AUCell_Function", group.by, col_trait)
      ps.map <- e(Seurat::DimPlot(
          object(ref), pt.size = if (dim(object(ref))[2] > 30000) .3 else .5,
          group.by = group,
          cols = color_set(), combine = FALSE
          ))
      if (!is.null(pal) && !is.null(use.trait)) {
        whichTrait <- which(group %in% use.trait)
        for (i in whichTrait) {
          ps.map[[i]] <- ps.map[[i]] + scale_color_manual(values = pal)
        }
      }
      if (!is.null(col_r_trait)) {
        ps2.map <- e(Seurat::FeaturePlot(object(ref), 
            features = col_r_trait, combine = FALSE
            ))
        if (!is(ps2.map, "list")) {
          ps2.map <- list(ps2.map)
        }
        for (i in seq_along(col_r_trait)) {
          ps2.map[[i]] <- ps2.map[[i]] + .scale_for_cor_palette("color")
        }
        group <- c(group, col_r_trait)
        ps.map <- c(ps.map, ps2.map)
      }
      legend_ex <- ""
      if (!is.null(use.trait)) {
        legend_ex <- glue::glue("其中，{bind(col_r_trait)} 对应为 AUCell_Function 与 {bind(use.trait)} 的关联分析的相关系数。")
      }
      layout <- z7(wrap_layout(NULL, length(group), ncol = 2), 1.7, 1)
      p.map <- add(layout, ps.map, TRUE)
      p.map <- set_lab_legend(p.map,
        glue::glue("{x@sig} Cell Function UMAP mapping"),
        glue::glue("AUCell 功能活性 UMAP 图|||依次对应为 {bind(group)} 的 UMAP 图。{legend_ex}")
      )
      x[[ glue::glue("map_{.name}") ]] <- namel(p.map, metadata = meta)
    }
    return(x)
  })

setMethod("clear", signature = c(x = "job_aucell"),
  function(x, save = TRUE, lite = TRUE, suffix = NULL, name = substitute(x, parent.frame(1)))
  {
    eval(name)
    if (save) {
      callNextMethod(
        x, save = save, lite = FALSE, suffix = suffix, name = name
      )
    }
    object(x) <- NULL
    x$res_aucell <- NULL
    if (lite) {
      callNextMethod(
        x, save = FALSE, lite = TRUE, suffix = suffix, name = name
      )
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_aucell"),
  function(x, ref, type = "AUC", scale = FALSE){
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    if (type == "AUC") {
      res <- t(ref$res_aucell)
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

