# ==========================================================================
# workflow of hdwgcna
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_hdwgcna <- setClass("job_hdwgcna", 
  contains = c("job"),
  prototype = prototype(
    pg = "hdwgcna",
    info = c(
      "https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html",
      "https://smorabit.github.io/hdWGCNA/articles/module_trait_correlation.html"
    ),
    cite = "",
    method = "",
    tag = "hdwgcna",
    analysis = "hdWGCNA 加权基因共表达分析"
    ))

setGeneric("asjob_hdwgcna",
  function(x, ...) standardGeneric("asjob_hdwgcna"))

setMethod("asjob_hdwgcna", signature = c(x = "job_seurat"),
  function(x, workers = 5, assay = SeuratObject::DefaultAssay(object(x)),
    name = "wgcna")
  {
    object <- object(x)
    if (FALSE) {
      object <- e(Seurat::DietSeurat(
          object, assays = assay, layers = c("counts", "data", "scale.data")
        ))
    }
    e(WGCNA::enableWGCNAThreads(workers))
    object <- e(hdWGCNA::SetupForWGCNA(
      object, gene_select = "fraction",
      fraction = 0.05, wgcna_name = name
    ))
    SeuratObject::DefaultAssay(object) <- assay
    pr <- params(x)
    x <- .job_hdwgcna(object = object)
    x <- methodAdd(x, "**hdWGCNA** 是一种针对单细胞或高维转录组数据构建加权基因共表达网络的分析方法，其主要目的是在复杂数据中识别具有协同表达特征的基因模块，并解析其与细胞类型、状态或表型之间的关联。")
    x <- methodAdd(x, "为了系统鉴定关键细胞功能异质性的核心共表达模块及靶点基因，以 R 包 `hdWGCNA` ⟦pkgInfo('hdWGCNA')⟧对 Seurat 处理过的单细胞数据集展开加权基因共表达网络分析。")
    x@params <- append(x@params, pr)
    return(x)
  })

setMethod("step0", signature = c(x = "job_hdwgcna"),
  function(x){
    step_message("Prepare your data with function `job_hdwgcna`.")
  })

setMethod("step1", signature = c(x = "job_hdwgcna"),
  function(x, min_cells = 100, k = 25, max_shared = 15,
    reduction = NULL, group.by = x$group.by, auto = FALSE, 
    ..., debug = FALSE)
  {
    step_message("Metacells.")
    if (is.null(reduction)) {
      reduction <- names(object(x)@reductions)
      reduction <- tail(reduction[ reduction != "umap" ], n = 1)
      message(glue::glue("Use reduction: {reduction}"))
    }
    if (auto) {
      params <- .auto_metacell_params(object(x), group.by, ...)
      k <- params$k
      max_shared <- params$max_shared
      min_cells <- params$min_cells
      x <- methodAdd(x, "{params$snap}")
    } else {
      x <- methodAdd(
        x, "使用 ConstructMetacells 函数，设置 k = {k}, 最少细胞数量为 {min_cells}，基于 K 近邻算法将相似的细胞聚合为元细胞（metacells）。"
      )
    }
    if (!debug) {
      object(x) <- e(hdWGCNA::MetacellsByGroups(
          seurat_obj = object(x),
          group.by = group.by, ident.group = group.by,
          reduction = 'HarmonyIntegration', k = k, 
          min_cells = min_cells, max_shared = max_shared,
          verbose = TRUE
          ))
      object(x) <- e(hdWGCNA::NormalizeMetacells(object(x)))
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_hdwgcna"),
  function(x, celltypes = NULL, cut.r = .8, group.by = x$group.by, debug = FALSE){
    step_message("Soft power.")
    if (is.null(celltypes)) {
      celltypes <- unique(object(x)@meta.data[[x$group.by]])
    }
    x$celltypes <- celltypes
    ncells <- nrow(dplyr::filter(object(x)@meta.data, !!rlang::sym(x$group.by) %in% celltypes))
    if (!debug) {
      object(x) <- e(hdWGCNA::SetDatExpr(
          object(x), 
          group_name = celltypes, group.by = group.by,
          assay = SeuratObject::DefaultAssay(object(x))
        ))
      object(x) <- e(hdWGCNA::TestSoftPowers(
        object(x), networkType = "signed"
      ))
    }
    x <- methodAdd(x, "以 `SetDatExpr` 选择 {bind(celltypes)} 的表达矩阵为输入数据 (共包含 {ncells} 个细胞)。")
    x$power_table <- e(hdWGCNA::GetPowerTable(object(x)))
    x$use.power <- min(dplyr::filter(x$power_table, SFT.R.sq >= !!cut.r & Power > 3)$Power)
    message(glue::glue("Select power to plot: {x$use.power}"))
    p.sft <- patchwork::wrap_plots(
      e(hdWGCNA::PlotSoftPowers(object(x), x$use.power)), ncol = 2
    )
    p.sft <- set_lab_legend(
      wrap(p.sft, 7, 6),
      glue::glue("{x@sig} hdWGCNA soft power"),
      glue::glue("hdWGCNA 软阈值（soft power）筛选诊断图|||四个子图分别展示对应阈值下的无标度拓扑拟合指数（Scale-free Topology Model Fit）、平均连接度（Mean Connectivity）、中位连接度（Median Connectivity）及最大连接度（Max Connectivity）。每个点代表一个候选阈值，其数值标注在点旁；虚线用于指示经验阈值或参考标准。该图通过同时考察网络的无标度特性与连接性变化，帮助理解不同软阈值对共表达网络结构的影响，从而为后续网络构建提供依据。")
    )
    x <- methodAdd(x, "通过 pickSoftThreshold 函数计算不同软阈值下的无尺度网络拟合指数 R² 和平均连接度，⟦mark$blue('选择使得 R² 首次达到 {cut.r} 以上的最小 β 值 (即 β = {x$use.power}) 作为软阈值')⟧，确保网络符合无尺度拓扑结构{aref(p.sft)}。")
    x <- plotsAdd(x, p.sft)
    return(x)
  })

setMethod("step3", signature = c(x = "job_hdwgcna"),
  function(x, min.gene = 50, cut.height = .2, debug = FALSE){
    step_message("Network")
    message(glue::glue("Use power: {x$use.power}"))
    x$dir_cache <- create_job_cache_dir(x)
    if (!debug) {
      object(x) <- e(hdWGCNA::ConstructNetwork(
          object(x), x$use.power, overwrite_tom = TRUE,
          tom_outdir = x$dir_cache,
          minModuleSize = min.gene,
          mergeCutHeight = cut.height,
          tom_name = bind(x$celltypes, co = "_")
          ))
    }
    wgcna_name <- object(x)@misc$active_wgcna
    x$modules <- hdWGCNA::GetModules(object(x), wgcna_name)
    x$nm <- length(unique(dplyr::filter(x$modules, module != "grey")$module))
    p.dg <- funPlot(WGCNA::plotDendroAndColors,
      list(
        dendro = hdWGCNA::GetNetworkData(object(x), wgcna_name)$dendrograms[[1]],
        colors = as.character(x$modules$color),
        groupLabels = "Module colors", dendroLabels = FALSE,
        hang = .03, addGuide = TRUE, guideHang = .05
        ))
    p.dg <- set_lab_legend(
      wrap(p.dg, 7, 4),
      glue::glue("{x@sig} Gene co expression network module"),
      glue::glue("hdWGCNA 基因共表达网络模块|||上方为基于基因表达相似性进行层次聚类得到的树状结构（dendrogram），每个分支代表一组表达模式相近的基因；下方的彩色条带对应动态剪切（dynamic tree cut）后识别的不同共表达模块，不同颜色表示不同模块。")
    )
    x <- methodAdd(x, "使用 ConstructNetwork 函数，利用 TOM 矩阵进行层次聚类，采用动态树切割算法识别初始基因模块，设定最小模块基因数为 {min.gene} 以过滤过小模块，并通过合并模块特征基因（ME）相关性高于 {cut.height} 的高度相似模块，⟦mark$red('最终获得 {x$nm} 个稳定的基因模块')⟧{aref(p.dg)}。")
    x <- plotsAdd(x, p.dg)
    return(x)
  })

setMethod("step4", signature = c(x = "job_hdwgcna"),
  function(x, sample.by = "orig.ident", debug = FALSE)
  {
    step_message("Module eigen")
    if (!debug) {
      object(x) <- e(Seurat::ScaleData(
        object(x), features = Seurat::VariableFeatures(object(x))
      ))
      object(x) <- e(hdWGCNA::ModuleEigengenes(
        object(x), group.by.vars = sample.by
      ))
    }
    x <- methodAdd(x, "以 ModuleEigengenes 函数计算模块特征基因 (Module Eigengenes，MEs，如经过批次矫正，则为 hMEs)。")
    x$hMEs <- hMEs <- e(hdWGCNA::GetMEs(object(x)))
    hMEs <- dplyr::select(hMEs, -grey)
    x$cor_hMEs <- cor_hMEs <- safe_fortify_cor(hMEs)
    snap_corhMEs <- .stat_ggcor_table_list(
      cor_hMEs, "MEs", "MEs"
    )
    p.cor_hMEs <- .ggcor_add_general_style(ggcor::quickcor(cor_hMEs))
    p.cor_hMEs <- set_lab_legend(
      wrap_scale_heatmap(p.cor_hMEs, x$nm, x$nm, raw = FALSE),
      glue::glue("{x@sig} "),
      glue::glue("模块间相关性热图|||热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。P 值以 * 标注 ({.md_p_significant})。")
    )
    x <- plotsAdd(x, p.cor_hMEs)
    x <- snapAdd(x, "对 hMEs 之间关联分析，如图{aref(p.cor_hMEs)}，{snap_corhMEs}")
    return(x)
  })

setMethod("step5", signature = c(x = "job_hdwgcna"),
  function(x, global = TRUE, debug = FALSE){
    step_message("Module membership.")
    if (!debug) {
      object(x) <- e(
        hdWGCNA::ModuleConnectivity(
          object(x), group.by = x$group.by,
          group_name = x$celltypes
          ))
    }
    x <- methodAdd(x, "以 `ModuleConnectivity` 计算基因与模块特征基因 (hMEs) 之间的相关性 (即，WGCNA 中的 Module Membership) 得到 kME 值，以该值代表 Module Membership (MM) (算法有所不同，以 hdWGCNA 的指南为依据，筛选 Hub Genes 时不对 MM 取绝对值)。")
    modules <- e(hdWGCNA::GetModules(object(x)))
    modules <- dplyr::filter(tibble::as_tibble(modules), module != "grey")
    x$modules <- dplyr::mutate(modules, module = droplevels(module))
    x$name_modules <- levels(x$modules$module)
    layout <- wrap_layout(NULL, x$nm, 2)
    p.kme <- e(hdWGCNA::PlotKMEs(object(x), ncol = layout$ncol))
    # hdWGCNA::GetHubGenes
    p.kme <- set_lab_legend(
      add(layout, p.kme),
      glue::glue("{x@sig} module membership ranking"),
      glue::glue(
        "模块特征基因相关性（module membership）排序|||每个小图对应一个共表达模块，横轴为模块内基因，纵轴为对应基因的 kME 值，反映其与模块特征基因的相关程度；右侧标注为部分代表性高 kME 基因。"
      )
    )
    object <- object(x)
    object@meta.data <- cbind(object@meta.data, x$hMEs)
    snap_ex <- ""
    if (global) {
      p.umap_hMEs <- e(hdWGCNA::ModuleFeaturePlot(
        object, features = "hMEs", order = TRUE
      ))
    } else {
      celltypes <- x$celltypes
      subObj <- SeuratObject:::subset.Seurat(object, !!rlang::sym(x$group.by) %in% !!celltypes)
      p.umap_hMEs <- e(Seurat::FeaturePlot(
        subObj, features = x$name_modules, combine = FALSE
      ))
      snap_ex <- glue::glue("在 {celltypes} 中的")
    }
    p.umap_hMEs <- set_lab_legend(
      add(layout, p.umap_hMEs),
      glue::glue("{x@sig} hMEs in UMAP"),
      glue::glue("模块特征基因（module eigengene）{snap_ex} UMAP 分布图|||每个点代表一个细胞，点颜色表示该模块的特征基因 (hMEs) 数值。")
    )
    p.dot_hMEs <- e(
      Seurat::DotPlot(object, features = x$name_modules, group.by = x$group.by)
    )
    p.dot_hMEs <- p.dot_hMEs + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
    p.dot_hMEs <- set_lab_legend(
      wrap_scale(
        p.dot_hMEs, length(unique(object@meta.data[[x$group.by]])), x$nm,
        pre_width = 4.5, pre_height = 3.5
      ),
      glue::glue("{x@sig} dotplot of hMEs"),
      glue::glue("模块特征基因（hMEs）在不同细胞类型中表达水平与表达比例|||每个圆点代表一个模块在某一细胞类型中的表达情况，点的大小表示该模块在该细胞类型中表达的细胞百分比，点的颜色表示平均表达量（Average Expression）。")
    )
    x <- plotsAdd(x, p.kme, p.umap_hMEs, p.dot_hMEs)
    return(x)
  })

setMethod("step6", signature = c(x = "job_hdwgcna"),
  function(x, use.trait = c("group"),
    celltypes = x$celltypes,
    group.by = x$group.by, debug = FALSE, show = 5)
  {
    step_message("Module Trait Correlation")
    if (any(which <- use.trait == "group")) {
      # for debug, hdWGCNA During the calculation process, a "group" column will also be generated
      use.trait[which] <- "Group"
      object(x)@meta.data <- dplyr::mutate(object(x)@meta.data, Group = group)
    }
    for (i in use.trait) {
      if (is.character(object(x)@meta.data[[i]])) {
        if (i == "Group") {
          object(x)@meta.data[[i]] <- factor(
            object(x)@meta.data[[i]], rev(x$levels)
          )
        } else {
          stop('`i` is character, i != "Group", Please set the factor manually.')
        }
      }
    }
    if (!debug) {
      object(x) <- e(
        hdWGCNA::ModuleTraitCorrelation(
          object(x), traits = use.trait, group.by = group.by
        )
      )
    }
    x <- methodAdd(
      x, "使用 ModuleTraitCorrelation 计算各模块与性状 (trait, 即 {bind(use.trait)}) 之间的相关性。"
    )
    lst_cor <- e(hdWGCNA::GetModuleTraitCorrelation(object(x)))
    cor_module_traits <- sapply(as.character(celltypes), simplify = FALSE,
      function(type) {
        cp <- list(pvalue = lst_cor[["pval"]][[type]], cor = lst_cor[["cor"]][[type]])
        cp <- sapply(names(cp), simplify = FALSE,
          function(name) {
            obj <- cp[[name]]
            if (is.numeric(obj) && length(use.trait) == 1) {
              obj <- do.call(data.frame, as.list(obj))
              rownames(obj) <- use.trait
            }
            obj
          })
        safe_as_cor_tbl(cp, "cor", "pvalue")
      })
    x$cor_module_traits <- cor_module_traits
    snap_moTr <- .stat_ggcor_table_list(
      cor_module_traits, "Trait", "Module"
    )
    ts.cor_module_traits <- set_lab_legend(
      dplyr::bind_rows(cor_module_traits, .id = "Cell_Type"),
      glue::glue("{x@sig} All cells module traits correlation"),
      glue::glue("模块与性状关联性统计")
    )
    x <- tablesAdd(x, ts.cor_module_traits)
    ps.cor_module_traits <- lapply(names(cor_module_traits), 
      function(name) {
        data <- cor_module_traits[[name]]
        p <- .ggcor_add_general_style(ggcor::quickcor(data, type = "lower"))
        p + ggtitle(name)
        # wrap_scale_heatmap(ps, x$nm, length(use.trait), raw = FALSE)
      })
    layout <- wrap_layout(NULL, length(ps.cor_module_traits))
    ps.cor_module_traits <- patchwork::wrap_plots(
      ps.cor_module_traits, guides = "collect", ncol = layout$ncol
    )
    ps.cor_module_traits <- set_lab_legend(
      add(layout, ps.cor_module_traits),
      glue::glue("{x@sig} All cells module traits correlation heatmap"),
      glue::glue("模块与性状 (trait) 关联热图|||hdWGCNA 对 module 与 trait 的关联分析。热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。P 值以 * 标注 ({.md_p_significant})。")
    )
    x <- snapAdd(x, "如图{aref(ps.cor_module_traits)}，{snap_moTr}")
    x <- plotsAdd(x, ps.cor_module_traits)
    if (length(use.trait) > 1) {
      p.cor <- e(hdWGCNA::PlotModuleTraitCorrelation(
          object(x), label = "fdr", label_symbol = "stars",
          text_size = 2, text_digits = 2, text_color = "white",
          high_color = "yellow", mid_color = "black",
          low_color = "purple", plot_max = 0.2, combine = TRUE
          ))
      p.cor <- set_lab_legend(p.cor,
        glue::glue("{x@sig} module trait correlation heatmap"),
        glue::glue("")
      )
      x <- plotsAdd(x, p.cor)
    }
    return(x)
  })

setMethod("step7", signature = c(x = "job_hdwgcna"),
  function(x, cut.r = NULL, cut.mm = NULL,
    celltypes = x$celltypes, use.module = NULL,
    range_mm = c(100, 300), range_gs = c(100, 400),
    use.trait = "group", top = 1L, nlabel = 10, group.by = x$group.by)
  {
    step_message("Module Membership and Gene Significant (for module significant corrlated with trait).")
    if (length(use.trait) > 1) {
      stop('length(use.trait) > 1, not support now.')
    }
    if (top != 1L) {
      stop('top != 1L, not support now.')
    }
    if (use.trait == "group") {
      # this has been changed in previous step.
      use.trait <- "Group"
    }
    if (is.null(use.module)) {
      if (length(celltypes) > 1) {
        data <- dplyr::bind_rows(x$cor_module_traits, .id = ".celltypes")
        data <- dplyr::filter(
          data, .row.names == !!use.trait, p.value < .05
        )
        data <- dplyr::group_by(data, .celltypes)
        data <- dplyr::mutate(data, nSig = length(.row.names))
        data <- dplyr::arrange(data, dplyr::desc(nSig), dplyr::desc(abs(r)))
        use.module <- head(data$.col.names, n = top)
        x <- methodAdd(
          x, "根据 Module 与 Trait ({use.trait}) 的相关性分析，选择与最多细胞类型显著相关的模块 (若存在相同，则按最高关联的相关系数排序)，即 {use.module} 进一步筛选 Hub genes。"
        )
      } else {
        data <- dplyr::filter(
          x$cor_module_traits[[celltypes]], .row.names == !!use.trait
        )
        data <- dplyr::arrange(data, p.value)
        use.module <- head(data$.col.names, n = top)
        x <- methodAdd(
          x, "根据 Module 与 Trait ({use.trait}) 的相关性分析，选择最显著相关的模块，即 {use.module} 进一步筛选 Hub genes。"
        )
      }
    } else {
      x <- methodAdd(x, "选择 {use.module} 进一步筛选 Hub genes。")
    }
    message(
      glue::glue("Use module: {use.module}")
    )
    # do not filter now, for plot scatter.
    data_mm <- dplyr::filter(x$modules, module %in% !!use.module)
    data_mm <- dplyr::arrange(
      data_mm, dplyr::desc(!!rlang::sym(paste0("kME_", use.module)))
    )
    col_kme <- paste0("kME_", use.module)
    if (!any(colnames(data_mm) == col_kme)) {
      stop('!any(colnames(data_mm) == col_kme), no `ModuleConnectivity` running?')
    }
    if (is.null(cut.mm)) {
      cut.mm <- get_table_adapt_threshold(data_mm, col_kme, "high", range_mm, .1)
    }
    nMm <- nrow(dplyr::filter(data_mm, !!rlang::sym(col_kme) > cut.mm))
    message(glue::glue("Use `cut.mm` get number: {nMm}"))
    genes <- data_mm$gene_name
    message(
      glue::glue("The module contains gene number: {length(genes)}")
    )
    data_expr <- SeuratObject::GetAssayData(
      object(x), assay = SeuratObject::DefaultAssay(object(x)),
      layer = "data"
    )
    isTheCell <- object(x)@meta.data[[group.by]] %in% celltypes
    isTheGene <- rownames(data_expr) %in% genes
    gene_expr <- data_expr[isTheGene, isTheCell]
    message(
      glue::glue("Get expression data of dim: {bind(dim(gene_expr))}")
    )
    gene_expr <-  e(as.matrix(Seurat::ScaleData(gene_expr)))
    trait_data <- dplyr::select(object(x)@meta.data, use.trait)
    trait_data <- trait_data[isTheCell, , drop = FALSE]
    if (!identical(rownames(trait_data), colnames(gene_expr))) {
      stop('The `gene_expr` and `trait_data` not match cells.')
    }
    trait_data <- dplyr::mutate(
      trait_data, dplyr::across(dplyr::everything(),
        function(x) {
          if (is.factor(x)) {
            levels <- levels(x)
            message(glue::glue("Detected factor variable, the order: {bind(levels)}"))
            as.integer(x)
          } else x
        })
    )
    cor_gs <- safe_fortify_cor(t(gene_expr), trait_data)
    if (is.null(cut.r)) {
      .cor_gs <- dplyr::mutate(cor_gs, abs_r = abs(r))
      cut.r <- get_table_adapt_threshold(
        .cor_gs, "abs_r", "high", range_gs, .1
      )
    }
    cor_gs_filter <- dplyr::filter(cor_gs, abs(r) > !!cut.r)
    message(glue::glue("Use `cut.r` get number: {nrow(cor_gs_filter)}"))
    data_mm_gs <- merge(
      dplyr::select(cor_gs, gene = .row.names, gs_p.value = p.value, gs_r = r),
      dplyr::select(data_mm, gene = gene_name, mm = !!rlang::sym(col_kme)),
      by = "gene"
    )
    if (nrow(data_mm_gs) != nrow(data_mm) || nrow(data_mm_gs) != nrow(cor_gs)) {
      stop('After merge, the row number not match each other.')
    }
    data_mm_gs <- dplyr::mutate(
      data_mm_gs, type = ifelse(
        mm > cut.mm & abs(gs_r) > cut.r, "Hub_gene", "Normal"
      )
    )
    hubgenes <- dplyr::filter(data_mm_gs, type == "Hub_gene")$gene
    feature(x) <- as_feature(
      list(module_genes = data_mm$gene_name, hub_genes = hubgenes),
      "hdWGCNA Genes"
    )
    data_mm_gs <- set_lab_legend(
      data_mm_gs,
      glue::glue("{x@sig} GS MM data"),
      glue::glue("GS MM data")
    )
    x <- tablesAdd(x, data_mm_gs)
    p.scatter <- ggplot(data_mm_gs, aes(x = mm, y = abs(gs_r))) +
      geom_point(aes(color = type)) +
      scale_color_manual(values = c(Normal = "grey70", Hub_gene = use.module)) +
      theme_classic() +
      geom_vline(xintercept = cut.mm, linetype = 4) +
      geom_hline(yintercept = cut.r, linetype = 4) +
      labs(x = "Module Membership (kME)", y = "Gene Significant", color = "Type")
    p.scatter <- set_lab_legend(
      wrap(p.scatter, 4.7, 2.8),
      glue::glue("{x@sig} MM GS scatter plot"),
      glue::glue("Gene Significance 与 Module Membership|||图中 x 轴代表 Module Membership（由 `ModuleConnectivity` 计算的 kME），衡量基因表达谱与模块特征基因的相关性；y 轴代表基因显著性 (GS)，通常反映基因表达与外部性状（如疾病状态）的相关程度。")
    )
    x <- snapAdd(x, "根据 Module 与 Trait ({use.trait}) 的相关性分析，选择最显著的模块，即 {use.module} 之后，⟦mark$red('以 MM 与 GS 为条件共筛选到 {length(hubgenes)} 个基因')⟧{aref(p.scatter)}。")
    igraph <- e(hdWGCNA::HubGeneNetworkPlot(object(x),
      mods = use.module, sample_edges = FALSE, return_graph = TRUE,
      n_other = 0, n_hubs = nMm
    ))
    data_network <- tidygraph::as_tbl_graph(igraph)
    data_network <- dplyr::mutate(
      data_network, geneType = ifelse(
        name %in% !!hubgenes, "Hub_gene", "Normal"
      )
    )
    x$data_network <- data_network
    data_nodes <- tibble::as_tibble(data_network)
    coords <- get_coords.spiral(
      nrow(data_nodes), 1.3, minRad = .5, rank.by = data_nodes[[col_kme]]
    )
    data_graph <- create_layout(data_network, coords)
    data_label <- dplyr::filter(data_graph, name %in% hubgenes)
    if (!is.null(nlabel)) {
      data_label <- dplyr::slice_max(data_label, !!rlang::sym(col_kme), n = nlabel)
      snap_label <- glue::glue("标注了其中按 {col_kme} 排序的 Top {nlabel} 个基因。")
    }
    p.network <- ggraph(data_graph) +
      geom_edge_link(edge_color = "grey90") +
      geom_node_point(aes(color = !!rlang::sym(col_kme), shape = geneType, size = geneType)) +
      scale_shape_manual(values = c(Normal = 16, Hub_gene = 18)) +
      scale_size_manual(values = c(Normal = 3, Hub_gene = 10)) +
      ggrepel::geom_label_repel(
        data = data_label,
        aes(x = x, y = y, label = name), force_pull = -.1) +
      scale_color_gradient(low = "Orange", high = "Red4") +
      theme_void()
    p.network <- set_lab_legend(
      p.network,
      glue::glue("{x@sig} module {use.module} member networking"),
      glue::glue("螺旋型布局的 {use.module} 模块成员网络图|||仅展示 Module Membership (kME) 大于 {cut.mm} 的模块基因。用螺旋或同心圆状排列，从外向内 kME 值逐渐增大，即靠近中心的基因与该模块的特征基因（module eigengene）相关性更强、模块成员地位更高。不同的形状标注了是否为同时满足 Gene Significant (GS) 大于 {cut.r} 的 Hub Genes。{snap_label}")
    )
    x <- methodAdd(x, "⟦mark$blue('筛选标准为模块成员度（MM, 即 kME &gt; {cut.mm}）和基因显著性（GS, 即 |GS| &gt; {cut.r}）双重排序')⟧。")
    x <- snapAdd(x, "Hub Genes 在网络中的连接度如图所示{aref(p.network)}。")
    x <- plotsAdd(x, p.scatter, p.network)
    return(x)
  })

.stat_ggcor_table_list <- function(object, label.x, label.y, show = 3) {
  if (!is(object, "list") && is(object, "cor_tbl")) {
    lst <- list(object)
    leader <- ""
  } else if (is(object, "list")) {
    lst <- object
    leader <- glue::glue("{names(lst)} 的 ")
  }
  snaps <- vapply(lst, FUN.VALUE = character(1),
    function(data) {
      leader <- ""
      data <- dplyr::filter(data, .row.names != .col.names)
      data <- dedup_rows_by_unordered_cols(
        data, c(".row.names", ".col.names")
      )
      n_pre <- nrow(data)
      data <- dplyr::filter(data, p.value < .05)
      n_aft <- nrow(data)
      if (n_pre == n_aft && n_aft > 1) {
        if (label.x == label.y) {
          label <- label.x
        } else {
          label <- glue::glue("{label.x}与{label.y}")
        }
        leader <- glue::glue("所有{label}之间显著相关，")
      }
      data <- dplyr::arrange(data, dplyr::desc(abs(r)))
      snap <- .stat_correlation_table(
        data, ".row.names", ".col.names", "r", "p.value", label.x = label.x, label.y = label.y,
        maxShow = show
      )
      glue::glue("{leader}{bind(snap)}")
    })
  bind(glue::glue("{leader}{label.x} 与 {label.y} 关联分析结果，{snaps}。"), co = "\n\n")
}

safe_fortify_cor <- function(x, y, cor.test = TRUE, ...) {
  if (missing(y)) {
    y <- x
  }
  cor <- ggcor::correlate(x, y, cor.test = cor.test, ...)
  safe_as_cor_tbl(cor, "r", "p.value")
}

safe_as_cor_tbl <- function(lst, name_cor, name_pvalue) {
  if (!requireNamespace("ggcor", quietly = TRUE)) {
    stop('!requireNamespace("ggcor").')
  }
  if (any(!c(name_pvalue, name_cor) %in% names(lst))) {
    stop('any(!c(name_pvalue, name_cor) %in% names(lst)).')
  }
  mode <- character(1)
  lst <- lapply(lst,
    function(data) {
      if (!is.matrix(data) && !is.data.frame(data)) {
        rlang::abort(
          glue::glue("!is.matrix(data) && !is.data.frame(data), the class is: {bind(class(data))}")
        )
      }
      if (is.matrix(data)) {
        data <- as.data.frame(data)
      }
      dim <- dim(data)
      if (dim[1] == 1) {
        mode <<- "row"
        ex <- setNames(rep(.5, ncol(data)), colnames(data))
        data <- rbind(data, ex)
        rownames(data)[2] <- "...placeHolder"
      } else if (dim[2] == 1) {
        mode <<- "col"
        data[[ "...placeHolder" ]] <- .5
      }
      return(data)
    })
  data <- ggcor::cor_tbl(lst[[name_cor]], lst[[name_pvalue]])
  if (mode != character(1)) {
    type <- glue::glue(".{mode}.names")
    id <- glue::glue(".{mode}.id")
    data <- data[data[[type]] != "...placeHolder", ]
    data[[id]] <- 1L
    names <- attr(data, type)
    attr(data, type) <- names[ names != "...placeHolder" ]
  }
  return(data)
}

.auto_metacell_params <- function(
  seurat_obj, group.by, target_metacells = 10,
  k_min = 5, k_max = 30)
{
  # Get cell counts per group
  counts <- table(seurat_obj[[group.by]][,1])
  counts <- as.numeric(counts)
  names(counts) <- names(table(seurat_obj[[group.by]][,1]))
  
  # Smallest group size determines parameter upper bounds
  n_min <- min(counts)
  
  # ---- 1. Estimate k ----
  # Use sqrt(n) heuristic with constraints to avoid over-connection
  k <- floor(sqrt(n_min))
  k <- max(k, k_min)
  k <- min(k, floor(n_min / 2), k_max)
  
  # ---- 2. Estimate max_shared ----
  # Ensure enough metacells can be constructed
  # (n * max_shared) / k ≥ target_metacells
  max_shared <- ceiling((target_metacells * k) / n_min)
  max_shared <- max(max_shared, 2)
  max_shared <- min(max_shared, 10)
  
  # ---- 3. Set min_cells ----
  # Filter out extremely small groups while retaining most data
  min_cells <- floor(0.9 * n_min)
  
  # ---- 4. Estimate achievable metacell number ----
  est_metacells <- floor((n_min * max_shared) / k)
  
  message("---- Auto parameters ----")
  message("Min group size: ", n_min)
  message("k: ", k)
  message("max_shared: ", max_shared)
  message("min_cells: ", min_cells)
  message("Estimated metacells per smallest group: ", est_metacells)
  
  lines <- readLines(file.path(.expath, "description", "hdwgcna_metacell.md"))
  snap <- glue_ex(paste0(lines, collapse = "\n"))
  snap <- glue_ex(
    "⟦snap⟧\n\n\n\n最终，设置 $k = ⟦k⟧$，$min_{cells} = ⟦min_cells⟧$，$max_{shared} = ⟦max_shared⟧$。"
  )

  list(k = k, max_shared = max_shared,
    min_cells = min_cells, group_sizes = counts,
    est_metacells = est_metacells, snap = snap)
}

