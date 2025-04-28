# ==========================================================================
# workflow of WGCNA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_wgcna <- setClass("job_wgcna", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html"),
    cite = "[@WgcnaAnRPacLangfe2008]",
    method = "R package `WGCNA` used for gene co-expression analysis",
    tag = "wgcna",
    analysis = "WGCNA 分析"
    ))

setGeneric("asjob_wgcna", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_wgcna"))

setMethod("asjob_wgcna", signature = c(x = "job_seurat"),
  function(x, features = NULL, cells = NULL){
    step_message("Use SeuratObject:::subset.Seurat to subset the data.")
    hasIt <- features %in% rownames(object(x))
    message("Features found:")
    print(prop.table(table(hasIt)))
    sub <- e(suppressWarnings(SeuratObject:::subset.Seurat(object(x),
        features = features[ hasIt ], cells = cells
        )))
    log_counts <- as_tibble(sub[[ SeuratObject::DefaultAssay(sub) ]]@scale.data)
    metadata <- as_tibble(sub@meta.data)
    gene_annotation <- tibble::tibble(gene = rownames(sub))
    job_wgcna(metadata, log_counts, gene_annotation, "gene")
  })

job_wgcna <- function(metadata, log_counts,
  gene_annotation)
{
  elist <- new_elist(metadata, log_counts, gene_annotation)
  datExpr0 <- as_wgcData(elist)
  .job_wgcna(object = elist, params = list(datExpr0 = datExpr0))
}

setMethod("step0", signature = c(x = "job_wgcna"),
  function(x){
    step_message("Prepare your data with function `job_wgcna`. ",
      "Note that the ", crayon::red("first column"),
      " in each data (metadata, counts, genes annotation)",
      " were used as ID column of ",
      "corresponding content. \n",
      crayon::red("metadata:"), " with 'sample' and 'group' (optional).\n",
      crayon::red("counts:"), " with id column and then expression columns.\n",
      "genes: ", crayon::red("Traits data"), " could be placed herein; ",
      "a step would serve all numeric data in `genes` as trait data then ",
      "perform statistic."
    )
  })

setMethod("step1", signature = c(x = "job_wgcna"),
  function(x, mutate_name = TRUE){
    step_message("Cluster sample tree.",
    "This do:",
    "generate `x@params$raw_sample_tree`; `x@plots[[ 1 ]]`"
    )
    if (mutate_name) {
      dat <- params(x)$datExpr0
      rownames(dat) <- paste0(object(x)$targets$group, "_", seq_along(object(x)$targets$group))
      raw_sample_tree <- draw_sampletree(dat)
    } else {
      raw_sample_tree <- draw_sampletree(params(x)$datExpr0)
    }
    raw_sample_tree.p <- wrap(
      recordPlot(), 
      min(nrow(params(x)$datExpr0) * .4, 20), 
      min(nrow(params(x)$datExpr0), 10)
    )
    raw_sample_tree.p <- .set_lab(raw_sample_tree.p, sig(x), "sample clustering")
    raw_sample_tree.p <- setLegend(raw_sample_tree.p, "为样本聚类树")
    x@params$raw_sample_tree <- raw_sample_tree
    x@plots[[ 1 ]] <- list(raw_sample_tree = raw_sample_tree.p)
    x <- methodAdd(x, "以 R 包 `WGCNA` ({packageVersion('WGCNA')}) {cite_show('WgcnaAnRPacLangfe2008')} 对数据作共表达分析。分析方法参考 <{x@info}>。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_wgcna"),
  function(x, height = NULL, size = NULL){
    step_message("Cut sample tree with `height` and `size`. ",
      "This do: ",
      "clip `x@object`; generate `x@params$datExpr`; ",
      "generate `x@params$allTraits`. "
    )
    if (!is.null(size) || !is.null(height)) {
      iskeep <- cut_tree(params(x)$raw_sample_tree, height, size)
      message(glue::glue("Keep: {try_snap(iskeep)}\nDrop:\n{showStrings(which(!iskeep))}"))
      datExpr <- exclude(params(x)$datExpr0, iskeep)
      x@params$datExpr <- datExpr
      object(x) <- clip_data(object(x), datExpr)
      x <- snapAdd(
        x, "以 `WGCNA::cutreeStatic` (cutHeight = {height}, minSize = {size}) 剪切聚类树，滤掉样本 {showStrings(rownames(x$datExpr0)[!iskeep])}。"
      )
    } else {
      x$datExpr <- x$datExpr0
    }
    x@params$allTraits <- as_wgcTrait(object(x))
    return(x)
  })

setMethod("step3", signature = c(x = "job_wgcna"),
  function(x, cores = 4, powers = 1:50, ...)
  {
    step_message("Analysis of network topology for soft-thresholding powers. ",
      "This do: ",
      "Generate x@params$sft; plots in `x@plots[[ 3 ]]`. "
    )
    if (is.remote(x)) {
      object <- object(x)
      object(x) <- NULL
      x <- run_job_remote(x, wait = 3L, ...,
        {
          x <- step3(x, powers = seq_len("{max(powers)}"), cores = "{cores}")
        }
      )
      object(x) <- object
    } else {
      e(WGCNA::enableWGCNAThreads(cores))
      sft <- cal_sft(params(x)$datExpr, powers = powers)
      x@params$sft <- sft
      p.sft <- wrap(plot_sft(sft), 10, 5)
      p.sft <- .set_lab(p.sft, sig(x), "soft thresholding powers")
      p.sft <- setLegend(p.sft, "WGCNA 软阈值筛选曲线。")
      x@plots[[ 3 ]] <- list(sft = p.sft)
      x <- methodAdd(x, "以 `WGCNA::pickSoftThreshold` 预测最佳 soft thresholding powers。")
    }
    return(x)
  })

setMethod("step4", signature = c(x = "job_wgcna"),
  function(x, cores = 4, power = x@params$sft$powerEstimate, 
    inherit = TRUE, ...)
  {
    step_message("One-step network construction and module detection.
      Extra parameters would passed to `cal_module`.
      This do: Generate `x@params$MEs`; plots (net) in `x@plots[[ 4 ]]`.
      By default, red{{x@params$sft$powerEstimate}} is used
      as `power` for WGCNA calculation.
      "
    )
    if (is.remote(x)) {
      message(glue::glue("Note: `...` can not passed to remote."))
      object <- object(x)
      object(x) <- NULL
      x <- run_job_remote(
        x, wait = 3, inherit_last_result = inherit, ...,
        {
          x <- step4(x, power = "{power}", cores = "{cores}")
        }
      )
      object(x) <- object
    } else {
      e(WGCNA::enableWGCNAThreads(cores))
      net <- cal_module(params(x)$datExpr, power, ...)
      if (!is(net, "wgcNet")) {
        net <- .wgcNet(net)
      }
      x@params$MEs <- get_eigens(net)
      ME_genes <- net$colors
      ME_genes <- split(names(ME_genes), unname(ME_genes))
      names(ME_genes) <- paste0("ME", names(ME_genes))
      x$ME_genes <- ME_genes
      x <- snapAdd(x, "以 power {power} (soft thresholding powers) 创建基因共表达模块 (各模块基因数：{try_snap(ME_genes)})。")
      net <- .set_lab(
        wrap(net, 7, 6), sig(x), "co-expression module"
      )
      net <- setLegend(net, "为 WGCNA 创建的网络的基因共表达模块。")
      x@plots[[ 4 ]] <- list(net = net)
      x <- methodAdd(x, "选择 power 为 {power}, 以 `WGCNA::blockwiseModules` 创建共表达网络，检测基因模块。")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_wgcna"),
  function(x, traits = NULL, group_levels = NULL, cut.p = .05, cut.cor = .3)
  {
    step_message("Correlation test for modules with trait data. ",
      "This do:",
      "Generate plots in `x@plots[[ 5 ]]`; ",
      "tables in `x@tables[[ 5 ]]`"
    )
    if (!is.null(x$traits)) {
      message("Use `x$traits` for correlation.")
      traits <- x$traits
    }
    if (is.null(traits) && !is.null(group_levels) && !is.null(object(x)$targets[[ "group" ]])) {
      message(glue::glue("Use 'group' in `object(x)$targets`: {bind(group_levels)}."))
      traits <- object(x)$targets
      traits$group <- as.integer(factor(traits$group, levels = group_levels))
      x <- snapAdd(x, "将 'group' 设置为数值变量 ({bind(group_levels)} 依次为 {bind(seq_along(group_levels))}) 与基因共表达模块关联分析。")
    }
    if (!is.null(traits)) {
      .check_columns(traits, c("sample"))
      message("Match rownames in expression data.")
      traits <- traits[match(rownames(x@params$datExpr), traits$sample), ]
      rownames <- traits$sample
      message("The numeric columns will calculate correlation with expression data.")
      traits <- dplyr::select_if(traits, is.numeric)
      traits <- data.frame(traits)
      rownames(traits) <- rownames
      x$allTraits <- .wgcTrait(traits)
    }
    if (is.null(params(x)$allTraits)) {
      stop("is.null(params(x)$allTraits) == TRUE")
    }
    if (ncol(params(x)$allTraits) == 0) {
      stop("ncol(params(x)$allTraits) == 0, no data in `allTraits`.")
    }
    if (ncol(params(x)$allTraits) == 1) {
      traitName <- colnames(params(x)$allTraits)
      cor <- e(WGCNA::cor(params(x)$MEs, params(x)$allTraits, use = "p"))
      pvalue <- e(WGCNA::corPvalueStudent(cor, nrow(params(x)$MEs)))
      if (!identical(rownames(cor), rownames(pvalue))) {
        stop('!identical(rownames(cor), rownames(pvalue))')
      }
      x$corp_group <- dplyr::bind_cols(cor, pvalue)
      colnames(x$corp_group) <- c("cor", "pvalue")
      x$corp_group <- dplyr::mutate(x$corp_group, MEs = rownames(!!cor), .before = 1)
      x$corp_group <- dplyr::arrange(x$corp_group, dplyr::desc(abs(cor)))
      x$corp_group <- .set_lab(
        x$corp_group, sig(x), "correlation of module with", traitName 
      )
      x$corp_group <- setLegend(x$corp_group, "为共表达模块与 {traitName} 的关联性。")
      data <- dplyr::mutate(x$corp_group, group = !!traitName)
      fun_palette <- fun_color(
        values = data$cor, category = "div", rev = TRUE
      )
      p.corhp <- e(
        tidyHeatmap::heatmap(data, MEs, group, cor, palette_value = fun_palette)
      )
      p.corhp <- tidyHeatmap::layer_text(
        p.corhp, .value = signif(pvalue, 4)
      )
      p.corhp <- set_lab_legend(
        wrap(p.corhp, 5),
        glue::glue("{x@sig} correlation heatmap"),
        glue::glue("为模块与疾病表型 ({traitName}) 关联分析热图。")
      )
      x <- plotsAdd(x, p.corhp)
      sigModules <- dplyr::filter(
        x$corp_group, abs(cor) > cut.cor, pvalue < cut.p
        )$MEs
      x <- snapAdd(
        x, "筛选显著关联的共表达模块的基因 (pvalue &lt; {cut.p}, cor &gt; {cut.cor})。"
      )
      x$.feature <- as_feature(
        x$ME_genes[ names(x$ME_genes) %in% sigModules ], x,
        analysis = glue::glue("WGCNA 与 {traitName} 显著关联的共表达模块的基因")
      )
    } else {
      hps_corp <- new_heatdata(params(x)$MEs, params(x)$allTraits)
      hps_corp <- callheatmap(hps_corp)
      x@plots[[ 5 ]] <- list(hps_corp = hps_corp)
      x@tables[[ 5 ]] <- list(corp = hps_corp@data_long)
    }
    return(x)
  })

setMethod("step6", signature = c(x = "job_wgcna"),
  function(x, use.trait = NULL, use = c("adj.pvalue", "pvalue")){
    step_message("Calculate gene significance (GS) and module membership (MM).",
      "This do:",
      "Generate `x@params$mm`, `x@params$gs`; ",
      "tables (filter by pvalue < 0.05) `x@tables[[ 6 ]]`"
    )
    use <- match.arg(use)
    mm <- cal_corp(params(x)$datExpr, params(x)$MEs, "gene", "module")
    mm.s <- mutate(as_tibble(mm), adj.pvalue = p.adjust(pvalue, "BH"))
    mm.s <- dplyr::filter(mm.s, !!rlang::sym(use) < .05)
    gs <- cal_corp(params(x)$datExpr, params(x)$allTraits, "gene", "trait")
    gs.s <- dplyr::mutate(tibble::as_tibble(gs), adj.pvalue = p.adjust(pvalue, "BH"))
    gs.s <- dplyr::filter(gs.s, !!rlang::sym(use) < .05)
    if (!is.null(use.trait)) {
      gs.s <- dplyr::filter(gs.s, trait %in% dplyr::all_of(use.trait))
    }
    p.mm_gs <- new_upset(gs = gs.s$gene, mm = mm.s$gene)
    show(p.mm_gs)
    p.mm_gs <- wrap(recordPlot(), 3, 3)
    dev.off()
    x@params$mm <- mm
    x@params$gs <- gs
    x@params$ins.mm_gs <- intersect(gs.s$gene, mm.s$gene)
    x@tables[[ 6 ]] <- list(mm = mm.s, gs = gs.s)
    x@plots[[ 6 ]] <- namel(p.mm_gs)
    return(x)
  })

cut_tree <- function(tree, height, size) {
  clust <- e(WGCNA::cutreeStatic(tree, height, size))
  clust > 0
}

cal_sft <- function(data, powers = c(c(1:10), seq(12, 20, by = 2))) 
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == FALSE")
  }
  sft <- e(WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5))
  sft
}

plot_sft <- function(sft) 
{
  p1 <- ggplot(sft$fitIndices, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
    geom_line(color = "darkred", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2") +
    theme_classic()
  p2 <- ggplot(sft$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_line(color = "darkgreen", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Mean Connectivity") +
    theme_classic()
  require(patchwork)
  p1 + p2
}

cal_module <- function(data, power, save_tom = "tom", ...)
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == FALSE")
  }
  require(WGCNA)
  net <- e(WGCNA::blockwiseModules(
      data, power = power,
      TOMType = "unsigned", reassignThreshold = 0,
      numericLabels = TRUE, pamRespectsDendro = FALSE, loadTOM = TRUE,
      saveTOMs = TRUE, saveTOMFileBase = save_tom, verbose = 3, ...
      ))
  .wgcNet(net)
}

setMethod("set_remote", signature = c(x = "job_wgcna"),
  function(x, wd = glue::glue("~/wgcna_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

