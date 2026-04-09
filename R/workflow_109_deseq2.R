# ==========================================================================
# workflow of deseq2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_deseq2 <- setClass("job_deseq2", 
  contains = c("job"),
  prototype = prototype(
    pg = "deseq2",
    info = c("https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html"),
    cite = "[@Moderated_estim_Love_2014]",
    method = "",
    tag = "deseq2",
    analysis = "DESeq2 差异分析"
    ))

setGeneric("asjob_deseq2",
  function(x, ...) standardGeneric("asjob_deseq2"))

setMethod("asjob_deseq2", signature = c(x = "job_geo"),
  function(x, metadata, use.col = NULL, use_as_id = TRUE, 
    use = 1L, ...)
  {
    if (!x$rna) {
      stop('!x$rna, not RNA-seq data!')
    }
    counts <- as_tibble(data.frame(x@params$about[[ use ]]@assays@data$counts, check.names = FALSE))
    genes <- as_tibble(data.frame(x@params$about[[ use ]]@elementMetadata, check.names = FALSE))
    genes <- dplyr::mutate(genes, rownames = as.character(GeneID), .before = 1)
    if (missing(use.col)) {
      use.col <- .guess_symbol(colnames(genes))
      message(glue::glue("Guess `use.col` (gene symbol) is: {use.col}."))
    } else {
      if (!is.character(genes[[ use.col ]])) {
        stop('`use.col` should character (symbol).')
      }
    }
    genes <- dplyr::relocate(genes, rownames, symbol = !!rlang::sym(use.col))
    use.col <- "symbol"
    if (!is.null(x$ncbiNotGot)) {
      message(glue::glue("Filter out NCBI generated data (Missing samples): {bind(x$ncbiNotGot)}"))
      metadata <- dplyr::filter(metadata, !rownames %in% x$ncbiNotGot)
    }
    if (use_as_id) {
      message(crayon::red(glue::glue("Use `{use.col}` as ID columns.")))
      counts[[1]] <- genes[[ use.col ]]
      genes <- dplyr::relocate(genes, !!rlang::sym(use.col))
      keep <- !is.na(genes[[ use.col ]]) & genes[[ use.col ]] != ""
      genes <- genes[keep, ]
      counts <- counts[keep, ]
    }
    x <- job_deseq2(counts, metadata, genes, x$project, ...)
    return(x)
  })

job_deseq2 <- function(counts, metadata, genes, project = NULL, 
  formula = ~ 0 + group, recode = NULL)
{
  # Sort the average gene expression levels and retain the gene symbol with the highest average expression level
  if (!is.null(recode)) {
    if (!all(names(recode) %in% counts[[1]])) {
      stop('!all(names(recode) %in% counts[[1]]), no that genes involved?')
    }
    counts[[1]] <- dplyr::recode(counts[[1]], !!!recode)
    genes[[1]] <- counts[[1]]
  }
  object <- prepare_expr_data(metadata, counts, genes)
  genes <- object$genes
  object <- DESeq2::DESeqDataSetFromMatrix(
    object$counts, data.frame(object$metadata, row.names = object$metadata[[1]]),
    mx(formula, data = metadata)
  )
  x <- .job_deseq2(object = object)
  x$genes <- genes
  x$idcol <- colnames(genes)[1]
  x$formula <- formula
  if (is.null(project)) {
    message("Note: no project name set in 'job'.")
  } else {
    x$project <- project
  }
  x <- methodAdd(x, "差异表达分析 (Differential Expression Analysis) 使用 DESeq2 包 ⟦pkgInfo('DESeq2')⟧ 进行基因差异表达分析。构建 'DESeqDataSet' 对象，将各样本的原始 Count 矩阵与分组信息导入。{prepare_expr_data()}。")
  return(x)
}

# setValidity("job_deseq2",
#   function(object){
#     if (is.null(object$project)) {
#       message('is.null(object$project), Please set "project".')
#       return(FALSE)
#     } else {
#       TRUE
#     }
#   })

setMethod("step0", signature = c(x = "job_deseq2"),
  function(x){
    step_message("Prepare your data with function `job_deseq2`.")
  })

setMethod("step1", signature = c(x = "job_deseq2"),
  function(x){
    step_message("Quality control (QC).")
    object(x) <- e(DESeq2::DESeq(object(x)))
    x$vst <- e(DESeq2::vst(object(x), blind = FALSE))
    p.pca <- DESeq2::plotPCA(x$vst, intgroup = "group")
    p.pca <- set_lab_legend(
      wrap(p.pca, 7, 6),
      glue::glue("{x@sig} PCA of VST data"),
      glue::glue("样本 PCA 聚类图")
    )
    cook_data <- log10(object(x)@assays@data[[ "cooks" ]])
    colnames(cook_data) <- colnames(object(x))
    p.boxplot <- as_grob(expression(boxplot(cook_data)))
    p.boxplot <- set_lab_legend(
      p.boxplot,
      glue::glue("{x@sig} boxplot of cook distance"),
      glue::glue("样本 cook 距离箱线图。")
    )
    x <- plotsAdd(x, p.pca, p.boxplot)
    x <- methodAdd(x, "使用 DESeq 函数进行标准化和差异分析。DESeq 函数会针对每个基因和每个样本计算一个称为 Cook 距离的异常值诊断测试。绘制 Cook 距离的箱线图，以查看是否存在某个样本始终高于其他样本。为进一步质控数据，另一方面，对数据集以 vst 函数处理，随后绘制 PCA 图检查样本是否存在批次效应或离群样本。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_deseq2"),
  function(x, ..., cut.p = .05, cut.fc = 1, use = c("padj", "pvalue"), 
    order.by = "log2FoldChange")
  {
    step_message("stat.")
    use <- match.arg(use)
    if (!...length()) {
      stop('!...length().')
    }
    contrast <- limma::makeContrasts(..., levels = object(x)@design)
    cli::cli_alert_info("DESeq2::results")
    t.results <- apply(contrast, 2, simplify = FALSE,
      FUN = function(numVec) {
        data <- as_tibble(
          data.frame(DESeq2::results(object(x), contrast = numVec)), idcol = x$idcol
        )
        dplyr::arrange(data, dplyr::desc(abs(!!rlang::sym(order.by))))
      })
    names(t.results) <- .fmt_contr_names(colnames(contrast))
    t.sigResults <- lapply(t.results,
      function(data) {
        dplyr::filter(data, !!rlang::sym(use) < cut.p, abs(log2FoldChange) > cut.fc)
      })
    t.sigResults <- set_lab_legend(
      t.sigResults,
      glue::glue("{x@sig} {names(t.sigResults)} Contrast significant result table"),
      glue::glue("差异分析 {names(t.results)} 的显著基因的结果表格")
    )
    t.results <- set_lab_legend(
      t.results,
      glue::glue("{x@sig} {names(t.results)} Contrast all result table"),
      glue::glue("差异分析 {names(t.results)} 的全部基因的结果表格")
    )
    res <- .collate_snap_and_features_by_logfc(t.sigResults, "log2FoldChange", get = x$idcol)
    x <- methodAdd(
      x, "对 DESeq 函数处理过后的数据获取差异基因。⟦mark$blue('筛选差异表达基因（DEGs）的标准为 {use} &lt; {cut.p} 且 |log2(FoldChange)| &gt; {cut.fc}')⟧。"
    )
    x <- snapAdd(x, "显著基因统计：\n⟦mark$red('{res$snap}')⟧\n")
    feature(x) <- res$sets
    x <- tablesAdd(x, t.results, t.sigResults)
    return(x)
  })

setMethod("step3", signature = c(x = "job_deseq2"),
  function(x, top = 10, group = "group", use = x$.args$step2$use[1])
  {
    step_message("Draw plots.")
    use.fc <- "log2FoldChange"
    fc <- x$.args$step2$cut.fc
    cut.p <- x$.args$step2$cut.p
    data_expr <- SummarizedExperiment::assay(x$vst)
    metadata <- SummarizedExperiment::colData(x$vst)
    top_tables <- x@tables$step2$t.results
    idcol <- x$idcol
    x <- .plot_DEG_and_add_snap(
      x, data_expr, metadata, top_tables,
      use = use, use.fc = use.fc, fc = fc, cut.p = cut.p, top = top, 
      group = group, features = feature(x), idcol = idcol
    )
    return(x)
  })

.plot_DEG_and_add_snap <- function(x, data_expr, metadata, 
  top_tables, use, use.fc, fc, cut.p, top, group, features, idcol)
{
  topDegs <- sapply(
    names(top_tables), simplify = FALSE,
    function(name) {
      tops <- unlist(lapply(features[[ name ]], head, n = top))
      p.volcano <- plot_volcano(
        top_tables[[name]], idcol, use = use,
        fc = fc, cut.p = cut.p, use.fc = use.fc, HLs = tops
      )
      p.volcano <- set_lab_legend(
        p.volcano,
        glue::glue("{x@sig} {name} volcano of DEGs"),
        glue::glue("{name} 差异表达基因的火山图|||横坐标为Log2(Fold Change)，纵坐标为-Log10({detail(use)})，每个点代表一个基因；横向参考线代表-Log10({cut.p})=1.3，纵向参考线代表log2FC = ± {fc}；以参考线为划分，右上角的基因为疾病组相较对照组显著上调的差异表达基因，左上角的基因为显著下调的差异表达基因，其余基因为不具有显著统计学意义的基因（用灰色表示）；图中标注基因为按 |log2FC| 从高到低上下调 Top 10 差异表达基因。")
      )
      data <- data_expr
      data <- data[ rownames(data) %in% tops, , drop = FALSE ]
      data <- t(scale(t(data)))
      data <- pmin(pmax(data, -2), 2)
      orderSample <- order(metadata[[ group ]])
      metadata <- metadata[orderSample, ]
      data <- data[, orderSample]
      if (ncol(data) > 20) {
        isShowSample <- FALSE
      } else {
        isShowSample <- TRUE
      }
      p.hp <- pheatmap::pheatmap(
        data, annotation_col = data.frame(metadata)[, group, drop = FALSE],
        show_rownames = TRUE, cluster_rows = TRUE, 
        cluster_cols = FALSE, scale = "none",
        show_colnames = isShowSample
      )
      p.hp <- set_lab_legend(
        wrap_scale(p.hp, ncol(data), nrow(data), size = .2),
        glue::glue("{x@sig} {name} heatmap of top {top} DEGs"),
        glue::glue("Top {top} 差异表达基因热图|||热图中纵坐标代表基因，横坐标代表每一个样本；每一个‘格子’代表该行对应的基因在该列对应的样本中的表达水平。")
      )
      return(namel(p.hp, p.volcano, tops))
    }
  )
  psVol <- lapply(topDegs, function(x) x$p.volcano)
  psHp <- lapply(topDegs, function(x) x$p.hp)
  x <- snapAdd(x, "为了直观显示差异基因的分布，以 R 包 `ggplot2` ⟦pkgInfo('ggplot2')⟧绘制火山图对差异基因进行可视化{aref(psVol)}。根据差异倍数 ∣log2FC| 从高到低排序的上下调 Top 10 差异表达基因，以 R 包 `ComplexHeatmap` ⟦pkgInfo('ComplexHeatmap')⟧ 绘制DEGs的表达热图{aref(psHp)}。")
  x <- plotsAdd(x, topDegs)
  return(x)
}

setMethod("filter", signature = c(x = "job_deseq2"),
  function(x, ..., trace = FALSE){
    # filter sample
    data <- data.frame(object(x)@colData)
    data$...id <- seq_len(nrow(data))
    if (trace) {
      data <- trace_filter(data, ...)
    } else {
      data <- dplyr::filter(data, ...)
    }
    object(x) <- object(x)[, data$...id]
    object(x)@design <- NULL
    return(x)
  })

setMethod("regroup", signature = c(x = "job_deseq2", ref = "job_deseq2"),
  function(x, ref, sample = "sample", group = "group", formula = x$formula, sort = TRUE)
  {
    meta.x <- object(x)@colData
    meta.ref <- object(ref)@colData
    meta.x[[group]] <- meta.ref[[group]][ match(meta.x[[sample]], meta.ref[[sample]]) ]
    counts <- object(x)@assays@data$counts
    if (sort) {
      order <- order(meta.x[[group]])
      meta.x <- meta.x[order, ]
      counts <- counts[, order]
    }
    object(x) <- e(DESeq2::DESeqDataSetFromMatrix(
        counts, meta.x, mx(formula, data = meta.x)
        ))
    x$formula <- formula
    return(x)
  })

setMethod("regroup", signature = c(x = "job_deseq2", ref = "character"),
  function(x, ref, mode = "median", names = c("High", "Low"), col = "group", .add = FALSE)
  {
    data <- DESeq2::counts(object(x), normalized = TRUE)
    data <- log2(data + 1)
    if (length(ref) > 1) {
      stop('length(ref) > 1.')
    }
    if (!any(rownames(data) == ref)) {
      stop('!any(rownames(data) == ref).')
    }
    if (mode == "median") {
      cutoff <- median(data[ref, ])
      object(x)@colData[[ col ]] <- ifelse(data[ref, ] > cutoff, names[1], names[2])
    }
    x <- snapAdd(
      x, "根据 {ref} 的表达量，按 '{mode}' 值将样本重新分为 {bind(names)} 组。",
      step = "regroup", add = .add
    )
    return(x)
  })

.guess_compare_deseq2 <- function(x, which) {
  .get_group_from_contrast_character(names(x@tables$step2$t.results)[which])
}

setMethod("focus", signature = c(x = "job_deseq2"),
  function(x, ref, ref.use = "guess", which = NULL, run_roc = TRUE,
    which.roc = 1L, level.roc = .guess_compare_deseq2(x, which.roc),
    .name = "m", use = c("adj.P.Val", "P.Value"), 
    sig = FALSE, clear = "auto", test = "wilcox.test", ...)
  {
    # if which set to NULL, wilcox.test will be used.
    # else, the test results in 'results' table will be extracted.
    use <- match.arg(use)
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    data <- SummarizedExperiment::assay(x$vst)
    fakeLmJob <- .job_limma(sig = x@sig, step = 2L)
    allData <- lapply(
      x@tables$step2$t.results, dplyr::rename,
      adj.P.Val = padj, P.Value = pvalue, logFC = log2FoldChange
    )
    fakeLmJob@tables$step2$tops <- allData
    if (!is.null(which)) {
      data.which <- allData[[ which ]]
    } else {
      data.which <- NULL
    }
    fakeLmJob$normed_data <- elist <- new_from_package("EList", "limma",
      list(E = data, targets = data.frame(object(x)@colData), genes = x$genes)
    )
    if (identical(ref.use, "guess")) {
      ref.use <- .guess_symbol(fakeLmJob)
    }
    fakeLmJob <- suppressMessages(focus(
      fakeLmJob, ref = ref, ref.use = ref.use,
      .name = .name, which = which, data.which = data.which, sig = sig, 
      , use = use, test = test, run_roc = run_roc, ...
    ))
    where <- paste0("focusedDegs_", .name)
    if (identical(clear, "auto")) {
      clear <- !is.null(x[[ where ]])
    }
    resStat <- x[[ where ]] <- fakeLmJob[[ where ]]
    snap <- fakeLmJob@snap[[ paste0("step", .name) ]]
    x <- snapAdd(x, snap, step = .name, add = FALSE)
    return(x)
  })

setMethod("asjob_iobr", signature = c(x = "job_deseq2"),
  function(x, idType = "Symbol", source = c("local", "biomart"),
    levels = NULL, guess_which_level = 1L, ...)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    if (is.null(levels)) {
      levels <- .get_group_from_contrast_character(names(x@tables$step2$t.results)[guess_which_level])
    }
    mtx <- SummarizedExperiment::assay(object(des.iua))
    message("Transform the count to TPM.")
    message(glue::glue("The data dim: {bind(dim(mtx))}"))
    require(IOBR)
    dir.create("tmp", FALSE)
    source <- match.arg(source)
    args <- list(countMat = mtx, idType = idType, source = source)
    mtx <- expect_local_data(
      "tmp", "iobr_tpm", IOBR::count2tpm, args
    )
    # mtx <- e(IOBR::count2tpm(mtx, idType = idType, source = source, ...))
    message(glue::glue("The data dim: {bind(dim(mtx))}"))
    metadata <- data.frame(x$vst@colData)
    vst <- SummarizedExperiment::assay(x$vst)
    x <- job_iobr(mtx, metadata = metadata)
    x$levels <- levels
    x$vst <- vst
    return(x)
  })

.evaluate_by_ROC <- function(data, metadata, gene, 
  control = "control", sig)
{
  if (length(gene) > 1) {
    stop('length(gene) > 1.')
  }
  expr <- data[rownames(data) == gene, ]
  data <- data.frame(
    expression = expr, group = metadata$group[match(colnames(data), metadata$sample)]
  )
  data$group <- ifelse(data$group == control, 0L, 1L)
  roc <- e(pROC::roc(data$group, data$expression))
  auc <- pROC::auc(roc)
  ci <- pROC::ci.auc(roc)
  best_coords <- pROC::coords(roc, "best", best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv", "accuracy"))
  p.roc <- wrap(as_grob(expression(pROC::plot.roc(roc,
    main = "ROC Curve",
    print.auc = TRUE, print.auc.pattern = "AUC: %.3f",
    auc.polygon = TRUE, auc.polygon.col = "#90EE9020",
    print.auc.x = ifelse(roc$percent, 30, .3),
    print.auc.y = ifelse(roc$percent, 10, .1),
    grid = c(0.1, 0.1), legacy.axes = TRUE
  ))))
  p.roc <- set_lab_legend(
    wrap(p.roc, 6, 4.5),
    glue::glue("{sig} ROC plot of {gene}"),
    glue::glue("为 {gene} ROC 曲线|||{detail('note_roc')}")
  )
  snap <- glue::glue(
    paste0(
      "{gene} 的 ROC 分析结果{aref(p.roc)}，AUC = {round(auc, 3)}，AUC 的 95% 置信区间 = [{round(ci[1], 3)}, {round(ci[3], 3)}]；",
      "最佳诊断临界值（基于约登指数），临界值 (threshold) = {round(best_coords$threshold, 2)}，敏感度 (Sensitivity) = {round(best_coords$sensitivity, 3)}，特异度 (Specificity) = {round(best_coords$specificity, 3)}。"
    )
  )
  t.best_coords <- best_coords
  t.best_coords <- .set_lab(
    t.best_coords, sig, "best coords of", gene
  )
  t.best_coords <- setLegend(t.best_coords, "{gene}的最佳临界值诊断结果表格。")
  namel(p.roc, t.best_coords, auc, ci, roc, snap)
}

.collate_snap_and_features_by_logfc <- function(lst, col = "log2FoldChange", get = "rownames", prefix = "", cut = 0)
{
  if (is.null(names(lst))) {
    stop('is.null(names(lst)).')
  }
  fun_stat <- function(data, fun = `>`) {
    which <- fun(data[[ col ]], cut)
    sum <- sum(which)
    alls <- data[[ get ]][ which ]
    list(alls = alls, sum = sum)
  }
  res <- lapply(lst,
    function(data) {
      ups <- fun_stat(data, `>`)
      downs <- fun_stat(data, `<`)
      snap <- glue::glue("显著上调基因数量为{ups$sum}，显著下调基因数量为{downs$sum}")
      sets <- list(up = ups$alls, down = downs$alls)
      list(snap = snap, sets = sets)
    })
  res.snap <- lapply(res, function(x) x$snap)
  sets <- lapply(res, function(x) x$sets)
  snap <- paste0(
    paste0(
      prefix, "对比组 ", names(res.snap), 
      "，", res.snap, "。"
    ), collapse = "\n"
  )
  return(list(snap = snap, sets = sets))
}

setMethod("set_remote", signature = c(x = "job_deseq2"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
