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

job_deseq2 <- function(counts, metadata, genes, project = NULL, formula = ~ 0 + group)
{
  # Sort the average gene expression levels and retain the gene symbol with the highest average expression level
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
  x <- methodAdd(x, "差异表达分析 (Differential Expression Analysis) 使用 DESeq2 包 ({packageVersion('DESeq2')}) 进行基因差异表达分析。构建 'DESeqDataSet' 对象，将各样本的原始 Count 矩阵与分组信息导入。对于基因名重复，以基因表达量平均值排序后保留平均值最大的基因。")
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
    p.pca <- .set_lab(
      wrap(p.pca, 7, 6), sig(x), "PCA of VST data"
    )
    p.pca <- setLegend(p.pca, "样本 PCA 聚类图。")
    cook_data <- log10(object(x)@assays@data[[ "cooks" ]])
    colnames(cook_data) <- colnames(object(x))
    p.boxplot <- as_grob(expression(boxplot(cook_data)))
    p.boxplot <- .set_lab(p.boxplot, sig(x), "boxplot of cook distance")
    p.boxplot <- setLegend(p.boxplot, "样本 cook 距离箱线图。")
    x <- plotsAdd(x, p.pca, p.boxplot)
    x <- methodAdd(x, "使用 'DESeq()' 函数进行标准化和差异分析。DESeq 函数会针对每个基因和每个样本计算一个称为 Cook 距离的异常值诊断测试。绘制 Cook 距离的箱线图，以查看是否存在某个样本始终高于其他样本。为进一步质控数据，另一方面，对数据集以 vst 函数处理，随后绘制 PCA 图检查样本是否存在批次效应或离群样本。")
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
    names(t.results) <- sub(" - ", " vs ", colnames(contrast))
    t.sigResults <- lapply(t.results, 
      function(data) {
        dplyr::filter(data, !!rlang::sym(use) < cut.p, abs(log2FoldChange) > cut.fc)
      })
    t.sigResults <- .set_lab(t.sigResults, sig(x), names(t.sigResults), "Contrast significant result table")
    t.sigResults <- setLegend(t.sigResults, glue::glue("差异分析 {names(t.results)} 的显著基因的结果表格。"))
    t.results <- .set_lab(t.results, sig(x), names(t.results), "Contrast all result table")
    t.results <- setLegend(t.results, glue::glue("差异分析 {names(t.results)} 的全部基因的结果表格。"))
    res <- .stat_by_logfc(t.sigResults, "log2FoldChange", get = x$idcol)
    x <- snapAdd(x, "显著基因统计：\n\n{res$snap}\n\n\n")
    feature(x) <- res$sets
    x <- tablesAdd(x, t.results, t.sigResults)
    x <- methodAdd(
      x, "对 DESeq 函数处理过后的数据获取差异基因。筛选差异表达基因（DEGs）的标准为 {use} &lt; {cut.p} 且 |log2(FoldChange)| &gt; {cut.fc}。"
    )
    return(x)
  })

setMethod("step3", signature = c(x = "job_deseq2"),
  function(x, top = 10, group = "group"){
    step_message("Draw plots.")
    plots <- sapply(
      names(x@tables$step2$t.results), simplify = FALSE,
      function(name) {
        tops <- unlist(lapply(feature(x)[[ name ]], head, n = top))
        p.volcano <- plot_volcano(
          x@tables$step2$t.results[[name]], x$idcol, use = "padj", fc = 1,
          use.fc = "log2FoldChange", HLs = tops
        )
        p.volcano <- .set_lab(
          p.volcano, sig(x), name, "volcano of DEGs"
        )
        p.volcano <- setLegend(p.volcano, "{name} 差异表达基因的火山图")
        data <- SummarizedExperiment::assay(x$vst)
        data <- data[ rownames(data) %in% tops, , drop = FALSE ]
        data <- t(scale(t(data)))
        data <- pmin(pmax(data, -2), 2)
        p.hp <- pheatmap::pheatmap(data,
          annotation_col = data.frame(SummarizedExperiment::colData(x$vst))[, group, drop = FALSE],
          show_rownames = TRUE,
          cluster_rows = TRUE,
          cluster_cols = FALSE,
          scale = "none")
        p.hp <- .set_lab(
          wrap_scale(
            p.hp, ncol(data), nrow(data), size = .2
          ), sig(x), name, glue::glue("heatmap of top {top} DEGs")
        )
        p.hp <- setLegend(p.hp, "Top {top} 差异表达基因热图。")
        return(namel(p.hp, p.volcano, tops))
      }
    )
    x@plots[[ 3 ]] <- list(topDegs = plots)
    return(x)
  })

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
  function(x, ref, ref.use = "guess", which = 1L, run_roc = TRUE,
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
    fakeLmJob <- focus(
      fakeLmJob, ref = ref, ref.use = ref.use,
      .name = .name, which = which, data.which = data.which, sig = sig, 
      , use = use, test = test, ...
    )
    where <- paste0("focusedDegs_", .name)
    if (identical(clear, "auto")) {
      clear <- !is.null(x[[ where ]])
    }
    resStat <- x[[ where ]] <- fakeLmJob[[ where ]]
    if (run_roc) {
      x[[where]]$roc <- roc <- sapply(ref, simplify = FALSE,
        function(gene) {
          .evaluate_by_ROC(
            elist$E, elist$targets, gene, control = level.roc[2], sig = sig(x)
          )
        })
    }
    snap <- fakeLmJob@snap[[ paste0("step", .name) ]]
    x <- snapAdd(x, snap, step = .name, add = !clear)
    if (run_roc) {
      snaps <- bind(
        paste0("- ", vapply(roc, function(x) x$snap, character(1)), "\n"), co = "\n"
      )
      x <- snapAdd(x, "\n\n{snaps}", step = .name)
    }
    return(x)
  })

setMethod("asjob_iobr", signature = c(x = "job_deseq2"),
  function(x, idType = "Symbol", source = c("local", "biomart"), ...)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
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
    glue::glue("为 {gene} ROC 曲线。")
  )
  snap <- glue::glue(
    paste0(
      "{gene} 的 ROC 分析结果，AUC = {round(auc, 3)}，AUC 的95%置信区间 = [{round(ci[1], 3)}, {round(ci[3], 3)}]；",
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

.stat_by_logfc <- function(lst, col = "log2FoldChange", get = "rownames", cut = 0)
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
    paste0("- 对比组 ", names(res.snap), "，", res.snap), collapse = "\n"
  )
  return(list(snap = snap, sets = sets))
}

setMethod("set_remote", signature = c(x = "job_deseq2"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
