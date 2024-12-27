# ==========================================================================
# workflow of limma
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_limma <- setClass("job_limma", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(
      "https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html"),
    cite = "[@LimmaPowersDiRitchi2015; @EdgerDifferenChen]",
    method = "R package `Limma` and `edgeR` used for differential expression analysis",
    params = list(isTcga = FALSE, normed = FALSE, rna = FALSE),
    tag = "rna:diff",
    analysis = "Limma 差异分析"
    ))

job_limma_normed <- function(data, metadata, genes = NULL) {
  .check_columns(metadata, c("sample", "group"))
  metadata <- dplyr::slice(metadata, match(colnames(data), sample))
  if (is(data, "tbl_df")) {
    if (is.character(data[[1]])) {
      message("Convert as data.frame, and use the first column set as rownames.")
      rownames <- data[[1]]
      data <- data.frame(data, check.names = FALSE)
      rownames(data) <- rownames
    } else {
      stop("The first column not seems be 'ID' (character).")
    }
  }
  data <- dplyr::select(data, dplyr::all_of(metadata$sample))
  if (!identical(colnames(data), metadata$sample)) {
    stop("!identical(colnames(data), metadata$sample)")
  }
  if (!is.null(genes)) {
    message("Be careful, the first columns of `genes` were considered as ID columns.")
    if (any(duplicated(genes[[1]]))) {
      stop("any(duplicated(genes[[1]]))")
    }
    genes <- dplyr::slice(genes, match(rownames(data), genes[[1]]))
  }
  .job_limma(object = data, params = list(metadata = metadata, isTcga = FALSE, normed = TRUE, genes = genes))
}

setMethod("regroup", signature = c(x = "job_limma"),
  function(x, ...){
    modify_job_limma_meta(x, ..., fun = dplyr::mutate)
  })

modify_job_limma_meta <- function (x, ..., fun) {
  if (!is.null(x@object$samples)) {
    x@object$samples <- fun(x@object$samples, ...)
  }
  if (!is.null(x$metadata)) {
    x$metadata <- fun(x$metadata, ...)
  }
  if (!is.null(x$.metadata)) {
    x$.metadata <- fun(x$.metadata, ...)
  }
  if (!is.null(x$normed_data$targets)) {
    x$normed_data$targets <- fun(x$normed_data$targets, ...)
  }
  x
}

setMethod("filter", signature = c(x = "job_limma"),
  function(x, ..., type = c("gene", "metadata")){
    type <- match.arg(type)
    if (type == "gene") {
      message("Filter genes via `x@object$genes`.")
      message("Before Dim: ", paste0(dim(object(x)), collapse = ", "))
      data <- object(x)$genes
      data$...seq <- seq_len(nrow(data))
      keep <- dplyr::filter(data, ...)$...seq
      object(x) <- object(x)[keep, ]
      message("After Dim: ", paste0(dim(object(x)), collapse = ", "))
      return(x)
    } else if (type == "metadata") {
      metadata <- if (x$normed) x$metadata else x@object$samples
      snap <- snap(trace_filter(metadata, ...))
      x <- modify_job_limma_meta(x, ..., fun = dplyr::filter)
      if (x@step < 1L) {
        message(
          glue::glue("before filter, dim of object(x)$counts: {bind(dim(object(x)$counts))}")
        )
        object(x)$counts <- object(x)$counts[ , colnames(object(x)$counts) %in% object(x)$samples$sample ]
        message(
          glue::glue("after filter, dim of object(x)$counts: {bind(dim(object(x)$counts))}")
        )
      }
      x <- snapAdd(x, "{snap}")
      return(x)
    }
  })

job_limma <- function(DGEList, rna = TRUE)
{
  if (!is(DGEList, 'DGEList'))
    stop("is(DGEList, 'DGEList') == FALSE")
  x <- .job_limma(object = DGEList)
  x$.metadata <- tibble::as_tibble(dplyr::relocate(DGEList$samples, sample, group))
  x$rna <- rna
  return(x)
}

setMethod("step0", signature = c(x = "job_limma"),
  function(x){
    step_message("Prepare your data with function `job_limma`.
      "
    )
  })

.get_meta <- function(x, what) {
  if (x$normed) x$metadata[[ what ]] else x@object$samples[[ what ]]
}

.guess_symbol <- function(x) {
  if (x$isTcga) "gene_name" else if (x$rna) "hgnc_symbol" else "GENE_SYMBOL"
}

.guess_formula <- function(envir = parent.frame(1)) {
  thisEnv <- environment()
  lapply(c("group", "batch", "pairs"),
    function(name) {
      assign(name, get(name, envir), envir = thisEnv)
  })
  if (is.null(batch) && is.null(pairs)) {
    "~ 0 + group"
  } else if (is.null(batch) && !is.null(pairs)) {
    "~ 0 + group + pairs"
  } else if (!is.null(batch) && is.null(pairs)) {
    "~ 0 + group + batch"
  } else if (!is.null(batch) && !is.null(pairs)) {
    "~ 0 + group + batch + pairs"
  }
}

setMethod("step1", signature = c(x = "job_limma"),
  function(x,
    group = .get_meta(x, "group"), batch = .get_meta(x, "batch"),
    pairs = .get_meta(x, "pairs"),
    formula = .guess_formula(), design = mx(as.formula(formula)),
    min.count = 10,
    no.rna_filter = if (x$normed) TRUE else FALSE,
    no.rna_norm = no.rna_filter,
    no.array_norm = "guess",
    norm_vis = FALSE, pca = FALSE, data_type = c(
      "count", "cpm", "tpm"
    ))
  {
    step_message("Preprocess expression data.")
    if (x$rna) {
      x <- methodAdd(x, "以 R 包 `edgeR` ({packageVersion('edgeR')}) {cite_show('EdgerDifferenChen')} 对数据预处理。")
      x <- snapAdd(x, "以 `edgeR` 将{x$project} RNA-seq 数据标准化 (详见方法章节)。")
    }
    message(glue::glue("Use formula: {formula}"))
    data_type <- match.arg(data_type)
    plots <- list()
    s.com <- try_snap(if (x$normed) x$metadata else object(x)$sample, "group", "sample")
    x <- snapAdd(x, "样本分组：{s.com}。", FALSE)
    if (x$rna || x$isTcga) {
      message("Data from RNA-seq.")
      if (!no.rna_filter && data_type == "count") {
        object(x) <- filter_low.dge(object(x), group, min.count = min.count)
        x <- methodAdd(x, "以 `edgeR::filterByExpr` 过滤 count 数量小于 {min.count} 的基因。", TRUE)
        p.filter <- wrap(attr(object(x), "p"), 8, 3)
        p.filter <- .set_lab(p.filter, sig(x), "Filter low counts")
        plots <- c(plots, namel(p.filter))
      } else {
        message("Skip from filtering.")
      }
      if (!no.rna_norm) {
        object(x) <- norm_genes.dge(object(x), design, vis = norm_vis, data_type = data_type)
        if (data_type == "count") {
          x <- methodAdd(x, "以 `edgeR::calcNormFactors`，`limma::voom` 转化 {data_type} 数据为 log2 counts-per-million (logCPM)。")
        } else if (data_type == "tpm") {
          x$tpm_use_trend <- TRUE
        }
        if (norm_vis && data_type == "count") {
          if (length(x@object$targets$sample) < 50) {
            p.norm <- wrap(attr(object(x), "p"), 6, max(c(length(x@object$targets$sample) * .6, 10)))
          } else {
            p.norm <- wrap(attr(object(x), "p"))
          }
          p.norm <- .set_lab(p.norm, sig(x), "Normalization")
          x@params$p.norm_data <- p.norm@data$data
        } else {
          p.norm <- NULL
        }
        plots <- c(plots, namel(p.norm))
        x@params$normed_data <- object(x)
      } else {
        if (x$rna && is(object(x), "df")) {
          x$normed_data <- list(
            genes = if (is.null(x$genes)) data.frame(rownames = rownames(object(x))) else x$genes,
            targets = x$metadata,
            E = object(x)
          )  
        } else {
          x$normed_data <- object(x)
        }
      }
    } else {
      message("Data from Microarray.")
      if (no.array_norm == "guess") {
        range <- range(object(x)$counts)
        if (all(range > 0) && range[2] > 100) {
          message("May be raw expression dataset.")
          no.array_norm <- FALSE
        } else if (range[1] == 0 && range[2] > 100) {
          message("Min expression equal to 0, add prior value: 1.")
          object(x)$counts <- object(x)$counts + 1
          no.array_norm <- FALSE
        } else {
          no.array_norm <- TRUE
        }
      }
      if (!no.array_norm) {
        x <- methodAdd(x, "使用 `log2` 和 `limma::normalizeBetweenArrays` 对数据标准化。")
        object(x)$counts <- e(limma::normalizeBetweenArrays(log2(object(x)$counts)))
        x$normed_data <- new("EList",
          list(genes = object(x)$genes, targets = object(x)$samples, E = object(x)$counts)
        )
      }
    }
    if (pca) {
      pca <- pca_data.long(as_data_long(object(x)))
      p.pca <- plot_andata(pca)
      plots <- c(plots, namel(p.pca))
    }
    if (length(plots)) {
      x@plots[[ 1 ]] <- plots
    }
    x@params$group <- group
    x@params$design <- design
    snap <- ""
    if (!is.null(batch)) {
      snap <- paste0("(Batch: ", try_snap(batch), ")")
    }
    x <- snapAdd(x, "以 公式 {formula} 创建设计矩阵 (design matrix) {snap}。")
    x$.metadata <- .set_lab(x$.metadata, sig(x), "metadata of used sample")
    x$metadata <- .set_lab(x$metadata, sig(x), "metadata of used sample")
    return(x)
  })

setMethod("step2", signature = c(x = "job_limma"),
  function(x, ..., contrasts = NULL, block = NULL, use = c("adj.P.Val", "P.Value"),
    use.cut = .05, cut.fc = 1,
    label = .guess_symbol(x), batch = FALSE, HLs = NULL)
  {
    step_message("Difference test.")
    use <- match.arg(use)
    if (is.null(contrasts)) {
      if (!length(alist(...))) {
        contr <- NULL
      } else {
        contr <- limma::makeContrasts(..., levels = x@params$design)
      }
    } else {
      contr <- limma::makeContrasts(contrasts = contrasts, levels = x@params$design)
    }
    s.contr <- paste0(gsub('-', 'vs', colnames(contr)), collapse = ', ')
    x <- methodAdd(x, "以 `limma` ({packageVersion('limma')}) {cite_show('LimmaLinearMSmyth2005')} 差异分析。")
    if (x$rna) {
      x <- methodAdd(x, "分析方法参考 <https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html>。")
    }
    x <- methodAdd(x, "创建设计矩阵，对比矩阵，差异分析：{s.contr}。")
    x <- snapAdd(x, "差异分析：{s.contr}。(若 A vs B，则为前者比后者，LogFC 大于 0 时，A 表达量高于 B)。")
    ## here, remove batch effect
    ## limma::removeBatchEffect
    if (batch) {
      message("Note: `limma::removeBatchEffect` would convert the object class.")
      object(x) <- e(limma::removeBatchEffect(object(x),
          batch = object(x)$targets$batch, design = x@params$design, group = x$targets$group
          ))
    }
    if (!is.null(x$tpm_use_trend) && x$tpm_use_trend) {
      # see: https://support.bioconductor.org/p/98820/
      message("As data_type == 'tpm', limma::eBayes will use params: trend = TRUE.")
      trend <- TRUE
    } else {
      trend <- FALSE
    }
    object(x) <- diff_test(object(x), x@params$design, contr, block, trend = trend)
    x <- methodAdd(x, "使用 `limma::lmFit`, `limma::contrasts.fit`, `limma::eBayes` 拟合线形模型。")
    plots <- list()
    if (!is.null(contr)) {
      tops <- extract_tops(object(x), use = use, use.cut = use.cut, cut.fc = cut.fc)
      tops <- lapply(
        tops, dplyr::relocate,
        logFC, !!rlang::sym(use)
      )
      x <- methodAdd(x, "以 `limma::topTable` 提取所有结果，并过滤得到 {use} 小于 {use.cut}，|Log2(FC)| 大于 {cut.fc} 的统计结果。")
      if (x$normed) {
        if (!is.null(x$genes)) {
          tops <- lapply(tops,
            function(obj) {
              obj <- map(obj, colnames(obj)[1], x$genes, colnames(x$genes)[1], label, col = label)
              if (!is.null(x$from_scfea)) {
                ## the genes belong to the module
                obj <- map(obj, colnames(obj)[1], x$genes, colnames(x$genes)[1], "gene", col = "gene")
                dic <- nl(x$compounds_annotation$name, x$compounds_annotation$kegg, TRUE)
                obj <- dplyr::mutate(obj, compounds = strsplit(name, " -> "),
                  kegg = lapply(compounds,
                    function(x) {
                      dplyr::recode(x, !!!dic)
                    })
                )
              }
              obj
            })
        }
      }
      tops <- .set_lab(tops, sig(x), paste("data", gs(names(tops), "-", "vs")))
      lab(tops) <- paste(sig(x), "DEGs data")
      if (length(tops) >= 2) {
        lst <- lapply(tops, function(x) x[[ label ]])
        names(lst) <- gs(names(lst), "\\s*-\\s*", " vs ")
        p.contrast_cols <- new_col(lst = lst)
        p.contrast_cols <- .set_lab(p.contrast_cols, sig(x), "All Difference Feature of contrasts")
        plots <- c(plots, namel(p.contrast_cols))
      }
    } else {
      tops <- NULL
    }
    print(glue::glue("use: {label}"))
    p.volcano <- lapply(tops, plot_volcano, label = label, use = use, fc = cut.fc, seed = x$seed, HLs = HLs)
    p.volcano <- lapply(p.volcano,
      function(p) {
        p <- wrap(p, 5, 4)
        attr(p, "lich") <- new_lich(nl(c(paste0(use, " cut-off"), "Log2(FC) cut-off"), c(use.cut, cut.fc)))
        p
      }
    )
    p.volcano <- .set_lab(p.volcano, sig(x), gs(names(p.volcano), "-", "vs"))
    lab(p.volcano) <- paste(sig(x), "volcano plot")
    plots <- c(plots, namel(p.volcano))
    tables <- namel(tops)
    if (!is.null(x$from_scfea)) {
      message("Gene names split by '-'.")
      belong.flux <- reframe_col(dplyr::select(ref, gene, name), "gene",
        function(x) unlist(strsplit(unlist(x), "-")))
      belong.flux <- dplyr::relocate(belong.flux, gene, Metabolic_flux = name)
      belong.flux <- .set_lab(belong.flux, sig(x), "Differential metabolic flux related Genes")
      tables <- c(tops, namel(belong.flux))
    }
    x@tables[[ 2 ]] <- tables
    x@plots[[ 2 ]] <- plots
    return(x)
  })

setMethod("group", signature = c(x = "job_limma"),
  function(x, group)
  {
    stopifnot(x@step >= 1L)
    object(x) <- x$normed_data
    group <- x$normed_data$targets[[ group ]]
    stopifnot(!is.null(group))
    object(x)$targets$group <- group
    x$normed_data$targets$group <- group
    design <- mx(~ 0 + group)
    x$design <- design
    x@step <- 1L
    return(x)
  })

setMethod("step3", signature = c(x = "job_limma"),
  function(x, names = NULL, use = "all", use.gene = .guess_symbol(x),
    fun_filter = rm.no, trunc = NULL)
  {
    step_message("Sets intersection.")
    tops <- x@tables$step2$tops
    if (use != "all") {
      tops <- tops[ use ]
    }
    if (!is.null(names)) {
      names(tops) <- names
    }
    tops <- lapply(tops,
      function(data){
        if (is.null(data[[ use.gene ]])) {
          stop(
            glue::glue('is.null(data[[ use.gene ]]), you specified {use.gene}, but found NULL')
          )
        }
        up <- dplyr::filter(data, logFC > 0)[[ use.gene ]]
        down <- dplyr::filter(data, logFC < 0)[[ use.gene ]]
        lst <- list(up = up, down = down)
        lapply(lst, fun_filter)
      })
    tops <- lapply(tops, function(x) lapply(x, gname))
    if (length(tops) == 1) {
      x <- snapAdd(x, "上调或下调 DEGs 统计：{try_snap(tops[[1]])}")
      x$.feature <- tops
    } else {
      names(tops) <- gs(names(tops), "-", "vs")
      s.com <- vapply(tops, try_snap, character(1))
      s.com <- paste0("- ", names(s.com), "：", s.com, "。")
      s.com <- paste0(s.com, collapse = "\n")
      x <- snapAdd(x, "各组差异分析 DEGs 统计：\n\n {s.com}\n\n")
      sfun <- function(which) length(unique(unlist(lapply(tops, function(x) x[[ which ]]))))
      x <- snapAdd(x, "所有上调 DEGs 共 {sfun('up')} 个，所有下调 DEGs 共 {sfun('down')} 个。")
      x <- snapAdd(x, "所有非重复基因共 {length(unique(unlist(tops)))} 个。")
      tops <- unlist(tops, recursive = FALSE)
      x$.feature <- x$sets_intersection <- tops
      message("The guess use dataset combination of:\n",
        "\t ", names(tops)[1], " %in% ", names(tops)[4], "\n",
        "\t ", names(tops)[2], " %in% ", names(tops)[3])
      x$guess_use <- unique(c(
          intersect(tops[[ 1 ]], tops[[ 4 ]]),
          intersect(tops[[ 2 ]], tops[[ 3 ]])
          ))
      p.hp <- plot_genes_heatmap.elist(x$normed_data, unlist(tops), use.gene)
      p.hp <- .set_lab(wrap(p.hp), sig(x), "Heatmap of DEGs")
      p.sets_intersection <- new_upset(lst = tops, trunc = trunc)
      p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "Difference", "intersection")
      x@plots[[ 3 ]] <- namel(p.sets_intersection, p.hp)
    }
    return(x)
  })

plot_genes_heatmap.elist <- function(normed_data, degs, use = "hgnc_symbol") {
  degs <- unique(degs)
  normed_data <- normed_data[ normed_data$genes[[ use ]] %in% degs, ]
  rownames(normed_data) <- normed_data$genes[[ use ]]
  metadata <- dplyr::select(as_tibble(normed_data$targets), sample, group)
  p.hp <- plot_genes_heatmap(normed_data$E, metadata)
  return(p.hp)
}

plot_genes_heatmap <- function(data, metadata) {
  data <- dplyr::rename(as_tibble(data), genes = rownames)
  data <- tidyr::pivot_longer(data, !genes, names_to = "sample", values_to = "expression")
  data <- tbmerge(data, metadata, by = "sample", all.x = TRUE)
  maxBreak <- max(ceiling(abs(range(data$expression))))
  # ComplexHeatmap::Heatmap
  tidyHeatmap::heatmap(dplyr::group_by(data, group), genes, sample, expression,
    palette_value = fun_color(-maxBreak, maxBreak), show_column_names = FALSE)
}

.guess_intersect <- function(tops) {
  unique(c(
      intersect(tops[[ 1 ]], tops[[ 4 ]]),
      intersect(tops[[ 2 ]], tops[[ 3 ]])
      ))
}

setMethod("clear", signature = c(x = "job_limma"),
  function(x, save = TRUE, suffix = NULL){
    if (save)
      saveRDS(x, paste0(substitute(x, parent.frame(1)), x@step, suffix, ".rds"))
    object(x) <- NULL
    x@params$normed_data <- NULL
    return(x)
  })

setMethod("map", signature = c(x = "job_limma"),
  function(x, ref, ref.use = .guess_symbol(x),
    group = NULL, group.use = "group", pvalue = TRUE){
    object <- x@params$normed_data
    if (identical(class(object), "list")) {
      object <- new("EList", object)
    }
    object <- object[!duplicated(object$genes[[ ref.use ]]), ]
    rownames(object) <- object$genes[[ ref.use ]]
    object <- object[rownames(object) %in% ref, ]
    if (any(duplicated(rownames(object)))) {
      object <- object[ !duplicated(rownames(object)), ]
    }
    if (!is.null(group)) {
      object <- object[, object$targets[[ group.use ]] %in% group]
    }
    data <- tibble::as_tibble(t(object$E))
    data$group <- object$targets$group
    data <- tidyr::gather(data, var, value, -group)
    p <- wrap(.map_boxplot2(data, pvalue))
    p <- .set_lab(p, sig(x), "wilcox test for some genes")
    snap(p) <- "以 wilcox.test 检验少量基因 ({bind(ref)}) 的表达。"
    feature(p) <- ref
    p
  })

plot_volcano <- function(top_table, label = "hgnc_symbol", use = "adj.P.Val",
  fc = .3, seed = 1, HLs = NULL, use.fc = "logFC", label.fc = "log2(FC)",
  label.p = paste0("-log10(", use, ")"))
{
  set.seed(seed)
  if (!any(label == colnames(top_table))) {
    if (any("rownames" == colnames(top_table)))
      label <- "rownames"
  }
  data <- dplyr::select(top_table, !!rlang::sym(label), !!rlang::sym(use.fc), !!rlang::sym(use))
  data <- dplyr::mutate(data,
    change = ifelse(!!rlang::sym(use.fc) > abs(fc), "up",
      ifelse(!!rlang::sym(use.fc) < -abs(fc), "down", "stable"))
  )
  pal <- color_set2()
  p <- ggplot(data, aes(x = !!rlang::sym(use.fc), y = -log10(!!rlang::sym(use)), color = change)) + 
    geom_point(alpha = 0.8, stroke = 0, size = 1.5) + 
    scale_color_manual(values = c("down" = pal[1], "stable" = "grey90", "up" = pal[2])) +
    geom_hline(yintercept = -log10(0.05), linetype = 4, size = 0.8) +
    geom_vline(xintercept = c(-abs(fc), abs(fc)), linetype = 4, size = 0.8) + 
    labs(x = label.fc, y = label.p) + 
    ggrepel::geom_text_repel(
      data = dplyr::distinct(rbind(
        dplyr::slice_min(data, !!rlang::sym(use), n = 10),
        dplyr::slice_max(data, abs(!!rlang::sym(use.fc)), n = 20)
        )), 
      aes(label = !!rlang::sym(label)), size = 3) +
    rstyle("theme") +
    geom_blank()
  if (!is.null(HLs)) {
    data <- dplyr::filter(data, !!rlang::sym(label) %in% HLs)
    p <- p +
      ggrepel::geom_label_repel(data = data,
        aes(x = !!rlang::sym(use.fc), y = -log10(!!rlang::sym(use)), label = !!rlang::sym(label))) +
      geom_point(data = data, aes(x = !!rlang::sym(use.fc), y = -log10(!!rlang::sym(use))),
        color = "darkred", shape = 21, fill = "transparent", size = 3, stroke = .5)
  }
  p
}

setMethod("asjob_wgcna", signature = c(x = "job_limma"),
  function(x, filter_genes = NULL, use = "hgnc_symbol"){
    step_message("Use `x@params$normed_data` converted as job_wgcna.")
    if (is.null(object <- x@params$normed_data))
      stop("is.null(x@params$normed_data)")
    if (is.null(use)) {
      if (any(duplicated(rownames(object$genes)))) {
        stop("The rownames of `object$genes` has duplicated value.")
      }
      rownames(object) <- rownames(object$genes)
    } else {
      if (any(duplicated(object$genes[[ use ]]))) {
        stop(glue::glue("The column of `use`: {use} has duplicated value."))
      }
      rownames(object) <- object$genes[[ use ]]
    }
    if (!is.null(filter_genes)) {
      object <- object[ object$genes[[ use ]] %in% filter_genes, ]
    }
    log_counts <- as_tibble(object$E)
    gene_annotation <- as_tibble(object$genes)
    log_counts[[1]] <- gene_annotation[[1]]
    job_wgcna(select(object$targets, sample, group),
      log_counts, gene_annotation)
  })

setMethod("meta", signature = c(x = "job_limma"),
  function(x, use = "group"){
    if (x@step < 1) {
      metadata <- object(x)$samples
    } else {
      metadata <- x$normed_data$targets
    }
    x@params$p.meta <- new_pie(metadata[[ use ]])
    return(x)
  })

setMethod("tops", signature = c(x = "job_limma"),
  function(x, key = 1L, col = "hgnc_symbol"){
    features <- x@tables$step2$tops[[key]][[col]]
    features <- features[!is.na(features) & features != ""]
    features <- unlist(strsplit(features, " /// "), use.names = FALSE)
    gs(features, "\\.[0-9]*$", "")
  })

setMethod("cal_corp", 
  signature = c(x = "job_limma", y = "job_limma"),
  function(x, y, from, to, names = NULL, use = .guess_symbol(x),
    theme = NULL, HLs = NULL, mode = c("heatmap", "linear"))
  {
    message("Filter out others.")
    eval(use)
    sigs <- c(x@sig, y@sig)
    x <- x$normed_data
    x <- x[x$genes[[ use ]] %in% from, ]
    if (!nrow(x)) {
      stop('!nrow(x), No any genes found')
    }
    y <- y$normed_data
    y <- y[y$genes[[ use ]] %in% to, ]
    if (!nrow(y)) {
      stop('!nrow(y), No any genes found')
    }
    message("Got common sample.")
    commonSample <- ins(x$targets$sample, y$targets$sample)
    if (!length(commonSample)) {
      stop("No any common sample.")
    }
    x <- x[, x$targets$sample %in% commonSample]
    y <- y[, y$targets$sample %in% commonSample]
    message("`rbind` for genes.")
    x$E <- rbind(x$E, y$E)
    x$genes <- rbind(x$genes, y$genes)
    newjob <- .job_limma(params = list(normed_data = x))
    x <- cal_corp(newjob, NULL, from, to, names, use, theme, HLs, mode)
    x$.snap <- c(
      list(step0 = glue::glue("将两组相同样品来源的数据集 (dataset: {bind(sigs)})) 关联分析。")),
      step1 = x$.snap[[1]])
    x
  })

setMethod("cal_corp", signature = c(x = "job_limma", y = "NULL"),
  function(x, y, from, to, names = NULL, use = .guess_symbol(x),
    theme = NULL, HLs = NULL, mode = c("heatmap", "linear"))
  {
    mode <- match.arg(mode)
    data <- as_tibble(x@params$normed_data$E)
    anno <- as_tibble(x@params$normed_data$genes)
    if (!any(colnames(anno) == use)) {
      if (!is.null(x$genes)) {
        if (any(colnames(x$genes) == use)) {
          message("Use '", use, "' in `x$genes`.")
          anno <- map(anno, colnames(anno)[1], x$genes, colnames(x$genes)[1], use, col = use)
        } else {
          stop("`use` not found.")
        }
      }
    }
    if (mode == "heatmap") {
      lst <- .cal_corp.elist(data, anno, use, unique(from), unique(to), names, HLs = HLs, fast = TRUE)
      if (length(unique(from)) > 1 && length(unique(to)) > 1) {
        lst$hp <- .set_lab(wrap(lst$hp), sig(x), theme, "correlation heatmap")
        lst$sig.corp <- .set_lab(lst$sig.corp, sig(x), theme, "significant correlation")
      }
    } else if (mode == "linear") {
      lst <- .cal_corp.elist(data, anno, use, unique(from), unique(to), names, HLs = HLs, fast = FALSE)
    }
    x <- .job(params = lst)
    fun <- function(x) length(unique(x))
    if (identical(from, to)) {
      x <- snapAdd(x, "将基因集 ({fun(from)}) 关联分析，")
    } else {
      x <- snapAdd(x, "将基因集 (A -> B) (A:{fun(from)}, B:{fun(to)}) 关联分析，")
    }
    x <- snapAdd(x, "共得到 {nrow(lst$sig.corp)} 个显著的基因对 (P &lt; 0.05)。")
    return(x)
  })

.cal_corp.elist <- function(data, anno, use, from, to, names, HLs = NULL, fast = TRUE)
{
  if (is.null(data$rownames)) {
    data <- as_tibble(data)
  }
  data$rownames <- anno[[ use ]]
  colnames(data)[1] <- use
  data <- dplyr::mutate(data, symbol = gname(!!rlang::sym(use)))
  if (use != "symbol") {
    data <- dplyr::select(data, -!!rlang::sym(use))
    data <- dplyr::relocate(data, symbol)
  }
  lst <- lapply(list(from, to),
    function(set) {
      set <- gname(set)
      data <- dplyr::filter(data, symbol %in% dplyr::all_of(set))
      dplyr::distinct(data, symbol, .keep_all = TRUE)
    })
  # heatmap
  if (is.null(names)) {
    corp <- cal_corp(lst[[1]], lst[[2]], "From", "To", trans = TRUE, fast = fast)
  } else {
    corp <- cal_corp(lst[[1]], lst[[2]], names[[1]], names[[2]], trans = TRUE, fast = fast)
  }
  if (fast) {
    sig.corp <- dplyr::filter(tibble::as_tibble(corp), sign != "-")
    if (length(from) > 1 && length(to) > 1) {
      hp <- new_heatdata(corp)
      hp <- callheatmap(hp, HLs = HLs)
      namel(corp, sig.corp, hp)
    } else {
      namel(corp, sig.corp)
    }
  } else {
    corp
    # regression line
  }
}

setMethod("vis", signature = c(x = "corp"),
  function(x, group = NULL, facet = ".id")
  {
    x <- as_tibble(x)
    x <- dplyr::mutate(x, .id = paste0(From, "_", To))
    facet <- match.arg(facet)
    theme <- rstyle("theme")
    fun_plot <- function(x) {
      models <- x$model
      names(models) <- x$.id
      data <- frbind(models, idcol = ".id")
      if (facet == ".id") {
        facet <- ggplot2::facet_wrap(~ .id, scales = "free")
      }
      anno <- dplyr::mutate(x, model = lapply(model,
          function(x) {
            c(j = max(zoRange(x[[ "j" ]], 1.3), na.rm = TRUE),
              i = min(x[[ "i" ]], na.rm = TRUE))
          }),
        y = vapply(model, function(x) x[[ "j" ]], numeric(1)),
        x = vapply(model, function(x) x[[ "i" ]], numeric(1))
      )
      p <- ggplot(data, aes(y = j, x = i)) +
        geom_point() +
        stat_smooth(method = "lm", col = "red") +
        facet +
        geom_text(data = anno,
          aes(x = x, y = y,
            label = paste0("Cor = ", round(cor, 2), "\n", "P-value = ", signif(pvalue, 2))),
          size = 3, hjust = 0, vjust = 1) +
        theme
      p <- as_grob(p)
      layout <- dplyr::filter(p$layout, grpl(name, "^panel"))
      f.w <- max(as.integer(strx(layout$name, "[0-9]+")))
      f.h <- max(as.integer(strx(layout$name, "[0-9]+$")))
      wrap(p, 4 * f.w, 4 * f.h)
    }
    if (!is.null(group)) {
      lst.p <- lapply(split(x, x[[group]]), fun_plot)
      lst.p <- .set_lab(lst.p, "Linear regression of", names(lst.p))
    } else {
      p <- fun_plot(x)
      .set_lab(p, "Linear regression")
    }
  })

setMethod("getsub", signature = c(x = "job_limma"),
  function(x, tnbc = FALSE){
    if (tnbc) {
      if (x$isTcga && identical(x$project, "TCGA-BRCA")) {
        cli::cli_alert_info("get_data.rot2016")
        pb.rot <- get_data.rot2016()
        isTnbc <- gs(object(x)$samples$sample, "[A-Z]$", "") %in%
          dplyr::filter(object(pb.rot), TNBC == "YES")$TCGA_SAMPLE
        isNormal <- object(x)$samples$isTumor == "normal"
        object(x) <- object(x)[, isTnbc | isNormal ]
        object(x)$samples$group <- object(x)$samples$isTumor
        message("Change group to `isTumor`.")
        x$pb <- pb.rot
        .add_internal_job(pb.rot)
      } else {
        stop("!x$isTcga && identical(x$project, 'TCGA-BRCA')")
      }
      message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    }
    return(x)
  })


expand.cons <- function(...) {
  apply(expand.grid(...), 1, simplify = TRUE,
    function(x){
      paste0(x[1], "-", x[2])
    })
}

filter_low.dge <- function(dge, group., min.count = 10, prior.count = 2) {
  dge$samples$group = group.
  mean <- mean(dge$samples$lib.size) * 1e-6
  median <- median(dge$samples$lib.size) * 1e-6
  cutoff <- log2(min.count / median + 2 / mean)
  ## raw
  raw_dge <- dge
  raw_dge$counts <- edgeR::cpm(raw_dge, log = TRUE, prior.count = prior.count)
  ## filter
  keep.exprs <- e(edgeR::filterByExpr(dge, group = group., min.count = min.count))
  pro_dge <- dge <- e(edgeR::`[.DGEList`(dge, keep.exprs, , keep.lib.sizes = FALSE))
  pro_dge$counts <- e(edgeR::cpm(pro_dge, log = TRUE))
  ## plot
  data <- list(Raw = as_data_long(raw_dge), Filtered = as_data_long(pro_dge))
  data <- data.table::rbindlist(data, idcol = TRUE)
  data <- dplyr::select(data, .id, sample, value)
  p <- ggplot(data) +
    geom_density(aes(x = value, color = sample), alpha = .7) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    facet_wrap(~ factor(.id, c("Raw", "Filtered")), scales = "free_y") +
    labs(x = "Log2-cpm", y = "Density") +
    rstyle("theme")
  if (length(unique(data$sample)) > 50) {
    p <- p + theme(legend.position = "none")
  }
  attr(dge, "p") <- p
  dge
}

norm_genes.dge <- function(dge, design, prior.count = 2, ..., vis = TRUE,
  data_type = c("count", "cpm", "tpm"))
{
  ## raw
  data_type <- match.arg(data_type)
  if (data_type == "count") {
    raw_dge <- dge
    raw_dge$counts <- edgeR::cpm(raw_dge, log = TRUE, prior.count = prior.count)
    ## pro
    dge <- e(edgeR::calcNormFactors(dge, method = "TMM"))
    pro_dge <- dge <- e(limma::voom(dge, design, ...))
  } else {
    dge$counts <- log2(dge$counts)
  }
  ## data long
  if (data_type == "count" && vis) {
    cli::cli_alert_info("as_data_long")
    data <- list(Raw = as_data_long(raw_dge), Normalized = as_data_long(pro_dge))
    data <- data.table::rbindlist(data, idcol = TRUE, fill = TRUE)
    data <- dplyr::select(data, .id, sample, value)
    if (length(unique(data$sample)) < 50) {
      p <- ggplot(data) +
        geom_boxplot(aes(x = sample, y = value),
          outlier.color = "grey60", outlier.size = .5) +
        coord_flip() +
        facet_wrap(~ factor(.id, c("Raw", "Normalized"))) +
        labs(x = "Sample", y = "Log2-cpm") +
        rstyle("theme")
    } else {
      p <- plot_median_expr_line(data)
    }
  } else {
    p <- NULL
  }
  attr(dge, "p") <- p
  dge
}

diff_test <- function(x, design, contr = NULL, block = NULL, trend = FALSE)
{
  if (!is.null(block)){
    dupcor <- e(limma::duplicateCorrelation(x, design, block = block))
    cor <- dupcor$consensus.correlation
    message("## Within-donor correlation:", cor)
  } else {
    cor <- NULL
  }
  if (is(x, "DGEList")) {
    x <- new(
      "EList", list(E = x$counts, targets = x$samples, genes = x$genes)
    )
  }
  fit <- e(limma::lmFit(x, design, block = block, correlation = cor))
  if (!is.null(contr)) {
    fit <- e(limma::contrasts.fit(fit, contrasts = contr))
  }
  # https://liuyujie0136.gitbook.io/sci-tech-notes/bioinformatics/p-value
  fit <- e(limma::eBayes(fit, trend = trend))
  fit
}

extract_tops <- function(x, use = "adj.P.Val", use.cut = 0.05, cut.fc = 0.3){
  res <- e(lapply(seq_len(ncol(x$contrasts)),
    function(coef){
      res <- limma::topTable(x, coef = coef, number = Inf)
      if (!is.null(use.cut) & !is.null(cut.fc)) {
        res <- dplyr::filter(res, !!rlang::sym(use) < use.cut, abs(logFC) > cut.fc)
      }
      as_tibble(res) 
    }))
  names(res) <- colnames(x$contrasts)
  res
}

prepare_expr_data <- function(metadata, counts, genes, message = TRUE) 
{
  ## sort and make name for metadata and counts
  checkDup <- function(x) {
    if (any(duplicated(x))) {
      stop("The value in ID column are duplicated.")
    }
  }
  lapply(list(counts[[ 1 ]], genes[[ 1 ]], metadata[[ 1 ]]), checkDup)
  colnames(metadata) %<>% make.names()
  if (any(colnames(metadata) == "sample")) {
    metadata <- dplyr::relocate(metadata, sample)
  } else {
    metadata <- dplyr::rename(metadata, sample = 1)
  }
  cli::cli_alert_info("Only keep samples record in `metadata`.")
  counts <- dplyr::select(counts, 1, dplyr::all_of(metadata$sample))
  if (ncol(counts) != nrow(metadata) + 1) {
    stop("ncol(counts) != nrow(metadata) + 1")
  }
  cli::cli_alert_info("`make.names` for sample names.")
  metadata$sample %<>% make.names()
  colnames(counts) %<>% make.names()
  ## sort genes
  data_id <- do.call(data.frame, nl(colnames(genes)[1], list(counts[[1]])))
  genes <- dplyr::distinct(genes, !!rlang::sym(colnames(genes)[1]), .keep_all = TRUE)
  genes <- tbmerge(
    data_id, genes,
    by.x = colnames(data_id)[1], by.y = colnames(genes)[1],
    sort = FALSE, all.x = TRUE
  )
  if (ncol(genes) > 1) {
    if (message) {
      message("## The missing genes in `genes`")
      check.na <- dplyr::filter(
        genes, is.na(!!rlang::sym(colnames(genes)[2]))
      )
      print(check.na)
    }
  }
  counts <- data.frame(dplyr::select(counts, -1))
  colnames(counts) <- metadata$sample
  rownames(counts) <- genes[[1]]
  cli::cli_alert_info("Please Check the head count data:")
  print(counts[1:10, 1:5])
  namel(counts, metadata, genes)
}

new_dge <- function(metadata, counts, genes, message = TRUE)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  data <- lst$counts
  if (any(is.na(data))) {
    message("NA value detected in maybe 'counts' data, mutate as zero.")
    data <- dplyr::mutate(data, dplyr::across(dplyr::where(is.numeric), 
      function(x) {
        ifelse(is.na(x), 0, x)
      }))
  }
  e(edgeR::DGEList(data, samples = lst$metadata, genes = lst$genes))
}

new_elist <- function(metadata, counts, genes, message = TRUE)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  elist(list(E = tibble::as_tibble(lst$counts), targets = lst$metadata, genes = lst$genes))
}

# setMethod("snap", 
  # signature = c(x = "job_limma"),
  # function(x, ref, group = "group"){
  #   if (missing(ref) || ref == 1) {
  #     meta <- x$normed_data$targets
  #     if (is.null(meta)) {
  #       meta <- x$metadata
  #     }
  #     if (!is.null(group)) {
  #       metalst <- split(meta$sample, meta[[ group ]])
  #       each <- vapply(names(metalst), function(x) paste0(x, " (", length(metalst[[x]]), ") "), character(1))
  #       group_text <- glue::glue("包含 {paste(each, collapse = ', ')}。")
  #     } else {
  #       group_text <- ""
  #     }
  #     glue::glue("共 {nrow(meta)} 个样本。{group_text}")  
  #   } else if (ref == 2) {
  #     tops <- x@tables$step2$tops
  #     vs <- paste0(gs(names(tops), "-", "vs"), collapse = ", ")
  #     glue::glue("差异分析 {vs} (若 A vs B，则为前者比后者，LogFC 大于 0 时，A 表达量高于 B)")
  #   } else if (ref == 3) {
  #     # tops <- x@tables$step2$tops
  #     # stats <- lapply(tops,
  #     # function(x) {
  #     #   up <- nrow(dplyr::filter(x, LogFC > 0))
  #     #   down <- nrow(dplyr::filter(x, LogFC < 0))
  #     # })
  #     raws <- x@plots$step3$p.sets_intersection$raw
  #     if (!is.null(raws)) {
  #       sums <- lapply(c("\\.up$", "\\.down$"),
  #         function(pattern) length(unique(unlist(raws[ grpl(names(raws), pattern)])))
  #       )
  #       alls <- length(unique(unlist(raws)))
  #     } else {
  #       raws <- x@tables$step2$tops[[1]]
  #       alls <- nrow(raws)
  #       sums <- lapply(1:2, function(n) {
  #         if (n == 1) {
  #           nrow(dplyr::filter(raws, logFC > 0))
  #         } else {
  #           nrow(dplyr::filter(raws, logFC < 0))
  #         }
  #       })
  #     }
  #     glue::glue("所有上调 DEGs 有 {sums[[1]]} 个，下调共 {sums[[2]]}；一共 {alls} 个 (非重复)。")
  #   }
  # })
