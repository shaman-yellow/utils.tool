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
    info = c("|log~2~(FC)| &gt; 0.03, P-value or adjusted P-value &lt; 0.05"),
    cite = "[@LimmaPowersDiRitchi2015; @EdgerDifferenChen]",
    method = "R package `Limma` and `edgeR` used for differential expression analysis",
    params = list(isTcga = F, normed = F)
    ))

job_limma_normed <- function(data, metadata) {
  .check_columns(metadata, c("sample", "group"))
  metadata <- dplyr::slice(metadata, match(colnames(data), sample))
  data <- dplyr::select(data, dplyr::all_of(metadata$sample))
  if (!identical(colnames(data), metadata$sample)) {
    stop("!identical(colnames(data), metadata$sample)")
  }
  .job_limma(object = data, params = list(metadata = metadata, isTcga = F, normed = T))
}

job_limma <- function(DGEList)
{
  if (!is(DGEList, 'DGEList'))
    stop("is(DGEList, 'DGEList') == F")
  .job_limma(object = DGEList)
}

setMethod("step0", signature = c(x = "job_limma"),
  function(x){
    step_message("Prepare your data with function `job_limma`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_limma"),
  function(x,
    group = if (x$normed) x$metadata$group else x@object$samples$group,
    batch = if (x$normed) x$metadata$batch else x@object$samples$batch,
    design = if (is.null(batch)) mx(~ 0 + group) else mx(~ 0 + group + batch),
    min.count = 10,
    no.filter = if (x$normed) T else F,
    no.norm = no.filter,
    norm_vis = F)
  {
    step_message("Preprocess expression data.")
    plots <- list()
    if (!no.filter) {
      object(x) <- filter_low.dge(object(x), group, min.count = min.count)
      p.filter <- wrap(attr(object(x), "p"), 8, 3)
      plots <- c(plots, namel(p.filter))
    }
    if (!no.norm) {
      object(x) <- norm_genes.dge(object(x), design, norm_vis)
      if (length(x@object$targets$sample) < 50) {
        p.norm <- wrap(attr(object(x), "p"), 6, length(x@object$targets$sample) * .6)
      } else {
        p.norm <- wrap(attr(object(x), "p"))
      }
      plots <- c(plots, namel(p.norm))
      x@params$p.norm_data <- p.norm@data$data
      x@params$normed_data <- object(x)
    } else {
      x$normed_data <- list(
        genes = data.frame(rownames = rownames(object(x))),
        targets = metadata,
        E = object(x)
      )
    }
    if (F) {
      pca <- pca_data.long(as_data_long(object(x)))
      p.pca <- plot_andata(pca)
    }
    if (length(plots)) {
      x@plots[[ 1 ]] <- plots
    }
    x@params$group <- group
    x@params$design <- design
    return(x)
  })

setMethod("step2", signature = c(x = "job_limma"),
  function(x, ..., contrasts = NULL, block = NULL, use = c("adj.P.Val", "P.Value"),
    use.cut = .05, cut.fc = .3,
    label = if (x$isTcga) "gene_name" else "hgnc_symbol", batch = F)
  {
    step_message("Difference test.")
    use <- match.arg(use)
    if (is.null(contrasts)) {
      if (length(alist(...)) == 0) {
        contr <- NULL
      } else {
        contr <- limma::makeContrasts(..., levels = x@params$design)
      }
    } else {
      contr <- limma::makeContrasts(contrasts = contrasts, levels = x@params$design)
    }
    ## here, remove batch effect
    ## limma::removeBatchEffect
    if (batch) {
      object(x) <- e(limma::removeBatchEffect(object(x),
          batch = object(x)$targets$batch, design = x@params$design, group = x$targets$group
          ))
    }
    object(x) <- diff_test(object(x), x@params$design, contr, block)
    plots <- list()
    if (!is.null(contr)) {
      tops <- extract_tops(object(x), use = use, use.cut = use.cut, cut.fc = cut.fc)
      tops <- .set_lab(tops, sig(x), paste("data", gs(names(tops), "-", "vs")), "DEGs")
      lab(tops) <- paste(sig(x), "data", "DEGs")
      if (length(tops) >= 2) {
        lst <- lapply(tops, function(x) x[[ label ]])
        names(lst) <- gs(names(lst), "\\s*-\\s*", " vs ")
        p.contrast_cols <- new_col(lst = lst)
        p.contrast_cols <- .set_lab(p.contrast_cols, sig(x), "All DEGs of contrasts")
        plots <- c(plots, namel(p.contrast_cols))
      }
    } else {
      tops <- NULL
    }
    p.valcano <- lapply(tops, plot_valcano, label = label, use = use, fc = cut.fc)
    p.valcano <- lapply(p.valcano, function(p) wrap(p, 5, 4))
    p.valcano <- .set_lab(p.valcano, sig(x), gs(names(p.valcano), "-", "vs"), "DEGs")
    lab(p.valcano) <- paste(sig(x), "volcano plot", "DEGs")
    plots <- c(plots, namel(p.valcano))
    x@tables[[ 2 ]] <- namel(tops)
    x@plots[[ 2 ]] <- plots
    return(x)
  })

setMethod("step3", signature = c(x = "job_limma"),
  function(x, names = NULL, use = "all", use.gene = "hgnc_symbol", fun_filter = rm.no){
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
        up <- dplyr::filter(data, logFC > 0)[[ use.gene ]]
        down <- dplyr::filter(data, logFC < 0)[[ use.gene ]]
        lst <- list(up = up, down = down)
        lapply(lst, fun_filter)
      })
    tops <- unlist(tops, recursive = F)
    x$sets_intersection <- tops
    message("The guess use dataset combination of:\n",
      "\t ", names(tops)[1], " %in% ", names(tops)[4], "\n",
      "\t ", names(tops)[2], " %in% ", names(tops)[3])
    x$guess_use <- unique(c(
      intersect(tops[[ 1 ]], tops[[ 4 ]]),
      intersect(tops[[ 2 ]], tops[[ 3 ]])
    ))
    p.sets_intersection <- new_upset(lst = tops)
    p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "DEGs", "intersection")
    x@plots[[ 3 ]] <- namel(p.sets_intersection)
    return(x)
  })

.guess_intersect <- function(tops) {
  unique(c(
      intersect(tops[[ 1 ]], tops[[ 4 ]]),
      intersect(tops[[ 2 ]], tops[[ 3 ]])
      ))
}

setMethod("clear", signature = c(x = "job_limma"),
  function(x, save = T, suffix = NULL){
    if (save)
      saveRDS(x, paste0(substitute(x, parent.frame(1)), x@step, suffix, ".rds"))
    object(x) <- NULL
    x@params$normed_data <- NULL
    return(x)
  })

setMethod("map", signature = c(x = "job_limma"),
  function(x, ref, ref.use = "hgnc_symbol", group = NULL, group.use = "group", pvalue = T){
    object <- x@params$normed_data
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
    p <- .map_boxplot(data, pvalue)
    p
  })

plot_valcano <- function(top_table, label = "hgnc_symbol", use = "adj.P.Val", fc = .3) {
  if (!any(label == colnames(top_table))) {
    if (any("rownames" == colnames(top_table)))
      label <- "rownames"
  }
  data <- dplyr::select(top_table, !!rlang::sym(label), logFC, !!rlang::sym(use))
  data <- dplyr::mutate(data,
    change = ifelse(logFC > abs(fc), "up",
      ifelse(logFC < -abs(fc), "down", "stable"))
  )
  p <- ggplot(data, aes(x = logFC, y = -log10(!!rlang::sym(use)), color = change)) + 
    geom_point(alpha = 0.8, stroke = 0, size = 1.5) + 
    scale_color_manual(values = c("down" = "#4DBBD5FF", "stable" = "#8491B4FF", "up" = "#DC0000FF")) +
    geom_hline(yintercept = -log10(0.05), linetype = 4, size = 0.8) +
    geom_vline(xintercept = c(-abs(fc), abs(fc)), linetype = 4, size = 0.8) + 
    labs(x = "log2(FC)", y = paste0("-log10(", use, ")")) + 
    ggrepel::geom_text_repel(
      data = distinct(rbind(
        slice_min(data, !!rlang::sym(use), n = 10),
        slice_max(data, abs(logFC), n = 20)
        )), 
      aes(label = !!rlang::sym(label)), size = 3) +
    geom_blank()
  p
}

setMethod("asjob_wgcna", signature = c(x = "job_limma"),
  function(x, filter_genes = NULL, use = "hgnc_symbol"){
    step_message("Use `x@params$normed_data` converted as job_wgcna.")
    if (is.null(object <- x@params$normed_data))
      stop("is.null(x@params$normed_data)")
    rownames(object) <- object$genes[[ use ]]
    log_counts <- as_tibble(object$E)
    gene_annotation <- select(as_tibble(object$genes), -1)
    log_counts[[1]] <- gene_annotation[[1]]
    if (!is.null(filter_genes)) {
      log_counts <- filter(log_counts, rownames %in% filter_genes)
    }
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
    features <- unlist(strsplit(features, " /// "), use.names = F)
    gs(features, "\\.[0-9]*$", "")
  })

setMethod("cal_corp", signature = c(x = "job_limma", y = "NULL"),
  function(x, y, from, to, names = NULL, use = if (x$isTcga) "gene_name" else "hgnc_symbol", theme = NULL)
  {
    data <- as_tibble(x@params$normed_data$E)
    anno <- x@params$normed_data$genes
    lst <- .cal_corp.elist(data, anno, use, from, to, names)
    lst$hp <- .set_lab(wrap(lst$hp), sig(x), theme, "correlation heatmap")
    lst$sig.corp <- .set_lab(lst$sig.corp, sig(x), theme, "significant correlation")
    return(lst)
  })

.cal_corp.elist <- function(data, anno, use, from, to, names)
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
      dplyr::distinct(data, symbol, .keep_all = T)
    })
  if (is.null(names)) {
    corp <- cal_corp(lst[[1]], lst[[2]], "From", "To", trans = T)
  } else {
    corp <- cal_corp(lst[[1]], lst[[2]], names[[1]], names[[2]], trans = T)
  }
  sig.corp <- filter(corp, sign != "-")
  hp <- new_heatdata(corp)
  hp <- callheatmap(hp)
  namel(corp, sig.corp, hp)
}

expand.cons <- function(...) {
  apply(expand.grid(...), 1, simplify = T,
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
  raw_dge$counts <- edgeR::cpm(raw_dge, log = T, prior.count = prior.count)
  ## filter
  keep.exprs <- e(edgeR::filterByExpr(dge, group = group., min.count = min.count))
  pro_dge <- dge <- e(edgeR::`[.DGEList`(dge, keep.exprs, , keep.lib.sizes = F))
  pro_dge$counts <- e(edgeR::cpm(pro_dge, log = T))
  ## plot
  data <- list(Raw = as_data_long(raw_dge), Filtered = as_data_long(pro_dge))
  data <- data.table::rbindlist(data, idcol = T)
  data <- dplyr::select(data, .id, sample, value)
  p <- ggplot(data) +
    geom_density(aes(x = value, color = sample), alpha = .7) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    facet_wrap(~ factor(.id, c("Raw", "Filtered"))) +
    labs(x = "Log2-cpm", y = "Density")
  if (length(unique(data$sample)) > 50) {
    p <- p + theme(legend.position = "none")
  }
  attr(dge, "p") <- p
  dge
}

norm_genes.dge <- function(dge, design, prior.count = 2, fun = limma::voom, ..., vis = T){
  ## raw
  raw_dge <- dge
  raw_dge$counts <- edgeR::cpm(raw_dge, log = T, prior.count = prior.count)
  ## pro
  dge <- e(edgeR::calcNormFactors(dge, method = "TMM"))
  pro_dge <- dge <- fun(dge, design, ...)
  ## data long
  if (vis) {
    cli::cli_alert_info("as_data_long")
    data <- list(Raw = as_data_long(raw_dge), Normalized = as_data_long(pro_dge))
    data <- data.table::rbindlist(data, idcol = T, fill = T)
    data <- dplyr::select(data, .id, sample, value)
    if (length(unique(data$sample)) < 50) {
      p <- ggplot(data) +
        geom_boxplot(aes(x = sample, y = value),
          outlier.color = "grey60", outlier.size = .5) +
        coord_flip() +
        facet_wrap(~ factor(.id, c("Raw", "Normalized"))) +
        labs(x = "Sample", y = "Log2-cpm")
    } else {
      p <- plot_median_expr_line(data)
    }
  } else {
    p <- NULL
  }
  attr(dge, "p") <- p
  dge
}

diff_test <- function(x, design, contr = NULL, block = NULL){
  if (!is.null(block)){
    dupcor <- e(limma::duplicateCorrelation(x, design, block = block))
    cor <- dupcor$consensus.correlation
    message("## Within-donor correlation:", cor)
  } else {
    cor <- NULL
  }
  fit <- e(limma::lmFit(x, design, block = block, correlation = cor))
  if (!is.null(contr)) {
    fit <- e(limma::contrasts.fit(fit, contrasts = contr))
  }
  # https://liuyujie0136.gitbook.io/sci-tech-notes/bioinformatics/p-value
  fit <- e(limma::eBayes(fit))
  fit
}

extract_tops <- function(x, use = "adj.P.Val", use.cut = 0.05, cut.fc = 0.3){
  res <- e(lapply(1:ncol(x$contrasts),
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

prepare_expr_data <- function(metadata, counts, genes, message = T) 
{
  ## sort and make name for metadata and counts
  checkDup <- function(x) {
    if (any(duplicated(x))) {
      stop("The value in ID column are duplicated.")
    }
  }
  lapply(list(counts[[ 1 ]], genes[[ 1 ]], metadata[[ 1 ]]), checkDup)
  colnames(metadata) %<>% make.names()
  metadata <- dplyr::rename(metadata, sample = 1)
  counts <- dplyr::select(counts, 1, dplyr::all_of(metadata$sample))
  metadata$sample %<>% make.names()
  colnames(counts) %<>% make.names()
  ## sort genes
  data_id <- do.call(data.frame, nl(colnames(genes)[1], list(counts[[1]])))
  genes <- dplyr::distinct(genes, !!rlang::sym(colnames(genes)[1]), .keep_all = T)
  genes <- tbmerge(
    data_id, genes,
    by.x = colnames(data_id)[1], by.y = colnames(genes)[1],
    sort = F, all.x = T
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
  counts <- dplyr::select(counts, -1)
  colnames(counts) <- metadata$sample
  namel(counts, metadata, genes)
}

new_dge <- function(metadata, counts, genes, message = T)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  e(edgeR::DGEList(lst$counts, samples = lst$metadata, genes = lst$genes))
}

new_elist <- function(metadata, counts, genes, message = T)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  elist(list(E = tibble::as_tibble(lst$counts), targets = lst$metadata, genes = lst$genes))
}


