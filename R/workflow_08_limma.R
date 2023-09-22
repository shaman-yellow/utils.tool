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
    info = c("...")
    ))

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
  function(x, group = x@object$samples$group, design = mx(~ 0 + group), min.count = 10){
    step_message("Preprocess expression data.
      "
    )
    object(x) <- filter_low.dge(object(x), group, min.count = min.count)
    p.filter <- wrap(attr(object(x), "p"), 8, 3)
    object(x) <- norm_genes.dge(object(x), design)
    if (length(x@object$targets$sample) < 50) {
      p.norm <- wrap(attr(object(x), "p"), 6, length(x@object$targets$sample) * .6)
    } else {
      p.norm <- wrap(attr(object(x), "p"))
    }
    if (F) {
      pca <- pca_data.long(as_data_long(object(x)))
      p.pca <- plot_andata(pca)
    }
    x@plots[[ 1 ]] <- namel(p.filter, p.norm)
    x@params$p.norm_data <- p.norm@data$data
    x@params$group <- group
    x@params$design <- design
    x@params$normed_data <- object(x)
    return(x)
  })

setMethod("step2", signature = c(x = "job_limma"),
  function(x, ..., contrasts = NULL, block = NULL, use = c("adj.P.Val", "P.Value"),
    use.cut = .05, cut.fc = .3, label = "hgnc_symbol")
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
    object(x) <- diff_test(object(x), x@params$design, contr, block)
    if (!is.null(contr)) {
      tops <- extract_tops(object(x), use = use, use.cut = use.cut, cut.fc = cut.fc)
    }
    p.valcano <- lapply(tops, plot_valcano, label = label, use = use, fc = cut.fc)
    p.valcano <- lapply(p.valcano, function(p) wrap(p, 5, 4))
    x@tables[[ 2 ]] <- namel(tops)
    x@plots[[ 2 ]] <- namel(p.valcano)
    return(x)
  })

setMethod("clear", signature = c(x = "job_limma"),
  function(x, save = T, suffix = NULL){
    if (save)
      saveRDS(x, paste0(substitute(x, parent.frame(1)), x@step, suffix, ".rds"))
    object(x) <- NULL
    x@params$normed_data <- NULL
    return(x)
  })

setMethod("map", signature = c(x = "job_limma"),
  function(x, ref, ref.use = "hgnc_symbol", group = NULL, pvalue = T){
    object <- x@params$normed_data
    rownames(object) <- object$genes[[ ref.use ]]
    object <- object[rownames(object) %in% ref, ]
    if (!is.null(group)) {
      object <- object[, object$targets$group %in% group]
    }
    data <- tibble::as_tibble(t(object$E))
    data$group <- object$targets$group
    data <- tidyr::gather(data, var, value, -group)
    p <- .map_boxplot(data, pvalue)
    p
  })

plot_valcano <- function(top_table, label = "hgnc_symbol", use = "adj.P.Val", fc = .3) {
  data <- select(top_table, !!rlang::sym(label), logFC, P.Value, adj.P.Val)
  data <- mutate(data,
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
  function(x, filter_genes = lm@tables$step2$tops[[1]][[1]]){
    step_message("Use `x@params$normed_data` converted as job_wgcna.")
    if (is.null(object <- x@params$normed_data))
      stop("is.null(x@params$normed_data)")
    log_counts <- as_tibble(object$E)
    gene_annotation <- select(as_tibble(object$genes), -1)
    log_counts[[1]] <- gene_annotation[[1]]
    if (!is.null(filter_genes)) {
      log_counts <- filter(log_counts, rownames %in% filter_genes)
    }
    job_wgcna(select(object$targets, sample, group),
      log_counts, gene_annotation)
  })

