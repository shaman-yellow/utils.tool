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
    info = c("..."),
    cite = "[@LimmaPowersDiRitchi2015; @EdgerDifferenChen]",
    method = "Limma and edgeR used for differential expression analysis"
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
  function(x, group = x@object$samples$group, batch = x@object$samples$batch,
    design = if (is.null(batch)) mx(~ 0 + group) else mx(~ 0 + group + batch),
    min.count = 10, no.filter = F, no.norm = F, norm_vis = F)
  {
    step_message("Preprocess expression data.
      "
    )
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
    }
    if (F) {
      pca <- pca_data.long(as_data_long(object(x)))
      p.pca <- plot_andata(pca)
    }
    x@plots[[ 1 ]] <- plots
    x@params$p.norm_data <- p.norm@data$data
    x@params$group <- group
    x@params$design <- design
    x@params$normed_data <- object(x)
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
    if (!is.null(contr)) {
      tops <- extract_tops(object(x), use = use, use.cut = use.cut, cut.fc = cut.fc)
      tops <- .set_lab(tops, sig(x), paste("data", gs(names(tops), "-", "vs")), "DEGs")
      lab(tops) <- paste(sig(x), "data", "DEGs")
    } else {
      tops <- NULL
    }
    p.valcano <- lapply(tops, plot_valcano, label = label, use = use, fc = cut.fc)
    p.valcano <- lapply(p.valcano, function(p) wrap(p, 5, 4))
    p.valcano <- .set_lab(p.valcano, sig(x), gs(names(p.valcano), "-", "vs"), "DEGs")
    lab(p.valcano) <- paste(sig(x), "volcano plot", "DEGs")
    x@tables[[ 2 ]] <- namel(tops)
    x@plots[[ 2 ]] <- namel(p.valcano)
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
  data <- dplyr::select(top_table, !!rlang::sym(label), logFC, P.Value, adj.P.Val)
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
  function(x, y, from, to, names = NULL, use = if (x$isTcga) "gene_name" else "hgnc_symbol")
  {
    data <- as_tibble(x@params$normed_data$E)
    anno <- x@params$normed_data$genes
    .cal_corp.elist(data, anno, use, from, to, names)
  })

.cal_corp.elist <- function(data, anno, use, from, to, names)
{
  if (is.null(data$rownames)) {
    data <- as_tibble(data)
  }
  data$rownames <- anno[[ use ]]
  colnames(data)[1] <- use
  data <- dplyr::mutate(data, symbol = gname(!!rlang::sym(use)))
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
