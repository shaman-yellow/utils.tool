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
    params = list(isTcga = F, normed = F),
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
      data <- data.frame(data, check.names = F)
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
  .job_limma(object = data, params = list(metadata = metadata, isTcga = F, normed = T, genes = genes))
}

job_limma <- function(DGEList)
{
  if (!is(DGEList, 'DGEList'))
    stop("is(DGEList, 'DGEList') == F")
  x <- .job_limma(object = DGEList)
  x$.metadata <- tibble::as_tibble(dplyr::relocate(DGEList$samples, sample, group))
  return(x)
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
      object(x) <- norm_genes.dge(object(x), design, vis = norm_vis)
      if (norm_vis) {
        if (length(x@object$targets$sample) < 50) {
          p.norm <- wrap(attr(object(x), "p"), 6, length(x@object$targets$sample) * .6)
        } else {
          p.norm <- wrap(attr(object(x), "p"))
        }
        x@params$p.norm_data <- p.norm@data$data
      } else {
        p.norm <- NULL
      }
      plots <- c(plots, namel(p.norm))
      x@params$normed_data <- object(x)
    } else {
      x$normed_data <- list(
        genes = if (is.null(x$genes)) data.frame(rownames = rownames(object(x))) else x$genes,
        targets = x$metadata,
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
    x$.metadata <- .set_lab(x$.metadata, sig(x), "metadata")
    return(x)
  })

setMethod("step2", signature = c(x = "job_limma"),
  function(x, ..., contrasts = NULL, block = NULL, use = c("adj.P.Val", "P.Value"),
    use.cut = .05, cut.fc = .5,
    label = if (x$isTcga) "gene_name" else "hgnc_symbol", batch = F, HLs = NULL)
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
      if (x$normed) {
        if (!is.null(x$genes)) {
          tops <- lapply(tops,
            function(obj) {
              obj <- map(obj, colnames(obj)[1], x$genes, colnames(x$genes)[1], label, col = label)
              if (!is.null(x$from_scfea)) {
                ## the genes belong to the module
                obj <- map(obj, colnames(obj)[1], x$genes, colnames(x$genes)[1], "gene", col = "gene")
                dic <- nl(x$compounds_annotation$name, x$compounds_annotation$kegg, T)
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
      lab(tops) <- paste(sig(x), "data")
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
    p.sets_intersection <- .set_lab(p.sets_intersection, sig(x), "Difference", "intersection")
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
    p <- .map_boxplot2(data, pvalue)
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
  function(x, y, from, to, names = NULL, use = if (x$isTcga) "gene_name" else "hgnc_symbol",
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
      lst <- .cal_corp.elist(data, anno, use, unique(from), unique(to), names, HLs = HLs, fast = T)
      if (length(unique(from)) > 1 && length(unique(to)) > 1) {
        lst$hp <- .set_lab(wrap(lst$hp), sig(x), theme, "correlation heatmap")
        lst$sig.corp <- .set_lab(lst$sig.corp, sig(x), theme, "significant correlation")
      }
    } else if (mode == "linear") {
      lst <- .cal_corp.elist(data, anno, use, unique(from), unique(to), names, HLs = HLs, fast = F)
    }
    return(lst)
  })

.cal_corp.elist <- function(data, anno, use, from, to, names, HLs = NULL, fast = T)
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
  # heatmap
  if (is.null(names)) {
    corp <- cal_corp(lst[[1]], lst[[2]], "From", "To", trans = T, fast = fast)
  } else {
    corp <- cal_corp(lst[[1]], lst[[2]], names[[1]], names[[2]], trans = T, fast = fast)
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
            c(j = max(zoRange(x[[ "j" ]], 1.3), na.rm = T),
              i = min(x[[ "i" ]], na.rm = T))
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
            label = paste0("Cor = ", round(cor, 2), "\n", "P-value = ", round(pvalue, 5))),
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
  function(x, tnbc = F){
    if (tnbc) {
      if (!is.null(x$isTcga) && identical(x$project, "TCGA-BRCA")) {
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
        stop("!is.null(x$isTcga) && identical(x$project, 'TCGA-BRCA')")
      }
      message("Dim: ", paste0(dim(object(x)), collapse = ", "))
    }
    return(x)
  })


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
  if (any(colnames(metadata) == "sample")) {
    metadata <- dplyr::relocate(metadata, sample)
  } else {
    metadata <- dplyr::rename(metadata, sample = 1)
  }
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


