# ==========================================================================
# workflow of enrich
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_enrich <- setClass("job_enrich", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@ClusterprofilerWuTi2021]",
    method = "R package `ClusterProfiler` used for gene enrichment analysis",
    tag = "enrich:clusterProfiler",
    analysis = "富集分析"
    ))

setGeneric("asjob_enrich", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_enrich"))

setMethod("asjob_enrich", signature = c(x = "feature"),
  function(x, ...){
    x <- resolve_feature_snapAdd_onExit("x", x)
    x <- job_enrich(x)
    return(x)
  })

job_enrich <- function(ids, annotation, from = "hgnc_symbol", to = "entrezgene_id")
{
  if (!is(ids, "list")) {
    ids <- list(ids = rm.no(ids))
  }
  if (is.null(names(ids))) {
    stop("is.null(names(ids))")
  }
  if (missing(annotation)) {
    mart <- new_biomart()
    annotation <- filter_biomart(
      mart, general_attrs(), from, unique(unlist(ids))
    )
  }
  maps <- lapply(ids,
    function(id) {
      res <- unique(filter(annotation, !!rlang::sym(from) %in% !!id)[[ to ]])
      res[!is.na(res)]
    })
  en <- .job_enrich(object = maps)
  en@params$raw <- ids
  en@params$annotation <- annotation
  en$from <- from
  en
}

setMethod("step0", signature = c(x = "job_enrich"),
  function(x){
    step_message("Prepare your data with function `job_enrich`.
      "
    )
  })

# Biological Process, Molecular Function, and Cellular Component groups
setMethod("step1", signature = c(x = "job_enrich"),
  function(x, organism = c("hsa", "mmu", "rno"),
    orgDb = switch(organism,
      hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db", rno = "org.Rn.eg.db"),
    cl = 4, maxShow.kegg = 10,
    maxShow.go = 10, use = c("p.adjust", "pvalue"))
  {
    step_message("Use clusterProfiler for enrichment.")
    cli::cli_alert_info("clusterProfiler::enrichKEGG")
    organism <- match.arg(organism)
    orgDb <- match.arg(orgDb)
    use <- match.arg(use)
    res.kegg <- multi_enrichKEGG(object(x), organism = organism)
    p.kegg <- vis_enrich.kegg(res.kegg, maxShow = maxShow.kegg, use = use)
    use.p <- attr(p.kegg, "use.p")
    p.kegg <- .set_lab(p.kegg, sig(x), names(p.kegg), "KEGG enrichment")
    fun <- function(sets) {
      lapply(sets,
        function(set) {
          from_ids <- x@params$annotation$entrezgene_id
          to_names <- x@params$annotation[[ x$from ]]
          to_names[ match(set, from_ids) ]
        })
    }
    res.kegg <- lapply(res.kegg, mutate, geneName_list = fun(geneID_list))
    res.kegg <- .set_lab(res.kegg, sig(x), names(res.kegg), "KEGG enrichment data")
    cli::cli_alert_info("clusterProfiler::enrichGO")
    res.go <- multi_enrichGO(object(x), orgDb = orgDb, cl = cl)
    p.go <- vis_enrich.go(res.go, maxShow = maxShow.go, use = use.p)
    p.go <- .set_lab(p.go, sig(x), names(p.go), "GO enrichment")
    res.go <- lapply(res.go,
      function(data) {
        if (all(vapply(data, is.data.frame, logical(1)))) {
          data <- as_tibble(data.table::rbindlist(data, idcol = TRUE))
          data <- dplyr::mutate(data, geneName_list = fun(geneID_list))
          dplyr::relocate(data, ont = .id)
        }
      })
    res.go <- .set_lab(res.go, sig(x), names(res.go), "GO enrichment data")
    x@tables[[ 1 ]] <- namel(res.kegg, res.go)
    x@plots[[ 1 ]] <- namel(p.kegg, p.go)
    x@params$check_go <- check_enrichGO(res.go)
    x$organism <- organism
    methodAdd_onExit("x", "以 ClusterProfiler R 包 ({packageVersion('clusterProfiler')}) {cite_show('ClusterprofilerWuTi2021')}进行 KEGG 和 GO 富集分析。以 {use} 表示显著水平。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_enrich"),
  function(x, pathways, which.lst = 1, species = x$organism,
    name = paste0("pathview", gs(Sys.time(), " |:", "_")),
    search = "pathview",
    external = FALSE, gene.level = NULL, gene.level.name = "hgnc_symbol")
  {
    step_message("Use pathview to visualize reults pathway.")
    require(pathview)
    data <- x@tables$step1$res.kegg[[ which.lst ]]
    if (is.null(x$pathview_dir)) {
      x$pathview_dir <- name
    } else {
      name <- x$pathview_dir
    }
    if (!is.null(gene.level)) {
      if (is(gene.level, "data.frame")) {
        message("Use first (symbol) and second (logFC) columns of `gene.level`.")
        gene.level <- nl(gene.level[[1]], gene.level[[2]], FALSE)
      } else if (is.numeric(gene.level)) (
        if (is.null(names(gene.level))) {
          stop("is.null(names(gene.level))")
        }
      )
      message("Note that only hgnc_symbol support for this feature: `gene.level`")
      names(gene.level) <- x$annotation$entrezgene_id[ match(names(gene.level),
        x$annotation[[ gene.level.name ]]) ]
      snap <- "通路图中的基因的映射颜色表示基因显著富集，并在数据集中有上调或下调变化趋势。"
    } else {
      snap <- "通路图中的基因的映射颜色表示是否显著富集。"
    }
    dir.create(name, FALSE)
    setwd(name)
    cli::cli_alert_info("pathview::pathview")
    tryCatch({
      res.pathviews <- sapply(pathways, simplify = FALSE,
        function(pathway) {
          if (!external) {
            data <- dplyr::filter(data, ID == !!pathway)
            pathway <- gs(data$ID, "^[a-zA-Z]*", "")
            genes <- as.character(unlist(data$geneID_list))
          } else {
            pathway <- gs(pathway, "^[a-zA-Z]*", "")
            genes <- x@object$ids
          }
          if (!is.null(gene.level)) {
            genes <- gene.level[ match(genes, names(gene.level)) ]
            discrete <- FALSE
            bins <- 10
          } else {
            discrete <- TRUE
            bins <- 1
          }
          res.pathview <- try(
            pathview::pathview(
              gene.data = genes,
              pathway.id = pathway, species = species,
              keys.align = "y", kegg.native = TRUE, same.layer = FALSE,
              key.pos = "topright", bins = list(gene = bins),
              na.col = "grey90", discrete = list(gene = discrete)
            )
          )
          if (inherits(res.pathview, "try-error")) {
            try(dev.off(), silent = TRUE)
          }
          return(res.pathview)
        })
    }, finally = {setwd("../")})
    x@tables[[ 2 ]] <- namel(res.pathviews)
    p.pathviews <- .pathview_search(
      name, search, x, res.pathviews, snap = snap
    )
    feature(x) <- lapply(p.pathviews, function(x) x$genes)
    x@plots[[ 2 ]] <- namel(p.pathviews)
    methodAdd_onExit("x", "以 `pathview` R 包 ({packageVersion('pathview')}) 对选择的 KEGG 通路可视化。")
    snapAdd_onExit("x", "以 `pathview` 探究基因集在通路 {bind(pathways)} 中的上下游关系。")
    .add_internal_job(.job(method = "R package `pathview` used for KEGG pathways visualization", cite = "[@PathviewAnRLuoW2013]"))
    return(x)
  })

setMethod("feature", signature = c(x = "job_enrich"),
  function(x, ref = 1L){
    if (identical(ref, "all")) {
      feas <- x$.feature
      as_feature(
        feas, x, analysis = "通路 ({bind(names(x$.feature))}) 中的富集基因"
      )
    } else {
      feas <- x$.feature[[ ref ]]
      as_feature(
        feas, x, analysis = "通路 ({names(x$.feature)[ref]}) 中的富集基因"
      )
    }
  })

setMethod("res", signature = c(x = "job_enrich", ref = "character"),
  function(x, ref = c("id", "des", "cate", "sub"), which = 1, key = 1, from = c("kegg", "go"))
  {
    type <- match.arg(ref)
    type <- switch(type, id = "ID", des = "Description", cate = "category", sub = "subcategory")
    from <- match.arg(from)
    data <- x@tables$step1[[ paste0("res.", from) ]][[ key ]]
    data[[ type ]][ which ]
  })

setMethod("asjob_enrich", signature = c(x = "job_seurat"),
  function(x, exclude.pattern = NULL, exclude.use = NULL,
    ignore.case = TRUE, marker.list = x@params$contrasts, geneType = "hgnc_symbol")
  {
    if (is.null(marker.list)) {
      if (x@step < 5) {
        stop("x@step < 5")
      }
      data <- x@tables$step5[[ "all_markers" ]]
      if (!is.null(exclude.pattern)) {
        exclude.cluster <- dplyr::filter(object(x)@meta.data,
          grepl(!!exclude.pattern, !!rlang::sym(exclude.use), ignore.case))$seurat_clusters
        exclude.cluster <- unique(exclude.cluster)
        message("Exclude clasters:\n  ", paste0(exclude.cluster, collapse = ", "))
        data <- dplyr::filter(data, !cluster %in% exclude.cluster)
      }
    } else {
      data <- marker.list
    }
    data <- dplyr::mutate(data, gene = gs(gene, "\\.[0-9]*$", ""))
    mart <- new_biomart()
    anno <- filter_biomart(mart, general_attrs(), geneType, unique(data$gene))
    data <- dplyr::filter(data, gene %in% anno[[ !!geneType ]])
    if (is.null(data$contrast)) {
      ids <- split(data$gene, data$cluster)
    } else {
      ids <- split(data$gene, data$contrast)
    }
    ids <- lst_clear0(ids)
    job_enrich(ids, anno)
  })

setMethod("focus", signature = c(x = "job_enrich"),
  function(x, symbols, data = x@tables$step1$res.kegg[[1]])
  {
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    isThat <- vapply(data$geneName_list, FUN.VALUE = logical(1),
      function(genes) {
        any(genes %in% symbols)
      })
    data <- dplyr::filter(data, !!isThat)
    data
  })

multi_enrichKEGG <- function(lst.entrez_id, organism = 'hsa')
{
  res <- pbapply::pblapply(lst.entrez_id,
    function(ids) {
      res.kegg <- clusterProfiler::enrichKEGG(ids, organism = organism)
      res.path <- tibble::as_tibble(res.kegg@result)
      res.path <- dplyr::mutate(res.path, geneID_list = lapply(strsplit(geneID, "/"), as.integer))
      res.path
    })
  res
}

multi_enrichGO <- function(lst.entrez_id, orgDb = 'org.Hs.eg.db', cl = NULL)
{
  res <- pbapply::pblapply(lst.entrez_id, cl = cl,
    function(ids) {
      onts <- c("BP", "CC", "MF")
      res <- sapply(onts, simplify = FALSE,
        function(ont) {
          res.go <- try(clusterProfiler::enrichGO(ids, orgDb, ont = ont), TRUE)
          if (inherits(res.go, "try-error")) {
            return("try-error of enrichment")
          }
          res.res <- try(res.go@result, TRUE)
          if (inherits(res.res, "try-error")) {
            value <- "try-error of enrichment"
            attr(value, "data") <- res.res
            return(value)
          }
          res.path <- tibble::as_tibble(res.res)
          res.path <- dplyr::mutate(res.path, geneID_list = lapply(strsplit(geneID, "/"), as.integer))
          res.path
        })
    })
}

check_enrichGO <- function(res.go) {
  isthat <- lapply(res.go,
    function(res) {
      !vapply(res, FUN.VALUE = logical(1), is.character)
    })
  isthat
}

vis_enrich.kegg <- function(lst, cutoff = .1, maxShow = 10,
  use = c("p.adjust", "pvalue"), least = 3L)
{
  use <- match.arg(use)
  use.p <- use
  res <- lapply(lst,
    function(data) {
      data <- dplyr::filter(raw <- data, !!rlang::sym(use) < !!cutoff)
      if (!nrow(data) | nrow(data) < least) {
        message("\n", "Too few of the results (", nrow(data), ")")
        if (use == "p.adjust") {
          message("Switch to use `pvalue`.")
          use <- "pvalue"
          use.p <<- use
          data <- dplyr::filter(raw, !!rlang::sym(use) < !!cutoff)
        }
      }
      data <- dplyr::arrange(data, !!rlang::sym(use))
      data <- head(data, n = maxShow)
      data <- dplyr::mutate(data, GeneRatio = as_double.ratioCh(GeneRatio))
      p <- .plot_kegg(data, use)
      p <- wrap(p, 7, 4 * (maxShow / 10))
      p <- setLegend(p, "KEGG 富集图展示了以 {use} 排序，前 {maxShow} 的富集通路。")
      p
    })
  attr(res, "use.p") <- use.p
  res
}

# for external use
plot_kegg <- function(data, cutoff = .1, maxShow = 10,
  use = c("p.adjust", "pvalue"), pattern = NULL, ...)
{
  use <- match.arg(use)
  data <- .format_enrich(data, use, cutoff, maxShow, pattern)
  p <- .plot_kegg(data, use, ...)
  p <- wrap(p, 8, 2 + nrow(data) * .2)
  p <- .set_lab(p, "KEGG-enrichment")
  p
}

.plot_kegg <- function(data, use, ratio = "GeneRatio", count = "Count") {
  p <- ggplot(data) +
    geom_point(aes(x = reorder(Description, !!rlang::sym(ratio)),
        y = !!rlang::sym(ratio), size = !!rlang::sym(count), fill = !!rlang::sym(use)),
      shape = 21, stroke = 0, color = "transparent") +
    scale_fill_gradient(high = "grey90", low = "darkred") +
    scale_size(range = c(4, 6)) +
    labs(x = "", y = "Hits Ratio") +
    guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
    coord_flip() +
    ylim(zoRange(data[[ ratio ]], 1.3)) +
    rstyle("theme")
  p
}

vis_enrich.go <- function(lst, cutoff = .1, maxShow = 10,
  use = c("p.adjust", "pvalue"), least = 3L)
{
  use <- match.arg(use)
  fun <- function(data) {
    data <- lapply(data,
      function(data) {
        if (is.character(data))
          return()
        data <- dplyr::filter(data, !!rlang::sym(use) < cutoff)
        data <- dplyr::arrange(data, !!rlang::sym(use))
        data <- head(data, n = maxShow)
        data
      })
    data <- data.table::rbindlist(data, idcol = TRUE)
    if (!nrow(data)) {
      return()
    }
    data <- dplyr::mutate(
      data, GeneRatio = as_double.ratioCh(GeneRatio),
      stringr::str_wrap(Description, width = 30)
    )
    p <- .plot_go_use_ggplot(data, use)
    p <- wrap(p, 7.5, 1 + nrow(data) * .2)
    p <- setLegend(p, "GO 富集图展示了基因集在 GO 的 BP (Biological Process), MF (Molecular Function), CC (Cellular Component) 组中的富集结果 (以 {use} 排序，各自展示前 {maxShow} 的富集通路) 。")
    p
  }
  res <- lapply(lst,
    function(x) {
      try(fun(x), silent = TRUE)
    })
  res
}

# this function is for external use. 
plot_go <- function(data, cutoff = .1, maxShow = 10,
  use = c("p.adjust", "pvalue"), facet = ".id", pattern = NULL)
{
  use <- match.arg(use)
  data <- lapply(split(data, data[[ facet ]]),
    function(data) {
      .format_enrich(data, use, cutoff, maxShow, pattern)
    })
  data <- frbind(data)
  p <- .plot_go_use_ggplot(data, use, facet)
  p <- wrap(p, 7.5, 1 + nrow(data) * .2)
  p <- .set_lab(p, "GO-enrichment")
  p <- setLegend(p, "GO 富集图展示了基因集在 GO 的 BP (Biological Process), MF (Molecular Function), CC (Cellular Component) 组中的富集结果 (以 {use} 排序，各自展示前 {maxShow} 的富集通路) 。")
  p
}

.format_enrich <- function(data, use, cutoff, maxShow, pattern = NULL)
{
  data <- dplyr::filter(data, !!rlang::sym(use) < cutoff)
  data <- dplyr::arrange(data, !!rlang::sym(use))
  less <- head(data, n = maxShow)
  if (!is.null(pattern)) {
    extra <- dplyr::filter(data, grpl(Description, pattern, TRUE))
    if (nrow(extra)) {
      less <- dplyr::bind_rows(less, extra)
      less <- dplyr::distinct(less)
    }
  }
  if (!is.null(data$GeneRatio)) {
    less <- dplyr::mutate(less, GeneRatio = as_double.ratioCh(GeneRatio))
    less <- dplyr::arrange(less, GeneRatio)
  }
  less
}

.plot_go_use_ggplot <- function(data, use, facet = ".id")
{
  p <- ggplot(data) +
    geom_point(aes(x = reorder(Description, GeneRatio),
        y = GeneRatio, size = Count, fill = !!rlang::sym(use)),
      shape = 21, stroke = 0, color = "transparent") +
    scale_fill_gradient(high = "grey90", low = "darkred") +
    scale_size(range = c(4, 6)) +
    guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
    coord_flip() +
    ggplot2::facet_grid(rows = ggplot2::vars(!!rlang::sym(facet)), scales = "free") +
    rstyle("theme") +
    theme(axis.title.y = element_blank()) +
    geom_blank()
  p
}

setMethod("map", signature = c(x = "job_enrich", ref = "job_enrich"),
  function(x, ref, use = c("kegg", "go"), key = 1, cutoff = .05, use.cutoff = c("p.adjust", "pvalue"))
  {
    message("Find intersection pathways across two 'job_enrich'")
    use <- match.arg(use)
    use.cutoff <- match.arg(use.cutoff)
    fun_extract <- function(x) {
      dplyr::filter(x@tables$step1[[ paste0("res.", use) ]][[ key ]],
        !!rlang::sym(use.cutoff) < cutoff)
    }
    lst <- lapply(list(x, ref), fun_extract)
    data <- dplyr::filter(lst[[1]], ID %in% !!lst[[2]]$ID)
    data <- .set_lab(data, sig(x), "pathways intersection")
    x$intersect_paths <- data
    return(x)
  })

get_genes.keggPath <- function(name) {
  if (!is(name, "character")) {
    stop("is(name, 'character')")
  }
  lst <- e(KEGGREST::keggGet(name))
  x <- strx(lst[[1]]$GENE, "^[A-Za-z][^;]+")
  x[ !is.na(x) ]
}

as_double.ratioCh <- function(ch) {
  values <- stringr::str_extract_all(ch, "[0-9]{1,}")
  vapply(values, FUN.VALUE = double(1),
    function(values) {
      values <- as.double(values)
      values[ 1 ] / values[ 2 ]
    })
}
