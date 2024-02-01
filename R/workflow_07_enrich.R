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
    method = "R package `ClusterProfiler` used for gene enrichment analysis"
    ))

job_enrich <- function(ids, annotation, from = "hgnc_symbol", to = "entrezgene_id")
{
  if (!is(ids, "list")) {
    ids <- list(ids = ids)
  }
  if (is.null(names(ids))) {
    stop("is.null(names(ids))")
  }
  if (missing(annotation)) {
    mart <- new_biomart()
    annotation <- filter_biomart(mart, general_attrs(), from, unique(unlist(ids)))
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
  function(x, organism = 'hsa', orgDb = 'org.Hs.eg.db', cl = 4, maxShow = 10){
    step_message("Use clusterProfiler for enrichment.
      "
    )
    cli::cli_alert_info("clusterProfiler::enrichKEGG")
    res.kegg <- multi_enrichKEGG(object(x), organism = organism)
    p.kegg <- vis_enrich.kegg(res.kegg, maxShow = maxShow)
    use.p <- attr(p.kegg, "use.p")
    p.kegg <- lapply(p.kegg, function(x) wrap(x, 8, 4 * (maxShow / 10)))
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
    p.go <- vis_enrich.go(res.go, maxShow = maxShow, use = use.p)
    p.go <- lapply(p.go, function(x) wrap(x))
    p.go <- .set_lab(p.go, sig(x), names(p.go), "GO enrichment")
    res.go <- lapply(res.go,
      function(data) {
        data <- as_tibble(data.table::rbindlist(data, idcol = T))
        data <- dplyr::mutate(data, geneName_list = fun(geneID_list))
        dplyr::relocate(data, ont = .id)
      })
    res.go <- .set_lab(res.go, sig(x), names(res.go), "GO enrichment data")
    x@tables[[ 1 ]] <- namel(res.kegg, res.go)
    x@plots[[ 1 ]] <- namel(p.kegg, p.go)
    x@params$check_go <- check_enrichGO(res.go)
    x$organism <- organism
    return(x)
  })

setMethod("step2", signature = c(x = "job_enrich"),
  function(x, pathways, which.lst = 1, species = x$organism,
    name = paste0("pathview", gs(Sys.time(), " |:", "_")),
    search = "pathview",
    external = F)
  {
    step_message("Use pathview to visualize reults pathway.")
    require(pathview)
    data <- x@tables$step1$res.kegg[[ which.lst ]]
    if (is.null(x$pathview_dir)) {
      x$pathview_dir <- name
    } else {
      name <- x$pathview_dir
    }
    dir.create(name, F)
    setwd(name)
    cli::cli_alert_info("pathview::pathview")
    tryCatch({
      res.pathviews <- lapply(pathways,
        function(pathway) {
          if (!external) {
            data <- dplyr::filter(data, ID == !!pathway)
            pathway <- gs(data$ID, "^[a-zA-Z]*", "")
            genes <- as.character(unlist(data$geneID_list))
          } else {
            pathway <- gs(pathway, "^[a-zA-Z]*", "")
            genes <- x@object$ids
          }
          res.pathview <- try(
            pathview::pathview(gene.data = genes,
              pathway.id = pathway, species = species,
              keys.align = "y", kegg.native = T,
              key.pos = "topright", limit = list(gene = 1, cpd = 1),
              bins = list(gene = 1, cpd = 1),
              na.col = "grey90", discrete = list(gene = T))
          )
          if (inherits(res.pathview, "try-error")) {
            try(dev.off(), silent = T)
          }
          return(res.pathview)
        })
    }, finally = {setwd("../")})
    x@tables[[ 2 ]] <- namel(res.pathviews)
    figs <- list.files(name, search, full.names = T)
    p.pathviews <- lapply(figs, function(x) .file_fig(x))
    names(p.pathviews) <- get_realname(figs)
    x@plots[[ 2 ]] <- namel(p.pathviews)
    return(x)
  })

setGeneric("asjob_enrich", 
  function(x, ...) standardGeneric("asjob_enrich"))

setMethod("asjob_enrich", signature = c(x = "job_seurat"),
  function(x, exclude.pattern = "macroph", exclude.use = "scsa_cell",
    ignore.case = T, marker.list = x@params$contrasts, geneType = "hgnc_symbol")
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
      res <- sapply(onts, simplify = F,
        function(ont) {
          res.go <- try(clusterProfiler::enrichGO(ids, orgDb, ont = ont), T)
          if (inherits(res.go, "try-error")) {
            return("try-error of enrichment")
          }
          res.res <- try(res.go@result, T)
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
      p <- ggplot(data) +
        geom_point(aes(x = reorder(Description, GeneRatio),
            y = GeneRatio, size = Count, fill = !!rlang::sym(use)),
          shape = 21, stroke = 0, color = "transparent") +
        scale_fill_gradient(high = "yellow", low = "red") +
        scale_size(range = c(4, 6)) +
        labs(x = "", y = "Gene Ratio") +
        guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
        coord_flip() +
        ylim(zoRange(data$GeneRatio, 1.3)) +
        theme_minimal()
      p
    })
  if (length(lst) == 1) {
    attr(res, "use.p") <- use.p
  }
  res
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
    data <- data.table::rbindlist(data, idcol = T)
    if (!nrow(data)) {
      return()
    }
    data <- dplyr::mutate(
      data, GeneRatio = as_double.ratioCh(GeneRatio),
      stringr::str_wrap(Description, width = 30)
    )
    p <- ggplot(data) +
      geom_point(aes(x = reorder(Description, GeneRatio),
          y = GeneRatio, size = Count, fill = !!rlang::sym(use)),
        shape = 21, stroke = 0, color = "transparent") +
      scale_fill_gradient(high = "yellow", low = "red") +
      scale_size(range = c(4, 6)) +
      guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
      coord_flip() +
      facet_grid(.id ~ ., scales = "free") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey70")) +
      geom_blank()
    p
  }
  res <- lapply(lst,
    function(x) {
      try(fun(x), silent = T)
    })
  res
}
