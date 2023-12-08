# ==========================================================================
# workflow of gsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gsea <- setClass("job_gsea", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://www.gsea-msigdb.org/gsea/downloads.jsp"),
    method = "ClusterProfiler used for GSEA enrichment",
    cite = "[@ClusterprofilerWuTi2021]"
    ))

setGeneric("asjob_gsea", 
  function(x, ...) standardGeneric("asjob_gsea"))

setMethod("asjob_gsea", signature = c(x = "job_limma"),
  function(x, key = 1L, annotation = x@params$normed_data$genes,
    filter = NULL, from = colnames(annotation)[ grp(colnames(annotation), "_symbol")[1] ])
  {
    data <- x@tables$step2$tops[[ key ]]
    if (!is.null(filter)) {
      data <- dplyr::filter(data, !!rlang::sym(from) %in% dplyr::all_of(filter))
    }
    data <- dplyr::rename(data, symbol = !!rlang::sym(from))
    annotation <- dplyr::rename(annotation, symbol = !!rlang::sym(from))
    x$from <- from
    job_gsea(data, annotation, use = "symbol")
  })

setMethod("asjob_gsea", signature = c(x = "job_seurat"),
  function(x, contrast.pattern = NULL, marker.list = x@params$contrasts)
  {
    if (!is.null(contrast.pattern)) {
      topTable <- dplyr::filter(marker.list, grpl(contrast, !!contrast.pattern))
    } else {
      topTable <- marker.list
    }
    x <- job_gsea(dplyr::relocate(topTable, symbol = gene, logFC = avg_log2FC), use = "symbol")
    sig(x) <- gs(contrast.pattern, "_", " ")
    x
  })

job_gsea <- function(topTable, annotation, use)
{
  .check_columns(topTable, c("logFC"))
  rename <- T
  if (missing(use)) {
    if (any(colnames(topTable) == "symbol")) {
      message("Use column `symbol` as gene input.")
      rename <- F
      use <- "symbol"
    } else if (length(whichs <- grp(colnames(topTable), "_symbol")) >= 1) {
      use <- colnames(topTable)[ whichs[1] ]
      message("Use column `", use, "` as gene input.")
    }
  } else if (use == "symbol") {
    rename <- F
  }
  if (rename) {
    topTable <- dplyr::rename(topTable, symbol = !!rlang::sym(use))
    if (!missing(annotation))
      annotation <- dplyr::rename(annotation, symbol = !!rlang::sym(use))
    use <- "symbol"
  }
  topTable <- dplyr::select(topTable, symbol, logFC)
  topTable <- dplyr::filter(topTable, !is.na(symbol) & symbol != "")
  topTable <- dplyr::arrange(topTable, dplyr::desc(logFC))
  topTable <- dplyr::distinct(topTable, symbol, .keep_all = T)
  used.symbol <- nl(topTable$symbol, topTable$logFC, F)
  if (missing(annotation)) {
    message("Missing `annotation`, use 'hgnc_symbol' to get annotation.")
    mart <- new_biomart()
    annotation <- filter_biomart(mart,
      c("hgnc_symbol", "entrezgene_id"), "hgnc_symbol", topTable$symbol)
    annotation <- dplyr::rename(annotation, !!!nl(use, "hgnc_symbol"))
  }
  topTable <- map(topTable, "symbol", annotation, use, "entrezgene_id")
  used.entrezgene_id <- nl(topTable$entrezgene_id, topTable$logFC, F)
  object <- list(symbol = used.symbol, entrezgene_id = used.entrezgene_id)
  x <- .job_gsea(object = object)
  x$annotation <- annotation
  return(x)
}

setMethod("step0", signature = c(x = "job_gsea"),
  function(x){
    step_message("Prepare your data with function `job_gsea`.")
  })

setMethod("step1", signature = c(x = "job_gsea"),
  function(x, OrgDb = org.Hs.eg.db::org.Hs.eg.db, org = "hsa"){
    step_message("GSEA enrichment.")
    isNa <- is.na(names(object(x)$entrezgene_id))
    object(x)$entrezgene_id <- object(x)$entrezgene_id[!isNa]
    object(x)$symbol <- object(x)$symbol[!isNa]
    ## kegg
    res.kegg <- e(clusterProfiler::gseKEGG(object(x)$entrezgene_id, organism = org))
    fun <- function(sets) {
      lapply(sets,
        function(set) {
          from_ids <- x$annotation$entrezgene_id
          to_names <- x$annotation$symbol
          to_names[ match(set, from_ids) ]
        })
    }
    table_kegg <- try(dplyr::mutate(res.kegg@result,
        geneID_list = strsplit(core_enrichment, "/"),
        geneName_list = fun(geneID_list),
        GeneRatio = stringr::str_extract(leading_edge, "[0-9]+"),
        Count = lengths(geneName_list)
        ), T)
    if (!inherits(table_kegg, "try-error")) {
      table_kegg <- dplyr::as_tibble(table_kegg)
      p.kegg <- e(enrichplot::dotplot(res.kegg))
    } else {
      table_kegg <- NULL
      p.kegg <- NULL
    }
    ## go
    if (is.null(x$res.go)) {
      res.go <- e(clusterProfiler::gseGO(object(x)$symbol, ont = "ALL", OrgDb = OrgDb,
          keyType = "SYMBOL"))
    } else {
      res.go <- x$res.go
    }
    if (nrow(res.go@result) == 0) {
      message("No enrichment available for GO.")
      table_go <- NULL
      p.go <- NULL
    } else {
      table_go <- dplyr::mutate(res.go@result,
        geneName_list = strsplit(core_enrichment, "/"),
        GeneRatio = as.double(stringr::str_extract(leading_edge, "[0-9]+")),
        Count = lengths(geneName_list)
      )
      table_go <- dplyr::arrange(table_go, p.adjust)
      table_go <- as_tibble(table_go)
      data <- dplyr::mutate(split_lapply_rbind(table_go, ~ ONTOLOGY, head, n = 10),
        Description = stringr::str_trunc(Description, 50)
      )
      p.go <- ggplot(data) +
        geom_point(aes(x = reorder(Description, GeneRatio),
            y = GeneRatio, size = Count, fill = p.adjust),
          shape = 21, stroke = 0, color = "transparent") +
        scale_fill_gradient(high = "yellow", low = "red") +
        scale_size(range = c(4, 6)) +
        guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
        coord_flip() +
        facet_grid(ONTOLOGY ~ ., scales = "free") +
        theme_minimal() +
        theme(axis.title.y = element_blank(),
          strip.background = element_rect(fill = "grey90", color = "grey70")) +
        geom_blank()
    }
    x@params$res.go <- res.go
    x@params$res.kegg <- res.kegg
    x@tables[[ 1 ]] <- namel(table_go, table_kegg)
    p.go <- .set_lab(wrap(p.go), sig(x), "GO", "enrichment")
    p.kegg <- .set_lab(wrap(p.kegg), sig(x), "KEGG", "enrichment")
    x@plots[[ 1 ]] <- namel(p.go, p.kegg)
    x$org <- org
    return(x)
  })

setMethod("step2", signature = c(x = "job_gsea"),
  function(x, key, highlight = key[1], use = "res.kegg"){
    step_message("GSEA visualization for specific pathway")
    obj <- x@params[[ use ]]
    p.code <- wrap(e(enrichplot::gseaplot2(obj, key, pvalue_table = F)), 7.5, 6)
    p.code <- .set_lab(p.code, sig(x), "GSEA plot of the pathways")
    if (!is.null(highlight)) {
      p.highlight <- plot_highlight_enrich(x@tables$step1$table_kegg, highlight, object(x)$symbol)
      p.highlight <- .set_lab(p.highlight, sig(x), "KEGG enrichment with enriched genes")
    } else {
      p.highlight <- NULL
    }
    x@plots[[ 2 ]] <- namel(p.code, p.highlight)
    return(x)
  })

setMethod("step3", signature = c(x = "job_gsea"),
  function(x, pathways, species = x$org,
    name = paste0("pathview", gs(Sys.time(), " |:", "_")),
    search = "pathview")
  {
    step_message("Use pathview to visualize reults pathway.")
    require(pathview)
    data <- x@tables$step1$table_kegg
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
          data <- dplyr::filter(data, ID == !!pathway)
          pathway <- gs(data$ID, "^[a-zA-Z]*", "")
          genes <- as.character(unlist(data$geneID_list))
          genes.all <- x@object$entrezgene_id
          genes <- genes.all[ match(genes, names(genes.all)) ]
          res.pathview <- try(
            pathview::pathview(gene.data = genes,
              pathway.id = pathway, species = species,
              keys.align = "y", kegg.native = T,
              key.pos = "topright", na.col = "grey90")
          )
          if (inherits(res.pathview, "try-error")) {
            try(dev.off(), silent = T)
          }
          return(res.pathview)
        })
    }, finally = {setwd("../")})
    x@tables[[ 3 ]] <- namel(res.pathviews)
    figs <- list.files(name, search, full.names = T)
    p.pathviews <- lapply(figs, function(x) .file_fig(x))
    names(p.pathviews) <- get_realname(figs)
    x@plots[[ 3 ]] <- namel(p.pathviews)
    return(x)
  })

setMethod("step4", signature = c(x = "job_gsea"),
  function(x, db, cutoff = .05, map = NULL, pvalue = F){
    step_message("Custom database for GSEA enrichment.")
    ## general analysis
    insDb <- lapply(split(db, ~ term),
      function(data) intersect(data$symbol, names(object(x)$symbol)))
    p.pie_insDb <- new_pie(rep(names(insDb), lengths(insDb)))
    table_insDb <- dplyr::filter(db, symbol %in% names(object(x)$symbol))
    ## enrichment
    res.gsea <- e(clusterProfiler::GSEA(object(x)$symbol, TERM2GENE = db,
        pvalueCutoff = cutoff))
    table_gsea <- dplyr::as_tibble(res.gsea@result)
    if (!is.null(map)) {
      alls <- table_gsea$ID
      map <- alls[ grepl(map, alls) ]
      p.code <- wrap(e(enrichplot::gseaplot2(res.gsea, map, pvalue_table = pvalue)), 7.5, 6)
    } else {
      p.code <- NULL
    }
    table_gsea <- dplyr::mutate(table_gsea,
      geneName_list = strsplit(core_enrichment, "/"),
      Count = lengths(geneName_list),
      GeneRatio = round(as.double(stringr::str_extract(leading_edge, "[0-9]+")) / 100, 2)
    )
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x@tables[[ 4 ]] <- namel(table_gsea, table_insDb)
    p.code <- .set_lab(p.code, sig(x), "GSEA plot of pathway")
    x@plots[[ 4 ]] <- namel(p.code, p.pie_insDb)
    return(x)
  })

setMethod("filter", signature = c(DF_object = "job_gsea"),
  function(DF_object, ref, use = c("entrezgene_id", "symbol")){
    use <- match.arg(use)
    isThat <- names(object(DF_object)[[ use ]]) %in% ref
    object(DF_object) <- lapply(object(DF_object),
      function(x) {
        x[ isThat ]
      })
    return(DF_object)
  })

setClassUnion("jobn_enrich", c("job_gsea", "job_gsea"))

setMethod("map", signature = c(x = "job_monocle", ref = "jobn_enrich"),
  function(x, ref, pathways,
    data = if (is(ref, "job_gsea")) ref@tables$step1$table_kegg else
      ref@tables$step1$res.kegg[[1]],
    trunc = 20, ...)
  {
    data <- dplyr::filter(data, ID %in% !!pathways)
    refs <- data$geneName_list
    names(refs) <- data$Description
    if (is.numeric(trunc)) {
      names(refs) <- stringr::str_trunc(names(refs), trunc)
    }
    p <- wrap(vis(x, refs, ...))
    p <- .set_lab(p, sig(x), "show pathway genes in pseudotime")
    p
  })

plot_highlight_enrich <- function(table_enrich, highlight, lst_logFC,
  n = 10L, shift = .2, top_by = "p.adjust", sort_by = "GeneRatio", use = "p.adjust")
{
  table_enrich <- dplyr::arrange(table_enrich, !!rlang::sym(top_by))
  table_enrich <- head(table_enrich, n = n)
  reshape <- function(data) {
    split_lapply_rbind(data, ~ ID,
      function(x) {
        symbol <- x$geneName_list[[ 1 ]]
        x <- dplyr::select(x, -dplyr::contains("_list"))
        x <- dplyr::slice(x, rep(1, length(!!symbol)))
        x <- dplyr::mutate(x, Symbol = !!symbol)
        return(x)
      })
  }
  table_enrich <- reshape(table_enrich)
  table_enrich$log2fc <- do.call(dplyr::recode, c(list(table_enrich$Symbol), lst_logFC))
  table_enrich <- dplyr::arrange(table_enrich, dplyr::desc(!!rlang::sym(sort_by)))
  ## data for plot the enrichment score
  data <- dplyr::distinct(table_enrich, Description, !!rlang::sym(use), Count, GeneRatio)
  data <- dplyr::mutate(data, path.p = nrow(data):1)
  ## edges for plot left panel (network of genes with pathway name)
  edges <- dplyr::distinct(table_enrich, Symbol, Description, ID, log2fc)
  edges <- dplyr::filter(edges, ID %in% !!highlight)
  tt <<- edges
  ## nodes for location
  nodes <- tidyr::gather(edges, type, name, -log2fc, -ID)
  nodes <- split(nodes, ~type)
  nodes$Symbol <- dplyr::arrange(nodes$Symbol, log2fc)
  nodes$Symbol <- dplyr::mutate(
    nodes$Symbol, x = -max(data$GeneRatio) * (1 + shift),
    y = seq(1, n, length.out = length(name))
  )
  nodes$Description <- dplyr::distinct(nodes$Description, type, name)
  nodes$Description <- dplyr::mutate(
    nodes$Description, x = 0L - shift * max(data$GeneRatio),
    y = data$path.p[ match(name, data$Description) ]
  )
  nodes <- data.table::rbindlist(nodes, fill = T)
  nodes <- dplyr::relocate(nodes, name)
  ## custom layout
  graph <- fast_layout(edges, dplyr::select(nodes, x, y), nodes = dplyr::select(nodes, -x, -y))
  p <- ggraph(graph) +
    geom_edge_diagonal(aes(x = x, y = y, width = abs(log2fc)),
      color = sample(ggsci::pal_npg()(10), 1), strength = 1, flipped = T, alpha = .1) +
    geom_node_label(
      data = filter(nodes, type == "Symbol"),
      aes(label = name, x = x, y = y, fill = log2fc),
      hjust = 1, size = 2) +
    ## enrichment
    geom_segment(data = data,
      aes(x = 0, xend = GeneRatio, y = path.p, yend = path.p, color = !!rlang::sym(use))) +
    geom_point(data = data,
      aes(x = GeneRatio, y = path.p, color = !!rlang::sym(use), size = Count)) +
    geom_label(data = data,
      aes(x = - shift * max(GeneRatio) / 2, y = path.p, label = stringr::str_wrap(Description, 30)),
      hjust = 1, size = 4) +
    geom_vline(xintercept = 0L, linetype = 4) +
    labs(x = "", y = "", edge_width = "|log2(FC)|", fill = "log2(FC)") +
    rstyle("theme") +
    theme(axis.text.y = element_blank()) +
    scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF") +
    scale_color_gradientn(colours = c(sample(ggsci:::ggsci_db$uchicago$dark, 1), sample(color_set()[1:10], 1))) +
    scale_x_continuous(breaks = round(seq(0, max(data$GeneRatio), length.out = 4), 3),
      limits = c(-max(data$GeneRatio) * (1 + shift) * 1.2, max(data$GeneRatio))) +
    geom_blank()
  p <- wrap(p, 12, 8)
  p
}
