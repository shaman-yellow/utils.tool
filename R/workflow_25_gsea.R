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
    info = c("http://www.gsea-msigdb.org/gsea/downloads.jsp")
    ))

setGeneric("asjob_gsea", 
  function(x, ...) standardGeneric("asjob_gsea"))

setMethod("asjob_gsea", signature = c(x = "job_limma"),
  function(x, key = 1L, annotation = x@params$normed_data$genes,
    filter = NULL)
  {
    if (is.null(filter)) {
      data <- x@tables$step2$tops[[ key ]]
    } else {
      data <- dplyr::filter(data, hgnc_symbol %in% dplyr::all_of(filter))
    }
    job_gsea(data, annotation)
  })

setMethod("asjob_gsea", signature = c(x = "job_seurat"),
  function(x, contrast.pattern = NULL, marker.list = x@params$contrasts)
  {
    if (!is.null(contrast.pattern)) {
      topTable <- dplyr::filter(marker.list, grpl(contrast, !!contrast.pattern))
    } else {
      topTable <- marker.list
    }
    job_gsea(dplyr::relocate(topTable, hgnc_symbol = gene, logFC = avg_log2FC))
  })

job_gsea <- function(topTable, annotation)
{
  .check_columns(topTable, c("hgnc_symbol", "logFC"))
  topTable <- dplyr::select(topTable, hgnc_symbol, logFC)
  topTable <- dplyr::filter(topTable, !is.na(hgnc_symbol) & hgnc_symbol != "")
  topTable <- dplyr::arrange(topTable, dplyr::desc(logFC))
  topTable <- dplyr::distinct(topTable, hgnc_symbol, .keep_all = T)
  used.hgnc_symbol <- nl(topTable$hgnc_symbol, topTable$logFC, F)
  if (missing(annotation)) {
    mart <- new_biomart()
    annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol", topTable$hgnc_symbol)
  }
  topTable <- map(topTable, "hgnc_symbol", annotation, "hgnc_symbol", "entrezgene_id")
  used.entrezgene_id <- nl(topTable$entrezgene_id, topTable$logFC, F)
  object <- list(hgnc_symbol = used.hgnc_symbol, entrezgene_id = used.entrezgene_id)
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
    ## kegg
    res.kegg <- e(clusterProfiler::gseKEGG(object(x)$entrezgene_id, organism = org))
    fun <- function(sets) {
      lapply(sets,
        function(set) {
          from_ids <- x$annotation$entrezgene_id
          to_names <- x$annotation$hgnc_symbol
          to_names[ match(set, from_ids) ]
        })
    }
    table_kegg <- dplyr::mutate(res.kegg@result,
      geneID_list = strsplit(core_enrichment, "/"),
      geneName_list = fun(geneID_list),
      GeneRatio = stringr::str_extract(leading_edge, "[0-9]+"),
      Count = lengths(geneName_list)
    )
    table_kegg <- dplyr::as_tibble(table_kegg)
    p.kegg <- e(enrichplot::dotplot(res.kegg))
    ## go
    if (is.null(x$res.go)) {
      res.go <- e(clusterProfiler::gseGO(object(x)$hgnc_symbol, ont = "ALL", OrgDb = OrgDb,
          keyType = "SYMBOL"))
    } else {
      res.go <- x$res.go
    }
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
    x@params$res.go <- res.go
    x@params$res.kegg <- res.kegg
    x@tables[[ 1 ]] <- namel(table_go, table_kegg)
    x@plots[[ 1 ]] <- namel(p.go, p.kegg)
    return(x)
  })

setMethod("step2", signature = c(x = "job_gsea"),
  function(x, key, use = "res.kegg"){
    step_message("GSEA visualization for specific pathway")
    obj <- x@params[[ use ]]
    p.code <- wrap(e(enrichplot::gseaplot2(obj, key, pvalue_table = F)), 7.5, 6)
    x@plots[[ 2 ]] <- namel(p.code)
    return(x)
  })

setMethod("step3", signature = c(x = "job_gsea"),
  function(x, pathways, species = "hsa",
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
  function(x, db, cutoff = .05, map = NULL){
    step_message("Custom database for GSEA enrichment.")
    ## general analysis
    insDb <- lapply(split(db, ~ term),
      function(data) intersect(data$symbol, names(object(x)$hgnc_symbol)))
    p.pie_insDb <- new_pie(rep(names(insDb), lengths(insDb)))
    table_insDb <- dplyr::filter(db, symbol %in% names(object(x)$hgnc_symbol))
    ## enrichment
    res.gsea <- e(clusterProfiler::GSEA(object(x)$hgnc_symbol, TERM2GENE = db,
        pvalueCutoff = cutoff))
    table_gsea <- dplyr::as_tibble(res.gsea@result)
    if (!is.null(map)) {
      alls <- table_gsea$ID
      map <- alls[ grepl(map, alls) ]
      p.code <- wrap(e(enrichplot::gseaplot2(res.gsea, map)), 7.5, 6)
    } else {
      p.code <- NULL
    }
    table_gsea <- dplyr::mutate(table_gsea,
      geneName_list = strsplit(core_enrichment, "/"),
      GeneRatio = stringr::str_extract(leading_edge, "[0-9]+"),
      Count = lengths(geneName_list)
    )
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x@tables[[ 3 ]] <- namel(table_gsea, table_insDb)
    x@plots[[ 3 ]] <- namel(p.code, p.pie_insDb)
    return(x)
  })
