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
    method = "R package `ClusterProfiler` used for GSEA enrichment",
    cite = "[@ClusterprofilerWuTi2021]",
    tag = "enrich:gsea",
    analysis = "ClusterProfiler GSEA 富集分析"
    ))

setGeneric("asjob_gsea", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_gsea"))

setMethod("asjob_gsea", signature = c(x = "job_limma"),
  function(x, key = 1L, annotation = NULL,
    filter = NULL, from = .guess_symbol(x),
    data = NULL)
  {
    if (is.null(data)) {
      data <- attr(x@tables$step2$tops[[ key ]], "all")
      if (is.null(data)) {
        stop(
          'is.null(attr(x@tables$step2$tops[[ key ]], "all")), can not found all gene set.'
        )
      }
      resolve_feature_snapAdd_onExit("x", as_feature(data[[ from ]], x))
    }
    if (!is.null(filter)) {
      data <- dplyr::filter(data, !!rlang::sym(from) %in% dplyr::all_of(filter))
    }
    data <- dplyr::rename(data, symbol = !!rlang::sym(from))
    x$from <- from
    if (is.null(annotation)) {
      x <- job_gsea(data, use = "symbol")
    } else {
      x <- job_gsea(data, annotation, use = "symbol")
    }
    return(x)
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
  rename <- TRUE
  if (missing(use)) {
    if (any(colnames(topTable) == "symbol")) {
      message("Use column `symbol` as gene input.")
      rename <- FALSE
      use <- "symbol"
    } else if (length(whichs <- grp(colnames(topTable), "_symbol")) >= 1) {
      use <- colnames(topTable)[ whichs[1] ]
      message("Use column `", use, "` as gene input.")
    }
  } else if (use == "symbol") {
    rename <- FALSE
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
  topTable <- dplyr::distinct(topTable, symbol, .keep_all = TRUE)
  used.symbol <- nl(topTable$symbol, topTable$logFC, FALSE)
  if (missing(annotation)) {
    message("Missing `annotation`, use 'hgnc_symbol' to get annotation.")
    if (TRUE) {
      topTable <- dplyr::mutate(topTable, symbol = gname(symbol))
      annotation <- e(
        clusterProfiler::bitr(
          topTable$symbol, fromType = "SYMBOL", toType = "ENTREZID",
          OrgDb = org.Hs.eg.db::org.Hs.eg.db
        )
      )
      annotation <- dplyr::rename(
        annotation, symbol = SYMBOL, entrezgene_id = ENTREZID
      )
    } else {
      mart <- new_biomart()
      message(glue::glue("Use: {bind(head(topTable$symbol))} ..."))
      annotation <- filter_biomart(mart,
        c("hgnc_symbol", "entrezgene_id"), "hgnc_symbol", topTable$symbol)
      annotation <- dplyr::rename(annotation, !!!nl(use, "hgnc_symbol"))
    }
  }
  topTable <- map(topTable, "symbol", annotation, use, "entrezgene_id")
  used.entrezgene_id <- nl(topTable$entrezgene_id, topTable$logFC, FALSE)
  object <- list(symbol = used.symbol, entrezgene_id = used.entrezgene_id)
  keep <- !is.na(names(object$entrezgene_id)) & !duplicated(names(object$entrezgene_id))
  object <- lapply(object, function(x) x[keep])
  x <- .job_gsea(object = object)
  x$annotation <- annotation
  return(x)
}

setMethod("step0", signature = c(x = "job_gsea"),
  function(x){
    step_message("Prepare your data with function `job_gsea`.")
  })

setMethod("step1", signature = c(x = "job_gsea"),
  function(x, OrgDb = org.Hs.eg.db::org.Hs.eg.db, org = "hsa", show = 15,
    order = c("x", "p.adjust", "p.value"), keep_res_go = FALSE)
  {
    step_message("GSEA enrichment.")
    order <- match.arg(order)
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
        GeneRatio = as.double(stringr::str_extract(leading_edge, "[0-9]+")) / 100,
        Count = lengths(geneName_list)
        ), TRUE)
    if (!inherits(table_kegg, "try-error")) {
      table_kegg <- dplyr::as_tibble(table_kegg)
      table_kegg <- .set_lab(table_kegg, sig(x), "GSEA KEGG enrichment data")
      table_kegg <- setLegend(table_kegg, "为 GSEA KEGG 富集分析统计附表。")
      p.kegg <- e(enrichplot::dotplot(
          res.kegg, showCategory = show, orderBy = order,
          decreasing = order == "x"
          ))
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
      table_go <- .set_lab(table_go, sig(x), "GSEA GO enrichment data")
      table_go <- setLegend(table_go, "为 GSEA GO 富集分析统计附表。")
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
    if (keep_res_go) {
      x@params$res.go <- res.go
    }
    x@params$res.kegg <- res.kegg
    x@tables[[ 1 ]] <- namel(table_go, table_kegg)
    p.go <- .set_lab(wrap(p.go), sig(x), "GSEA GO enrichment")
    p.go <- setLegend(p.go, "GSEA GO 富集分析气泡图。")
    p.kegg <- .set_lab(
      wrap(p.kegg, 6.5, 6), sig(x), "GSEA KEGG enrichment"
    )
    p.kegg <- setLegend(p.kegg, "为 GSEA KEGG 富集分析气泡图。")
    x <- methodAdd(x, "以 ClusterProfiler R 包 ({packageVersion('clusterProfiler')}) {cite_show('ClusterprofilerWuTi2021')} 按 GSVA 算法 (`clusterProfiler::gseGO`, `ClusterProfiler::gseKEGG`)，进行 KEGG 和 GO 富集分析 (P-value Cutoff = 0.05) 。")
    x <- snapAdd(x, "以 KEGG、GO 数据集，对基因集富集分析。")
    x@plots[[ 1 ]] <- namel(p.go, p.kegg)
    x$org <- org
    return(x)
  })

setMethod("step2", signature = c(x = "job_gsea"),
  function(x, key = res(x, "id", 1:3),
    highlight = key[1], use = "res.kegg", ...)
  {
    step_message("GSEA visualization for specific pathway")
    if (is.null(key)) {
      return(x)
    }
    obj <- x@params[[ use ]]
    title <- if (length(key) > 1) {
      ""
    } else {
      dplyr::filter(obj@result, ID == key)$Description
    }
    p.code <- wrap(
      e(enrichplot::gseaplot2(obj, key, pvalue_table = FALSE, title = title)),
      7.5, 6
    )
    p.code <- .set_lab(p.code, sig(x), "GSEA plot of the pathways")
    p.code <- setLegend(p.code, "为 GSEA KEGG {title} 富集条码图 (以 `clusterProfiler` 仿 GSEA 软件绘图)。")
    if (missing(key)) {
      data <- head(
        dplyr::arrange(x@tables$step1$table_kegg, p.adjust), n = 3
      )
      sigPaths <- data$Description
      x <- snapAdd(
        x, "按矫正 P 值排序 (BH, p.adjust) ，最显著的三条通路为: {bind(sigPaths)}。"
      )
    } else {
      data <- dplyr::filter(x@tables$step1$table_kegg, ID %in% key)
    }
    x$.feature <- data$geneName_list
    names(x$.feature) <- key
    if (!is.null(highlight)) {
      p.highlight <- plot_highlight_enrich(
        x@tables$step1$table_kegg, highlight, object(x)$symbol, ...
      )
      p.highlight <- .set_lab(p.highlight, sig(x), "KEGG enrichment with enriched genes")
      p.highlight <- setLegend(p.highlight, "为通路富集图，兼部分通路的富集基因表达。")
    } else {
      p.highlight <- NULL
    }
    x@plots[[ 2 ]] <- namel(p.code, p.highlight)
    return(x)
  })

setMethod("step3", signature = c(x = "job_gsea"),
  function(x, db, cutoff = .05, map = NULL, pvalue = FALSE, 
    db_anno = NULL, db_filter = NULL, mode = c(
      "curated gene sets" = "C2",
      "hallmark gene sets" = "H",
      "positional gene sets" = "C1",
      "regulatory target gene sets" = "C3",
      "computational gene sets" = "C4",
      "ontology gene sets" = "C5",
      "oncogenic signature gene sets" = "C6"
      ))
  {
    step_message("Custom database for GSEA enrichment.")
    ## general analysis
    if (missing(db)) {
      mode <- match.arg(mode)
      db_anno <- e(msigdbr::msigdbr(species = "Homo sapiens", category = mode))
      if (!is.null(db_filter)) {
        db_anno <- dplyr::filter(
          db_anno, grpl(gs_description, db_filter, TRUE)
        )
      }
      db <- dplyr::select(db_anno, gs_id, symbol = gene_symbol)
      x <- methodAdd(x, "以 R 包 `msigdbr` ({packageVersion('msigdbr')}) 获取 MSigDB 数据库基因集，用于 clusterProfiler GSEA 富集分析。")
      x <- snapAdd(x, "以 `msigdbr` 获取 {mode} ({names(mode)}) 基因集。")
    }
    if (FALSE) {
      insDb <- lapply(split(db, ~ term),
        function(data) intersect(data$symbol, names(object(x)$symbol)))
      p.pie_insDb <- new_pie(rep(names(insDb), lengths(insDb)))
      table_insDb <- dplyr::filter(db, symbol %in% names(object(x)$symbol))
    }
    ## enrichment
    res.gsea <- e(clusterProfiler::GSEA(object(x)$symbol, TERM2GENE = db,
        pvalueCutoff = cutoff))
    x <- snapAdd(
      x, "使用 {mode} 数据集, 以 `clusterProfiler::GSEA` 对 {less(names(object(x)$symbol))} 富集分析。"
    )
    table_gsea <- dplyr::as_tibble(res.gsea@result)
    table_gsea <- dplyr::mutate(table_gsea,
      geneName_list = strsplit(core_enrichment, "/"),
      Count = lengths(geneName_list),
      GeneRatio = round(as.double(stringr::str_extract(leading_edge, "[0-9]+")) / 100, 2)
    )
    table_gsea <- set_lab_legend(table_gsea, glue::glue("GSEA pathway list of {mode} data"),
        glue::glue("为 GSEA 按 {mode} ({names(mode)}) 数据集富集附表。")
    )
    if (!is.null(db_anno) && all(c("gs_id", "gs_description") %in% colnames(db_anno))) {
      table_gsea <- map(
        table_gsea, "ID", db_anno, "gs_id", "gs_description", col = "Description"
      )
      p.gsea <- plot_kegg(table_gsea)
      p.gsea <- .set_lab(p.gsea, sig(x), glue::glue("GSEA pathway list of {mode}"))
      p.gsea <- setLegend(p.gsea, "为 GSEA 按 {mode} ({names(mode)}) 数据集富集图。")
    } else {
      p.gsea <- NULL
    }
    if (!is.null(map)) {
      p.code <- vis(x, map, res.gsea, table_gsea, pvalue = pvalue)
    } else {
      p.code <- NULL
    }
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x$db_anno <- db_anno
    x@tables[[ 3 ]] <- namel(table_gsea)
    p.code <- .set_lab(p.code, sig(x), "GSEA plot of pathway")
    x@plots[[ 3 ]] <- namel(p.code, p.gsea)
    return(x)
  })

setMethod("vis", signature = c(x = "job_gsea"),
  function(x, map, res.gsea = NULL, table_gsea = x@tables$step3$table_gsea, pvalue = FALSE)
  {
    if (x@step < 3L) {
      stop('x@step < 3L.')
    }
    if (is.null(res.gsea)) {
      res.gsea <- x$res.gsea
    }
    alls <- table_gsea$ID
    if (is.null(alls)) {
      stop('is.null(alls).')
    }
    whichMapped <- which(grepl(map, table_gsea$Description, ignore.case = TRUE))
    map <- alls[ whichMapped ]
    if (!length(map)) {
      message(crayon::red("Not match any pathway, skip plot of 'p.code'."))
      p.code <- NULL
    } else {
      p.code <- wrap(e(
          enrichplot::gseaplot2(
            x$res.gsea, map, pvalue_table = pvalue,
            title = bind(
              stringr::str_wrap(table_gsea$Description[whichMapped], 80),
              co = "\n"
              ))), 7.5, 6)
    }
    if (length(map) > 1) {
      title <- ""
    } else {
      title <- table_gsea$Description[ whichMapped ]
    }
    ids <- table_gsea$ID[whichMapped]
    p.code <- set_lab_legend(
      p.code, glue::glue("{sig(x)} GSEA plot {bind(ids, co = '_')}"),
      glue::glue("为 GSEA plot {bind(ids)} {title}。")
    )
    p.code
  })

# setMethod("filter", signature = c(x = "job_gsea"),
  # function(x, ref, use = c("entrezgene_id", "symbol")){
  #   use <- match.arg(use)
  #   isThat <- names(object(x)[[ use ]]) %in% ref
  #   object(x) <- lapply(object(x),
  #     function(x) {
  #       x[ isThat ]
  #     })
  #   return(x)
  # })

setClassUnion("jobn_enrich", c("job_gsea", "job_enrich"))

setMethod("filter", signature = c(x = "jobn_enrich"),
  function(x, pattern, ..., use = c("kegg", "go", "gsea"), 
    which = 1, genes = NULL, return_type = c("job", "data.frame"), step = x@step)
  {
    message("Search genes in enriched pathways.")
    return_type <- match.arg(return_type)
    if (return_type == "job") {
      if (is(genes, "feature")) {
        genes <- resolve_feature(genes)
      }
    }
    if ((missing(pattern) || is.null(pattern)) && !is.null(genes)) {
      pattern <- paste0(paste0("^", genes, "$"), collapse = "|")
    }
    if (is(x, "job_gsea")) {
      prefix <- "table_"
    } else {
      prefix <- "res."
    }
    if (is.character(use)) {
      use <- match.arg(use)
      if (use == "kegg") {
        data <- x@tables$step1[[ paste0(prefix, "kegg") ]]
      } else if (use == "go") {
        data <- x@tables$step1[[ paste0(prefix, "go") ]]
      } else if (use == "gsea") {
        data <- x@tables$step3[[ paste0(prefix, "gsea") ]]
      }
      if (is(x, "job_enrich")) {
        data <- data[[ which ]]
      }
    } else if (is(use, "data.frame")) {
      message("Custom passed `data` for searching.")
      data <- use
    }
    isThat <- vapply(data$geneName_list, FUN.VALUE = logical(1),
      function(x) {
        any(grpl(x, pattern, ...))
      })
    data <- dplyr::filter(data, !!isThat)
    if (!is.null(genes)) {
      data[[ "match_genes" ]] <- lapply(data[[ "geneName_list" ]],
        function(x) {
          x[x %in% genes]
        })
    }
    data <- .set_lab(data, sig(x), "filter by match genes ", step)
    data <- setLegend(data, "以{less(genes)}筛选到的通路。")
    if (return_type == "data.frame") {
      return(data)
    }
    if (return_type == "job") {
      if (nrow(data)) {
        init(snap(x)) <- TRUE
        x <- snapAdd(
          x, "从富集结果中筛选包含{less(genes)}的通路。", step = step
        )
        x <- snapAdd(
          x, "筛选到{less(data$Description)}。", step = step
        )
      }
      x[[ paste0("filtered_pathways_", use, "_", step) ]] <- data
      return(x)
    }
  })

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

setMethod("res", signature = c(x = "job_gsea", ref = "character"),
  function(x, ref = c("id", "des", "p", "adj"),
    which = 1, from = c("kegg", "go"))
  {
    type <- match.arg(ref)
    type <- switch(
      type, id = "ID", des = "Description", p = "pvalue", adj = "p.adjust"
    )
    from <- match.arg(from)
    data <- x@tables$step1[[ paste0("table_", from) ]]
    data[[ type ]][ which ]
  })


plot_highlight_enrich <- function(table_enrich, highlight, lst_logFC,
  n = 10L, shift = .2, top_by = "p.adjust", sort_by = "GeneRatio", use = "p.adjust")
{
  if (!all(highlight %in% head(table_enrich$ID, n = 10))) {
    message(glue::glue("All not in top {n}, so just keep the `highlight`."))
    table_enrich <- dplyr::filter(table_enrich, ID %in% !!highlight)
  }
  table_enrich <- dplyr::arrange(table_enrich, !!rlang::sym(top_by))
  table_enrich <- head(table_enrich, n = n)
  n <- nrow(table_enrich)
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
  table_enrich <- dplyr::arrange(table_enrich,
    if (sort_by == "GeneRatio") dplyr::desc(!!rlang::sym(sort_by))
    else !!rlang::sym(sort_by)
  )
  ## data for plot the enrichment score
  data <- dplyr::distinct(table_enrich, Description, !!rlang::sym(use), Count, GeneRatio)
  data <- dplyr::mutate(data, path.p = nrow(data):1)
  ## edges for plot left panel (network of genes with pathway name)
  edges <- dplyr::distinct(table_enrich, Symbol, Description, ID, log2fc)
  edges <- dplyr::filter(edges, ID %in% !!highlight)
  ## nodes for location
  nodes <- tidyr::gather(edges, type, name, -log2fc, -ID)
  nodes <- split(nodes, ~type)
  nodes$Symbol <- dplyr::distinct(nodes$Symbol, name, log2fc, type)
  nodes$Symbol <- dplyr::arrange(nodes$Symbol, log2fc)
  nodes$Symbol <- dplyr::mutate(
    nodes$Symbol, x = -max(data$GeneRatio) * (1 + shift),
    y = seq(1, n, length.out = length(name)),
    hjust = 1
  )
  if (nrow(nodes$Symbol) > 50) {
    nodes$Symbol <- dplyr::mutate(nodes$Symbol, hjust = seq_len(nrow(nodes$Symbol)) %% 2)
  }
  nodes$Description <- dplyr::distinct(nodes$Description, type, name)
  nodes$Description <- dplyr::mutate(
    nodes$Description, x = 0L - shift * max(data$GeneRatio),
    y = data$path.p[ match(name, data$Description) ]
  )
  nodes <- data.table::rbindlist(nodes, fill = TRUE)
  nodes <- dplyr::relocate(nodes, name)
  ## custom layout
  graph <- fast_layout(edges, dplyr::select(nodes, x, y), nodes = dplyr::select(nodes, -x, -y))
  p <- ggraph(graph) +
    geom_edge_diagonal(aes(x = x, y = y, width = abs(log2fc), edge_color = node2.name),
      strength = 1, flipped = TRUE, alpha = .25) +
    geom_node_label(
      data = filter(nodes, type == "Symbol"),
      aes(label = name, x = x, y = y, fill = log2fc, hjust = hjust),
      size = 2) +
    ## enrichment
    geom_segment(data = data,
      aes(x = 0, xend = GeneRatio, y = path.p, yend = path.p, color = !!rlang::sym(use))) +
    geom_point(data = data,
      aes(x = GeneRatio, y = path.p, color = !!rlang::sym(use), size = Count)) +
    geom_label(data = data,
      aes(x = - shift * max(GeneRatio) / 2, y = path.p, label = stringr::str_wrap(Description, 50)),
      hjust = 1, size = 4) +
    geom_vline(xintercept = 0L, linetype = 4) +
    labs(x = "GeneRatio", y = "", edge_width = "|log2(FC)|",
      fill = "log2(FC)", edge_color = "Highlight Pathways") +
    rstyle("theme") +
    theme(axis.text.y = element_blank()) +
    scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF") +
    scale_edge_color_manual(values = color_set()) +
      scale_color_gradientn(colours = color_set2()) +
      scale_x_continuous(breaks = round(seq(0, max(data$GeneRatio), length.out = 4), 3),
        limits = c(-max(data$GeneRatio) * (1 + shift) * 1.2, max(data$GeneRatio))) +
      geom_blank()
    p <- wrap(p, 12, 8)
    p
}

color_set2 <- function() {
  unname(c(sample(ggsci:::ggsci_db$uchicago$dark, 1), sample(color_set()[1:5], 1)))
}

# setMethod("step3", signature = c(x = "job_gsea"),
#   function(x, pathways, species = x$org,
#     name = paste0("pathview", gs(Sys.time(), " |:", "_")),
#     search = "pathview")
#   {
#     step_message("Use pathview to visualize reults pathway.")
#     require(pathview)
#     data <- x@tables$step1$table_kegg
#     if (is.null(x$pathview_dir)) {
#       x$pathview_dir <- name
#     } else {
#       name <- x$pathview_dir
#     }
#     dir.create(name, FALSE)
#     setwd(name)
#     cli::cli_alert_info("pathview::pathview")
#     tryCatch({
#       res.pathviews <- sapply(pathways, simplify = FALSE,
#         function(pathway) {
#           data <- dplyr::filter(data, ID == !!pathway)
#           pathway <- gs(data$ID, "^[a-zA-Z]*", "")
#           genes <- as.character(unlist(data$geneID_list))
#           genes.all <- x@object$entrezgene_id
#           genes <- genes.all[ match(genes, names(genes.all)) ]
#           res.pathview <- try(
#             pathview::pathview(gene.data = genes,
#               pathway.id = pathway, species = species,
#               keys.align = "y", kegg.native = TRUE, same.layer = FALSE,
#               key.pos = "topright", na.col = "grey90")
#           )
#           if (inherits(res.pathview, "try-error")) {
#             try(dev.off(), silent = TRUE)
#           }
#           return(res.pathview)
#         })
#     }, finally = {setwd("../")})
#     x@tables[[ 3 ]] <- namel(res.pathviews)
#     p.pathviews <- .pathview_search(name, search, x, res.pathviews)
#     x@plots[[ 3 ]] <- namel(p.pathviews)
#     .add_internal_job(.job(method = "R package `pathview` used for KEGG pathways visualization", cite = "[@PathviewAnRLuoW2013]"))
#     return(x)
#   })
#
