# ==========================================================================
# workflow of corgsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_corgsea <- setClass("job_corgsea", 
  contains = c("job"),
  prototype = prototype(
    pg = "corgsea",
    info = c("http://www.gsea-msigdb.org/gsea/downloads.jsp"),
    method = "R package `ClusterProfiler` used for GSEA enrichment",
    cite = "[@ClusterprofilerWuTi2021]",
    tag = "enrich:gsea",
    analysis = "GSEA分析"
    ))

setGeneric("asjob_corgsea",
  function(x, ...) standardGeneric("asjob_corgsea"))

setMethod("asjob_corgsea", signature = c(x = "job_deseq2"),
  function(x, ref, method = "spearman")
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    if (any(!ref %in% rownames(object(x)))) {
      stop('any(!ref %in% rownames(object(x))).')
    }
    data <- SummarizedExperiment::assay(x$vst)
    which.ref <- rownames(data) %in% ref
    data.ref <- t(data[which.ref, ])
    data.others <- t(data[!which.ref, ])
    cors <- e(stats::cor(data.ref, data.others, method = "spearman"))
    cors <- apply(cors, 1, sort, decreasing = TRUE, simplify = FALSE)
    x <- .job_corgsea(object = cors)
    return(x)
  })

setMethod("step0", signature = c(x = "job_corgsea"),
  function(x){
    step_message("Prepare your data with function `job_corgsea`.")
  })

setMethod("step1", signature = c(x = "job_corgsea"),
  function(x, db, cutoff = .05, pattern = NULL, pvalue = FALSE, 
    cutoff.nes = 1, db_anno = NULL, db_filter = NULL, rerun = FALSE, mode = c(
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
      if (packageVersion("msigdbr") < "10.0.0") {
        db_anno <- e(msigdbr::msigdbr(species = "Homo sapiens", category = mode))
      } else {
        db_anno <- e(msigdbr::msigdbr(species = "Homo sapiens", collection = mode))
      }
      if (!is.null(db_filter)) {
        db_anno <- dplyr::filter(
          db_anno, grpl(gs_description, db_filter, TRUE)
        )
      }
      db <- dplyr::select(db_anno, gs_id, symbol = gene_symbol)
      x <- methodAdd(x, "以 R 包 `msigdbr` ({packageVersion('msigdbr')}) 获取 MSigDB 数据库基因集，用于 clusterProfiler GSEA 富集分析。")
      x <- snapAdd(x, "以 `msigdbr` 获取 {mode} ({names(mode)}) 基因集。")
    }
    cli::cli_h1("clusterProfiler::GSEA")
    dir.create("tmp", FALSE)
    res.gsea <- pbapply::pbsapply(names(object(x)), simplify = FALSE,
      function(name) {
        glist <- object(x)[[name]]
        args <- list(geneList = glist, TERM2GENE = db, pvalueCutoff = cutoff)
        res.gsea <- expect_local_data(
          "tmp", "gsea", clusterProfiler::GSEA, args, rerun = rerun
        )
        table_gsea <- dplyr::as_tibble(res.gsea@result)
        table_gsea <- dplyr::filter(table_gsea, abs(NES) > cutoff.nes)
        table_gsea <- dplyr::mutate(table_gsea,
          geneName_list = strsplit(core_enrichment, "/"),
          Count = lengths(geneName_list),
          GeneRatio = round(as.double(stringr::str_extract(leading_edge, "[0-9]+")) / 100, 2)
        )
        table_gsea <- set_lab_legend(table_gsea, glue::glue("GSEA pathway list of {mode} data"),
          glue::glue("为基因 {name} 的 GSEA 按 {mode} ({names(mode)}) 数据集富集附表。")
        )
        if (!is.null(db_anno) && all(c("gs_id", "gs_description") %in% colnames(db_anno))) {
          table_gsea <- map(
            table_gsea, "ID", db_anno, "gs_id", "gs_description", col = "Description"
          )
          table_gsea <- dplyr::mutate(
            table_gsea, Description = stringr::str_wrap(Description, 80)
          )
          p.gsea <- plot_kegg(table_gsea)
          p.gsea <- .set_lab(p.gsea, sig(x), glue::glue("Gene {name} GSEA pathway list of {mode}"))
          p.gsea <- setLegend(p.gsea, "为基因 {name} 的 GSEA 按 {mode} ({names(mode)}) 数据集富集图。")
        } else {
          p.gsea <- NULL
        }
        return(namel(table_gsea, p.gsea))
      }
    )
    p.gsea <- lapply(res.gsea, function(x) x$p.gsea)
    table_gsea <- lapply(res.gsea, function(x) x$table_gsea)
    x <- methodAdd(x, "使用 {mode} 数据集, 以 `clusterProfiler::GSEA` 对基因列表富集分析。富集设定阈值 adjust P value (FDR) &lt; 0.05，|NES| &gt; 1。")
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x$db_anno <- db_anno
    x <- tablesAdd(x, table_gsea)
    x <- plotsAdd(x, p.gsea)
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_corgsea"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
