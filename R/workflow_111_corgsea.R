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
    x <- job_corgsea(data, ref, method = method)
    return(x)
  })

setMethod("asjob_corgsea", signature = c(x = "job_limma"),
  function(x, ref, method = "spearman"){
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    data <- x$normed_data$E
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    x <- job_corgsea(data, ref, method = method)
  })

job_corgsea <- function(data, ref, method = "spearman") {
  if (!is(ref, "feature")) {
    stop('!is(ref, "feature").')
  }
  which.ref <- rownames(data) %in% ref
  data.ref <- t(data[which.ref, ])
  data.others <- t(data[!which.ref, ])
  cors <- e(stats::cor(data.ref, data.others, method = method))
  cors <- apply(cors, 1, sort, decreasing = TRUE, simplify = FALSE)
  x <- .job_corgsea(object = cors)
  x <- methodAdd(
    x, "以相关系数构建全基因排序列表，在此基础上实施基因集富集分析（GSEA），以识别与诊断基因表达模式协同变化的功能通路。该策略能够将有限的诊断基因扩展至其相关的基因网络层面，从而揭示其潜在的生物学过程及分子机制，提高结果的生物学解释性与稳健性。"
  )
  x <- snapAdd(x, "计算{snap(ref)}与其他基因的 Spearman 相关性系数，并以该系数为排序依据对全基因进行从大到小的排序。")
}

setMethod("step0", signature = c(x = "job_corgsea"),
  function(x){
    step_message("Prepare your data with function `job_corgsea`.")
  })

setMethod("step1", signature = c(x = "job_corgsea"),
  function(x, db, cutoff = .05, pattern = NULL, pvalue = FALSE, 
    cutoff.nes = 1, db_anno = NULL, rerun = FALSE, mode = c(
      "curated gene sets" = "C2",
      "hallmark gene sets" = "H",
      "positional gene sets" = "C1",
      "regulatory target gene sets" = "C3",
      "computational gene sets" = "C4",
      "ontology gene sets" = "C5",
      "oncogenic signature gene sets" = "C6"
      ), mode_sub = "CP")
  {
    step_message("Custom database for GSEA enrichment.")
    ## general analysis
    if (missing(db)) {
      mode <- match.arg(mode)
      x <- .set_msig_db(x, mode, mode_sub)
      db <- x$msig_db
    }
    if (is.null(db_anno)) {
      db_anno <- x$db_anno
    }
    cli::cli_h1("clusterProfiler::GSEA")
    dir.create("tmp", FALSE)
    all.gsea <- pbapply::pbsapply(names(object(x)), simplify = FALSE,
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
        table_gsea <- set_lab_legend(table_gsea, glue::glue("GSEA pathway list of {name} data"),
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
          p.gsea <- setLegend(p.gsea, "基因 {name} 的 GSEA 按 {mode} ({names(mode)}) 数据集富集图。")
        } else {
          p.gsea <- NULL
        }
        return(namel(table_gsea, p.gsea, res.gsea))
      }
    )
    p.gsea <- lapply(all.gsea, function(x) x$p.gsea)
    table_gsea <- lapply(all.gsea, function(x) x$table_gsea)
    if (length(all.gsea) < 6) {
      snaps <- vapply(table_gsea, FUN.VALUE = character(1),
        function(data) {
          .stat_table_by_pvalue(data, n = 10, use.p = "p.adjust")
        })
      snaps <- glue::glue("基因 {names(table_gsea)} 相关富集一共富集到 {snaps}")
      x <- snapAdd(x, "{bind(snaps, co = '\n\n')}")
    }
    res.gsea <- lapply(all.gsea, function(x) x$res.gsea)
    x <- methodAdd(x, "使用 {mode} 数据集, 以 R 包 `clusterProfiler` ⟦pkgInfo('clusterProfiler')⟧ 对基因列表富集分析。富集设定阈值 ⟦mark$blue('adjust P value (FDR) &lt; {cutoff}，|NES| &gt; {cutoff.nes}')⟧。")
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x$db_anno <- db_anno
    x <- tablesAdd(x, table_gsea)
    x <- plotsAdd(x, p.gsea)
    return(x)
  })

setMethod("step2", signature = c(x = "job_corgsea"),
  function(x, top = 10, intersect = TRUE){
    step_message("Select and Visualization")
    ins <- lapply(x@tables$step1$table_gsea,
      function(data) {
        head(data$ID, n = top)
      })
    ins <- ins(lst = ins)
    if (!length(ins)) {
      stop('!length(ins), no intersect found.')
    }
    insName <- dplyr::filter(x@tables$step1$table_gsea[[1]], ID %in% !!ins)$Description
    x$.feature <- list()
    p.codes <- sapply(names(x$res.gsea), simplify = FALSE, 
      function(name) {
        data <- x@tables$step1$table_gsea[[name]]
        dataTop <- head(data, n = top)
        idTop <- dataTop$ID
        x$.feature[[name]] <<- dataTop$Description
        p.code <- vis(
          x, map = idTop, res.gsea = x$res.gsea[[name]],
          table_gsea = data, .name = name
        )
        p.code
      })
    alls <- names(x$res.gsea)
    x <- snapAdd(x, "选取 {bind(alls)} 的 Top {top} 富集通路{aref(p.codes)}。")
    x <- snapAdd(x, "对 {bind(alls)} 的 Top {top} 通路取交集，得到 {length(ins)} 个通路：{bind(insName)}。")
    x <- plotsAdd(x, p.codes)
    return(x)
  })

setClassUnion("job_gseaSet", c("job_corgsea", "job_gsea"))

setMethod("vis", signature = c(x = "job_gseaSet"),
  function(x, pattern, map = NULL, res.gsea = NULL, table_gsea = NULL,
    mode = c("kegg", "gsea"), pvalue = FALSE, .name = "", merge = TRUE)
  {
    mode <- match.arg(mode)
    if (is.null(res.gsea)) {
      res.gsea <- x[[ glue::glue("res.{mode}") ]]
    }
    if (is.null(table_gsea)) {
      if (mode == "kegg") {
        table_gsea <- x@tables$step1$table_kegg
      } else if (mode == "gsea") {
        table_gsea <- x@tables$step3$table_gsea
      }
    }
    if (is.null(res.gsea) || is.null(table_gsea)) {
      stop('is.null(res.gsea) || is.null(table_gsea).')
    }
    alls <- table_gsea$ID
    if (is.null(alls)) {
      stop('is.null(alls).')
    }
    if (is.null(map)) {
      whichMapped <- which(grepl(pattern, table_gsea$Description, ignore.case = TRUE))
      map <- alls[ whichMapped ]
    } else {
      whichMapped <- which(table_gsea$ID %in% map)
    }
    if (!length(map)) {
      message(crayon::red("Not match any pathway, skip plot of 'p.code'."))
      p.code <- NULL
    } else {
      if (merge) {
        fun_plot <- function() {
          res.gsea@result <- map(
            res.gsea@result, "ID", table_gsea, "ID", "Description", col = "Description"
          )
          plst <- enrichplot::gseaplot2(res.gsea, map, pvalue_table = pvalue)
          plst[[1]] <- plst[[1]] +
            guides(color = guide_legend(ncol = 2)) +
            theme(legend.position = "top")
          print(plst)
        }
        p.code <- grid.grabExpr(fun_plot())
      } else {
        p.code <- sapply(map, simplify = FALSE,
          function(key) {
            title <- dplyr::filter(table_gsea, ID == key)$Description
            grob <- grid.grabExpr(
              print(enrichplot::gseaplot2(res.gsea, key, pvalue_table = pvalue, title = title))
            )
            wrap(grob, 5, 4)
          })
        if (length(map) > 1) {
          layout <- calculate_layout(length(map))
          p.code <- patchwork::wrap_plots(
            lapply(p.code, function(x) x@data), ncol = layout[["cols"]]
          )
          p.code <- wrap_layout(p.code, layout, 3)
        } else {
          p.code <- p.code[[1]]
        }
      }
    }
    ids <- table_gsea$ID[whichMapped]
    des <- table_gsea$Description[whichMapped]
    p.code <- set_lab_legend(
      p.code,
      glue::glue("{sig(x)} GSEA plot {.name}"),
      glue::glue("GSEA 富集条码图 ({.name})|||第一部分是 ES 折线图，离垂直距离 x = 0 轴最远的峰值便是基因集的 ES 值，正值表示基因集在列表的顶部富集，负值表示基因集在列表的底部富集。第二部分为基因集成员位置图，用竖线标记了基因集中各成员出现在基因排序列表中的位置。第三部分是排序后所有基因 rank 值的分布，以灰色面积图显展示，左侧灰色 rank 值为正即与该关键基因呈正相关，右侧 rank 值为负是负相关。")
    )
    p.code
  })


.set_msig_db <- function(x, mode, sub) {
  if (packageVersion("msigdbr") < "10.0.0") {
    db_anno <- e(msigdbr::msigdbr(species = "Homo sapiens", category = mode))
  } else {
    db_anno <- e(msigdbr::msigdbr(species = "Homo sapiens", collection = mode))
  }
  x <- methodAdd(
    x, "以 R 包 `msigdbr` ⟦pkgInfo('msigdbr')⟧ 获取 MSigDB 数据库{mode}基因集。该基因集包含多个子集：{try_snap(db_anno, 'gs_subcat', 'gs_name')}。"
  )
  if (!is.null(sub)) {
    select <- c("CP:REACTOME", "CP:KEGG", "CP:WIKIPATHWAYS")
    db_anno <- dplyr::filter(db_anno, gs_subcat %in% !!select)
    x <- methodAdd(x, "选取 {bind(select)} 子集用于后续分析。")
  }
  x$db_anno <- db_anno
  x$msig_db <- dplyr::select(db_anno, gs_id, symbol = gene_symbol)
  return(x)
}

