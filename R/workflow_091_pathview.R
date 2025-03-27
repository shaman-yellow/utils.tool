# ==========================================================================
# workflow of pathview
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_pathview <- setClass("job_pathview", 
  contains = c("job"),
  prototype = prototype(
    pg = "pathview",
    info = c("https://github.com/datapplab/pathview"),
    cite = "[@PathviewAnRLuoW2013]",
    method = "",
    tag = "pathview",
    analysis = "Pathview 通路可视化"
    ))

setGeneric("asjob_pathview", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_pathview"))

setMethod("asjob_pathview", signature = c(x = "job_gsea"),
  function(x, gene.level = x@object$symbol, 
    gene.level.name = "symbol", ...)
  {
    methodAdd_onExit("x", "以 R 包 `pathview` {cite_show('PathviewAnRLuoW2013')} 将 `clusterProfiler` 以 GSEA 算法富集的 KEGG 通路可视化。")
    x <- job_pathview(
      data = x@tables$step1$table_kegg,
      annotation = x$annotation, gene.level = gene.level, 
      gene.level.name = gene.level.name,
      species = x$org, ...
    )
    return(x)
  })

setMethod("asjob_pathview", signature = c(x = "job_enrich"),
  function(x, which = 1L, gene.level = NULL, gene.level.name = "hgnc_symbol", ...)
  {
    methodAdd_onExit("x", "以 R 包 `pathview` {cite_show('PathviewAnRLuoW2013')} 将 `clusterProfiler` 富集的 KEGG 通路可视化。")
    x <- job_pathview(
      data = x@tables$step1$res.kegg[[ which ]],
      annotation = x$annotation, gene.level = gene.level, 
      gene.level.name = gene.level.name,
      species = x$organism, ...
    )
    return(x)
  })

setMethod("asjob_pathview", signature = c(x = "job_pathview"),
  function(x, gene.level, gene.level.name = "hgnc_symbol", filter = TRUE, from_which = 1L, ...){
    message("Show pathway again, with extra genes.")
    if (is(gene.level, "job_limma")) {
      snapAdd_onExit("x", "将{snap(feature(gene.level))}映射到 Pathview 通路可视化。")
      gene.level <- dplyr::select(
        gene.level@tables$step2$tops[[ from_which ]],
        !!rlang::sym(gene.level.name), logFC
      )
    }
    data <- x$data
    if (filter) {
      data <- dplyr::filter(data, ID %in% names(x@tables$step1$t.ResultsOfPathviews))
      sets <- unique(unlist(lapply(x@tables$step1$t.ResultsOfPathviews, 
        function(x) {
          x$plot.data.gene$labels
        })))
      if (!length(sets)) {
        stop('!length(sets), no any genes found in t.ResultsOfPathviews.')
      }
      if (is(gene.level, "data.frame")) {
        message("Use first (symbol) and second (logFC) columns of `gene.level`.")
        gene.level <- nl(gname(gene.level[[1]]), gene.level[[2]], FALSE)
      }
      gene.level <- gene.level[ names(gene.level) %in% sets ]
    }
    tmp <- job_enrich(names(gene.level), from = gene.level.name)
    data <- dplyr::mutate(
      data, geneID_list_raw = geneID_list,
      geneID_list = rep(
        list(tmp$annotation$entrezgene_id), nrow(data)
      )
    )
    x <- job_pathview(
      data = data,
      annotation = tmp$annotation, gene.level = gene.level, 
      gene.level.name = gene.level.name,
      species = x$species, ...
    )
    return(x)
  })

job_pathview <- function(data, annotation = NULL,
  gene.level = NULL, gene.level.name = "hgnc_symbol", species = "hsa",
  name = paste0("pathview", gs(Sys.time(), " |:", "_")))
{
  .check_columns(data, c("ID", "geneID_list"), "data")
  if (!is(data[[ "geneID_list" ]], "list")) {
    stop('!is(data[[ "geneID_list" ]], "list"), geneID_list should be a "list" of entrezgene_id.')
  }
  args <- as.list(environment())
  x <- .job_pathview(params = args)
  return(x)
}

setMethod("step0", signature = c(x = "job_pathview"),
  function(x){
    step_message("Prepare your data with function `job_pathview`.")
  })

setMethod("step1", signature = c(x = "job_pathview"),
  function(x, pathways, external = FALSE)
  {
    step_message("Use pathview to visualize reults pathway.")
    if (is.null(x$pathview_dir)) {
      name <- x$pathview_dir <- x$name
    } else {
      name <- x$name <- x$pathview_dir
    }
    data <- x$data
    gene.level <- x$gene.level
    gene.level.name <- x$gene.level.name
    species <- x$species
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
      if (is.null(x$annotation)) {
        stop('is.null(x$annotation), do not known how to match genes.')
      }
      .check_columns(
        x$annotation, c(gene.level.name, "entrezgene_id"), "x$annotation"
      )
      names(gene.level) <- x$annotation$entrezgene_id[ match(names(gene.level),
        x$annotation[[ gene.level.name ]]) ]
      snap <- "通路图中的基因的映射颜色表示基因显著富集，并在数据集中有上调或下调变化趋势。"
    } else {
      snap <- "通路图中的基因的映射颜色表示是否显著富集。"
    }
    dir.create(name, FALSE)
    setwd(name)
    cli::cli_alert_info("pathview::pathview")
    require(pathview)
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
    x <- tablesAdd(x, t.ResultsOfPathviews = res.pathviews)
    p.pathviews <- .pathview_search(
      name, "pathview", x, res.pathviews, snap = snap
    )
    feature(x) <- lapply(p.pathviews, function(x) unique(x$genes))
    x <- plotsAdd(x, p.pathviews)
    snapAdd_onExit("x", "以 `pathview` 探究基因集在通路 {bind(pathways)} 中的上下游关系。")
    .add_internal_job(.job(method = "R package `pathview` used for KEGG pathways visualization", cite = "[@PathviewAnRLuoW2013]"))
    return(x)
  })

setMethod("step2", signature = c(x = "job_pathview"),
  function(x, genes, f.size = 1.5, gp = gpar(col = "steelblue2", lty = "solid", lwd = 1.5))
  {
    genes <- resolve_feature(genes)
    pathways <- x@plots$step1$p.pathviews
    metas <- x@tables$step1$t.ResultsOfPathviews
    x <- snapAdd(x, "从KEGG {bind(names(metas))}中，查找包含{less(genes)}的通路。")
    file_figs <- sapply(names(metas), simplify = FALSE,
      function(name) {
        file <- pathways[[name]]@.Data
        newfile <- file.path(
          dirname(file), 
          paste0("Mark_", tools::file_path_sans_ext(basename(file)), ".pdf")
        )
        data <- metas[[name]]$plot.data.gene
        data <- dplyr::filter(data, labels %in% !!genes)
        if (nrow(data)) {
          img <- magick::image_read(file)
          info <- magick::image_info(img)
          width <- 5
          height <- 5 * info$height / info$width
          pdf(newfile, width, height)
          res <- tryCatch({
            par(mar = rep(0, 4))
            graphics:::plot.raster(grDevices::as.raster(img))
            lapply(seq_len(nrow(data)), 
              function(n) {
                location <- data[n, ]
                w <- location$width / info$width * f.size
                h <- location$height / info$width * f.size
                x <- location$x / info$width
                y <- 1 - (location$y / info$height)
                grid.draw(roundrectGrob(x, y, w, h, gp = gp))    
              })
          }, finally = dev.off())
          file_fig <- pathways[[name]]
          file_fig@.Data <- newfile
          file_fig$lich[["Marked genes"]] <- bind(data$labels)
          return(file_fig)
        } else {
          return(NULL)
        }
      })
    file_figs <- lst_clear0(file_figs)
    if (length(file_figs)) {
      file_figs <- setLegend(file_figs, glue::glue("为KEGG {names(file_figs)}通路图，包含基因上调或下调信息，并标记了{less(genes)}。"))
      x <- plotsAdd(x, p.marked_pathway = file_figs)
    }
    return(x)
  })

setMethod("feature", signature = c(x = "job_pathview"),
  function(x, ref = "all"){
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

.pathview_search <- function(name, search, x, lst, snap = "") {
  figs <- list.files(name, search, full.names = TRUE)
  names <- get_realname(figs)
  p.pathviews <- mapply(figs, names, SIMPLIFY = FALSE,
    FUN = function(x, name) {
      x <- .file_fig(x)
      genes <- dplyr::filter(lst[[ name ]]$plot.data.gene, all.mapped != "")$labels
      attr(x, "genes") <- genes
      attr(x, "lich") <- new_lich(
        list("Interactive figure" = paste0("\\url{https://www.genome.jp/pathway/", name, "}"),
          "Enriched genes" = paste0(unique(genes), collapse = ", ")
        )
      )
      x <- setLegend(x, "KEGG 通路可视化 ({name}) 展示了富集基因在该通路的上下游关系。{snap}")
      return(x)
    })
  names(p.pathviews) <- names
  p.pathviews <- .set_lab(p.pathviews, sig(x), names, "visualization")
  p.pathviews
}
