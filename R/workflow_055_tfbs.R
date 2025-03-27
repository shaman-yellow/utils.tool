# ==========================================================================
# workflow of tfbs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tfbs <- setClass("job_tfbs", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://tfbsdb.systemsbiology.net/"),
    cite = "[@CausalMechanisPlaisi2016]",
    method = "The `Transcription Factor Target Gene Database` (<https://tfbsdb.systemsbiology.net/>) was used for discovering relationship between transcription factors and genes. ",
    tag = "tf",
    analysis = "TFBS 转录因子结合位点数据"
    ))

job_tfbs <- function(genes)
{
  .job_tfbs(object = rm.no(genes))
}

setMethod("step0", signature = c(x = "job_tfbs"),
  function(x){
    step_message("Prepare your data with function `job_tfbs`.")
  })

setMethod("step1", signature = c(x = "job_tfbs"),
  function(x, db = .prefix("tfbs/tfbs.rds", "db"), db_anno = .prefix("tfbs/anno.rds", "db"),
    file_anno = .prefix("tfbs/humanTFMotifEntrezMappings.xlsx", "db"), cl = NULL)
  {
    step_message("Obtain Transcription Factor.")
    db <- new_db(db, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      res <- pbapply::pbsapply(db@query, simplify = FALSE,
        function(gene) {
          url <- paste0("https://tfbsdb.systemsbiology.net/searchgene?searchterm=", gene)
          res <- RCurl::getURL(url)
          res <- try(get_table.html(res), TRUE)
          if (inherits(res, "try-error")) {
            message("\nCan not get results of ", gene, ", continue to next.")
            return()
          }
          res <- lapply(res,
            function(x) {
              if (is(x, "df"))
                as_tibble(x)
            })
          res$tfbs
        })
      res <- frbind(res, idcol = TRUE, fill = TRUE)
      db <- upd(db, res)
    }
    res <- dplyr::filter(db@db, .id %in% object(x))
    if (TRUE) {
      if (!file.exists(file_anno)) {
        cont <- e(RCurl::getURLContent(
            "https://tfbsdb.systemsbiology.net/static/downloads/humanTFMotifEntrezMappings.xlsx"
            ))
        writeBin(as.vector(cont), file_anno)
      }
      anno <- fxlsx(file_anno)
      res <- map(res, "Motif", anno, "Motif.Name", "Gene", col = "Symbol")
      if (any(is.na(res$Symbol))) {
        ids <- unique(res$Motif[ is.na(res$Symbol) ])
        db <- new_db(db_anno, ".id")
        db <- not(db, ids)
        if (length(db@query)) {
          message("Need query: ", length(db@query), " (total: ", length(unique(res$Motif)), ")")
          anno <- pbapply::pbsapply(db@query, simplify = FALSE, cl = cl,
            function(id) {
              url <- paste0("https://tfbsdb.systemsbiology.net/searchtf?searchterm=", id)
              res <- try_get_url(url)
              lst <- suppressMessages(try(get_table.html(res, which = 1:2), TRUE))
              if (inherits(lst, "try-error")) {
                message("Skip: ", id, " as error occurred.")
                return(data.frame(Symbol = NA_character_))
              }
              res <- lapply(lst,
                function(x) {
                  if (is(x, "df"))
                    as_tibble(x)
                })
              des <- strx(res[[2]][1,1]$V1, "^[^\n]+")
              Sys.sleep(1)
              data.frame(Symbol = strx(des, "^[^ ]+"), des = des)
            })
          anno <- frbind(anno, idcol = TRUE, fill = TRUE)
          db <- upd(db, anno)
        }
        anno <- dplyr::filter(db@db, .id %in% ids)
        res$Symbol[ is.na(res$Symbol) ] <- anno$Symbol[match(res$Motif[ is.na(res$Symbol) ], anno$.id)]
      }
      if (TRUE) {
        res$Symbol[ is.na(res$Symbol) ] <- toupper(strx(res$Motif[is.na(res$Symbol)], "[^_]{2,}"))
      }
    }
    res <- dplyr::relocate(res, target = .id, TF_symbol = Symbol)
    res <- .set_lab(res, sig(x), "Transcription Factor binding sites")
    res <- setLegend(res, "为转录因子与基因启动子结合位点数据。")
    x@tables[[ 1 ]] <- namel(res)
    x <- methodAdd(x, "以 `RSelenium` 抓取 `Transcription Factor Target Gene Database` 数据库  (<https://tfbsdb.systemsbiology.net/>) {cite_show('CausalMechanisPlaisi2016')} 基因与转录因子结合数据。")
    x <- snapAdd(x, "从 `Transcription Factor Target Gene Database` 数据库获取与 {less(object(x))} 结合的转录因子数据。")
    res <- dplyr::distinct(res, target, TF_symbol)
    feature(x) <- split(res$TF_symbol, res$target)
    return(x)
  })

setMethod("map", signature = c(x = "job_tfbs", ref = "feature"),
  function(x, ref, name = "TF", region = FALSE){
    data <- x@tables$step1$res
    init(meth(x)) <- init(snap(x)) <- FALSE
    if (is(ref, "feature_list")) {
      data_ref <- as_df.lst(ref, "target", "TF")
      data <- tbmerge(
        data, data_ref, by.x = c("target", "TF_symbol"), 
        by.y = c("target", "TF")
      )
    } else {
      data <- dplyr::filter(data, TF_symbol %in% ref)
    }
    data <- .set_lab(data, sig(x), glue::glue("filtered {name} binding site"))
    data <- map_gene(data, "target", "SYMBOL", "ENTREZID", get = "target_entrez")
    ## hg19 chrom
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) {
      e(BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene"))
    }
    genes_txdb <- suppressMessages(GenomicFeatures::genes(
      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    ))
    genes_txdb <- data.frame(genes_txdb)
    data <- map(
      data, "target_entrez", genes_txdb, "gene_id", "seqnames", col = "chr"
    )
    if (!requireNamespace("Gviz", quietly = TRUE)) {
      e(BiocManager::install("Gviz"))
    }
    data <- dplyr::mutate(data, chr = droplevels(chr))
    if (TRUE) {
      ## plot TF motif sequence logo.
      if (!requireNamespace("JASPAR2022", quietly = TRUE)) {
        BiocManager::install(
          c(
            "JASPAR2022", "DirichletMultinomial", "TFBSTools", 
            # "BSgenome.Hsapiens.UCSC.hg19", "EnsDb.Hsapiens.v75",
            "TxDb.Hsapiens.UCSC.hg19.knownGene"
          )
        )
      }
      tfs <- unique(data$TF_symbol)
      p.logos <- sapply(tfs, simplify = FALSE, 
        function(tf) {
          pfm <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022,
            opts = list(species = "Homo sapiens", name = tf)
          )
          TFBSTools::seqLogo(TFBSTools::toICM(pfm[[1]]))
          wrap(recordPlot(), 9, 3.5)
        })
      p.logos <- setLegend(p.logos, glue::glue("为 {tfs} motif sequence logo。{.seq_logo_legend()}"))
      x[[ glue::glue("p.logos_{name}") ]] <- .set_lab(p.logos, sig(x), names(p.logos), "motif sequence logo")
      x <- methodAdd(x, "{transmute(ref, '结合位点可视化')}以 R 包 `JASPAR2022` ({packageVersion('JASPAR2022')}) 和 `TFBSTools` ({packageVersion('TFBSTools')}) 绘制 {less(tfs)} Motif sequence logo。")
    }
    if (TRUE) {
      p.tracks <- sapply(seq_len(nrow(data)), simplify = FALSE, 
        function(n) {
          data <- data[n, ]
          range_bindingSite <- as_GRanges(
            data, col.start = "Start", col.end = "Stop"
          )
          range_bindingSite$id <- data$MatchSequence
          atrack <- Gviz::AnnotationTrack(range_bindingSite)
          btrack <- e(Gviz::BiomartGeneRegionTrack(
              genome = "hg19", name = "ENS.", symbol = data$target
              ))
          btrack@range <- btrack@range[btrack@range$symbol == data$target, ]
          otrack <- e(Gviz::OverlayTrack(trackList = list(btrack, atrack)))
          gtrack <- e(Gviz::GenomeAxisTrack())
          if (region) {
            atracks_region <- lapply(c("promoters", "transcripts"),
              function(type) {
                fun_get <- get_fun(type, envir = asNamespace("GenomicFeatures"))
                range <- fun_get(
                  TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                  columns = c("gene_id"), filter = list(gene_id = data$target_entrez)
                )
                range$id <- Hmisc::capitalize(type)
                Gviz::AnnotationTrack(range, shape = "box", name = "Region")
              }
            )
            otracks_type <- e(Gviz::OverlayTrack(trackList = atracks_region))
          }
          ranges <- lapply(otrack@trackList,
            function(x) {
              c(IRanges::start(x), IRanges::end(x))
            })
          lim <- e(grDevices::extendrange(range(unlist(ranges)), f = .5))
          if (region) {
            Gviz::plotTracks(
              list(gtrack, otracks_type, otrack),
              transcriptAnnotation = "symbol", sizes = c(1, 2, 2), 
              featureAnnotation = "id", fontcolor.feature = 1, 
              from = lim[1], to = lim[2]
            )
          } else {
            Gviz::plotTracks(
              list(gtrack, otrack),
              transcriptAnnotation = "symbol", sizes = c(1, 2), 
              featureAnnotation = "id", fontcolor.feature = 1, 
              from = lim[1], to = lim[2]
            )
          }
          wrap(recordPlot(), 11, 2)
        })
      names(p.tracks) <- paste0(data$TF_symbol, "_", data$target)
      p.tracks <- setLegend(p.tracks, glue::glue("为转录因子结合基因启动子示意 ({names(p.tracks)}) (Transcript 源于 BioMart 获取 ENSEMBL 数据库对应基因的注释区域；转录因子结合位点源于 Transcription Factor Target Gene Database)。 "))
      x[[ glue::glue("p.tracks_{name}") ]] <- .set_lab(
        p.tracks, sig(x), names(p.tracks), "transcript factor binding Illustrate"
      )
      x <- methodAdd(x, "以 R 包 `Gviz` ({packageVersion('Gviz')}) 绘制 ENSEMBL 数据库对应基因的转录子 (transcript)。")
    }
    x[[ glue::glue("map_data_{name}") ]] <- data
    x$.map_heading <- "TFBS 结合位点可视化"
    return(x)
  })

plot_region <- function(data, chr, genome = c("hg37", "hg38", "mm10", "rn6"),
  col.chr = "chr", col.start = "start", col.end = "end", zo = 1.2,
  ucscRefseq = FALSE, plot_pvalue = FALSE, pvalue = grpf(colnames(data), "pvalue", TRUE)[1],
  pvalueTrans = TRUE, symbol = NULL)
{
  if (is.numeric(chr)) {
    seq <- chr
    chr <- paste0("chr", chr)
  } else if (is.character(chr)) {
    seq <- strx(chr, "[0-9]+")
  }
  genome <- match.arg(genome)
  if (is(data, "GRanges")) {
    allRanges <- data[ data@seqnames == chr ]
  } else {
    if (missing(pvalue)) {
      message("Be Carefully, the columns of p-value is:\n\t", pvalue)
    }
    allRanges <- as_GRanges(data, chr, col.chr, col.start, col.end, pvalue = pvalue)
  }
  targetRanges <- allRanges
  ## Chromosome annotation Track
  if (FALSE) {
    iTrack <- e(Gviz::IdeogramTrack(genome = genome, chromosome = seq, name = paste0("Chromosome ", seq)))
  } else {
    iTrack <- NULL
  }
  gTrack <- e(Gviz::GenomeAxisTrack(col = "black", cex = 1, name = "Axis", fontcolor = "black"))
  if (is.null(iTrack)) {
    tracks <- list(gTrack)
    sizes <- 2
  } else {
    tracks <- list(iTrack, gTrack)
    sizes <- c(1, 2)
  }
  if (!is.null(symbol)) {
    fun_find <- function () {
      which <- which(allRanges$symbol == symbol)
      if (!length(which)) {
        stop("Nothing.")
      }
      message("Found targetRanges: ")
      r <- allRanges[ which ]
      print(r)
      r
    }
    targetRanges <- fun_find()
    if (missing(zo)) {
      zo <- 50
    }
    ucscRefseq <- TRUE
  }
  if (ucscRefseq) {
    message("Show UCSC annotation")
    print(targetRanges)
    ucscRange <- range.IRanges(targetRanges)
    if (!is.null(zo)) {
      ucscRange <- pmax(floor(zoRange(ucscRange, zo)), 0)
    }
    message("UCSC Track range: ", paste0(ucscRange, collapse = ", "))
    ucscTrack <- track_ncbiRefseq(genome, chr, ucscRange[1], ucscRange[2])
    tracks <- c(tracks, list(ucscTrack))
    sizes <- c(sizes, 3)
  }
  if (plot_pvalue) {
    if (pvalueTrans) {
      values <- allRanges$pvalue
      minLimit <- min(values[ values != 0 ]) / 10
      values[ values == 0 ] <- minLimit
      values <- -log10(values)
      allRanges$pvalue <- values
    }
    pvalueTrack <- e(Gviz::DataTrack(genome = genome, chromosome = chr,
      range = allRanges[, "pvalue"],
      name = if (pvalueTrans) "-Log10(P-value)" else "P-value",
      type = "hist", fill.histogram = "darkred", col.histogram = "transparent"
    ))
    tracks <- c(tracks, list(pvalueTrack))
    sizes <- c(sizes, if (small) 1 else 3)
  }
  ## plot
  ranges <- range.IRanges(targetRanges)
  if (is.null(zo)) {
    zo <- 1.2
  }
  message("In this ranges:")
  print(IRanges::subsetByOverlaps(allRanges, pseudo_ranges(ranges, chr)))
  e(Gviz::plotTracks(tracks, from = ranges[1], to = ranges[2], sizes = sizes))
  p.tracks <- recordPlot()
  wrap(p.tracks, 7, 6)
}


setMethod("map", signature = c(x = "job_tfbs", ref = "character"),
  function(x, ref, mode = c("allu", "shiny"), axes = 1:3, ...)
  {
    mode <- match.arg(mode)
    if (mode == "shiny") {
      stop('mode == "shiny"')
    } else if (mode == "allu") {
      data <- dplyr::select(x@tables$step1$res, target, TF_symbol, Motif, PValue)
      data <- dplyr::mutate(data, PValue = as.double(PValue))
      data <- split_lapply_rbind(data, seq_len(nrow(data)), args = list(fill = TRUE),
        function(x) {
          if (grpl(x$TF_symbol, "/")) {
            symbols <- strsplit(x$TF_symbol, "/")[[1]]
            x <- do.call(dplyr::bind_rows, rep(list(x), length(symbols)))
            x$TF_symbol <- symbols
          }
          x
        })
      p.venn <- new_venn(DB_TFs = data$TF_symbol, Candidates = ref)
      p.venn <- .set_lab(p.venn, sig(x), "intersection of TFs with queried Candidates")
      data <- dplyr::filter(data, TF_symbol %in% ref)
      p.allu <- new_allu(data, axes = axes, col.fill = 4, ...) +
        theme_minimal() +
        labs(fill = "P-value") +
        theme(axis.text.y = element_blank(),
          axis.title = element_blank()) +
        geom_blank()
      p.allu <- .set_lab(wrap(p.allu, 10, 7), sig(x), "The genes and related TFs")
      x$mapped <- namel(p.venn, p.allu, data)
    }
    return(x)
  })

.seq_logo_legend <- function() {
  "字母高度通常代表信息量（以bits为单位）或频率。信息量高意味着该位置对结合很重要，变异较少；频率则直接显示不同碱基的出现比例。字母排列顺序：高度高的字母在上方，表示该碱基在该位置更常见。总高度（Y轴）：可能代表该位置的信息量，信息量越高，表示该位置越保守，对结合越关键。"
}
