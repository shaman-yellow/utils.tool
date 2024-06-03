# ==========================================================================
# workflow of dmr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_dmr <- setClass("job_dmr", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "dmr",
    info = c("https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html"),
    cite = "[@VisualizingGenHahne2016]",
    method = "R package `Gviz` were used for methylation data visualization",
    tag = "methyl:dmr"
    ))

job_dmr <- function(data, genome = c("hg38", "mm10", "rn6"),
  col.chr = "chr", col.start = "start", col.end = "end", col.gene = "symbol",
  use.delta = grpf(colnames(data), "_vs_")[1],
  use.qvalue = grpf(colnames(data), "qvalue", T)[1])
{
  data <- dplyr::filter(data, !is.na(!!rlang::sym(col.chr)))
  object <- as_GRanges(data, NULL, col.chr, col.start, col.end,
    delta = use.delta, qvalue = use.qvalue, symbol = col.gene)
  x <- .job_dmr(object = object)
  x$data <- data
  x$genome <- genome
  x
}

setMethod("step0", signature = c(x = "job_dmr"),
  function(x){
    step_message("Prepare your data with function `job_dmr`.")
  })

setMethod("step1", signature = c(x = "job_dmr"),
  function(x, plot_cpgIsland = T)
  {
    step_message("DMR distribution.")
    data <- dplyr::filter(as.data.frame(object(x)), !is.na(symbol))
    p.volcano <- plot_volcano(data, label = "symbol",
      use = "qvalue", label.p = "Q-value",
      use.fc = "delta", label.fc = "Delta",
    )
    p.volcano <- wrap(p.volcano, 5, 4)
    p.volcano <- .set_lab(p.volcano, sig(x), "All DMR volcano plot")
    p.distri <- new_pie(object(x)@seqnames, title = "DMR distribution")
    p.distri <- wrap(p.distri, 7, 5)
    p.distri <- .set_lab(p.distri, sig(x), "DMR distribution")
    subset <- object(x)[!is.na(object(x)$symbol)]
    p.distriGene <- new_pie(subset@seqnames, title = "DMR in Genes")
    p.distriGene <- wrap(p.distriGene, 7, 5)
    p.distriGene <- .set_lab(p.distriGene, sig(x), "DMR in Genes")
    ## CpG Island
    cpgIsland <- get_cpgisland_GRranges(x$genome)
    cpgIsland <- x$cpgIsland <- IRanges::subsetByOverlaps(cpgIsland, object(x))
    if (T) {
      ## stats overlap
      function_stats <- function() {
        metaOverlaps <- IRanges::findOverlaps(object(x), cpgIsland)
        obj <- object(x)[ metaOverlaps@from ]
        obj <- split(obj, as.character(obj@seqnames))
        res <- lapply(obj,
          function(x) {
            ifelse(is.na(x$symbol), "CpG Island", "CpG Island + Gene")
          })
        as_df.lst(res)
      }
      dat <- function_stats()
      p.viewOverlap <- ggplot(dat, aes(x = sortCh(type, T, T), y = 1, fill = name)) +
        ggchicklet::geom_chicklet(width = 0.7, radius = grid::unit(5, "pt")) + 
        scale_fill_brewer(palette = "Set2") +
        labs(x = "Chromosome", y = "Number", fill = "Type") +
        theme_minimal() +
        coord_flip()
      p.viewOverlap <- wrap(p.viewOverlap, 6, colSum(dat$type) * .3 + .5)
      p.viewOverlap <- .set_lab(p.viewOverlap, sig(x), "specific methylation")
    }
    if (plot_cpgIsland) {
      p.cpgIsland <- plot_multiChr(cpgIsland, x$genome, "cpgIsland")
      p.cpgIsland <- .set_lab(wrap(p.cpgIsland), sig(x), "CpG Island methylation")
    } else {
      p.cpgIsland <- NULL
    }
    #################################
    x@plots[[ 1 ]] <- namel(p.volcano, p.distri, p.viewOverlap, p.distriGene, p.cpgIsland)
    return(x)
  })

setMethod("step2", signature = c(x = "job_dmr"),
  function(x, chrs = NULL, geneOverlap = T, must = T)
  {
    step_message("Specify chromosomes to visualization.")
    if (is.null(chrs)) {
      message("Use chromosome which have CpG Island overlap.")
      chrs <- unique(as.character(x$cpgIsland@seqnames))
    }
    obj <- object(x)
    if (geneOverlap) {
      message("Plot DMRs overlap genes")
      note <- "DMRs overlap genes visualization"
      obj <- obj[ which(!is.na(obj$symbol)) ]
    } else {
      note <- "all DMRs visualization"
    }
    p.dmrs <- sapply(chrs, simplify = F,
      function(chr) {
        # if must is T, maybe return NULL
        p <- plot_dmr(obj, chr, x$genome, mustOverlapCpgIsland = must, gr.cpgIsland = x$cpgIsland)
        if (!is.null(p)) {
          wrap(p, 7, 4)
        }
      })
    if (length(p.dmrs)) {
      p.dmrs <- lst_clear0(p.dmrs)
      p.dmrs <- .set_lab(p.dmrs, sig(x), names(p.dmrs), note)
      lab(p.dmrs) <- "all DMRs"
      x@plots[[ 2 ]] <- namel(p.dmrs)
    }
    return(x)
  })

setMethod("vis", signature = c(x = "job_dmr"),
  function(x, chr, subset = NULL, ucscRefseq = F, zo = 1.2, symbol = NULL, ...)
  {
    if (length(symbol) > 1) {
      stop("`symbol` should be a length 1 character.")
    }
    if (!is.null(symbol) && missing(chr)) {
      chr <- as.character(object(x)[ which(object(x)$symbol == symbol) ]@seqnames)
      message("Found chr of ", symbol, " in ", chr)
    }
    if (is.numeric(chr)) {
      chr <- paste0("chr", chr)
    }
    p <- plot_dmr(object(x), chr, x$genome, subset.cpgIsland = subset,
      ucscRefseq = ucscRefseq, zo = zo, gr.cpgIsland = x$cpgIsland, symbol = symbol, ...)
    p <- .set_lab(p, sig(x), chr, paste0(symbol, " DMR annotation"))
    p
  })

plot_dmr <- function(data, chr, genome = c("hg38", "mm10", "rn6"),
  col.chr = "chr", col.start = "start", col.end = "end",
  dmrRegion = F, cpgIsland = T, cpgIsland.overlap = T, subset.cpgIsland = NULL, zo = 1.2,
  dnaseHyper = T, ucscRefseq = F,
  plot_qvalue = T, plot_delta = T,
  use.delta = grpf(colnames(data), "_vs_")[1], use.qvalue = grpf(colnames(data), "qvalue", T)[1],
  qvalueTrans = T, mustOverlapCpgIsland = F, gr.cpgIsland = NULL, symbol = NULL)
{
  if (is.numeric(chr)) {
    chr <- paste0("chr", chr)
  }
  genome <- match.arg(genome)
  if (is(data, "GRanges")) {
    allRanges <- data[ data@seqnames == chr ]
  } else {
    if (missing(use.delta) || missing(use.qvalue)) {
      message("Be Carefully, the columns of delta and q-value is:\n\t", use.delta, ", ", use.qvalue)
    }
    allRanges <- as_GRanges(data, chr, col.chr, col.start, col.end,
      delta = use.delta, qvalue = use.qvalue)
  }
  targetRanges <- allRanges
  ## Chromosome annotation Track
  iTrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr, name = paste0("Chromosome ", chr))
  gTrack <- Gviz::GenomeAxisTrack(col = "black", cex = 1, name = "Axis", fontcolor = "black")
  tracks <- list(iTrack, gTrack)
  sizes <- c(1, 2)
  if (!is.null(symbol) || !is.null(subset.cpgIsland)) {
    small <- T
  } else small <- F
  ## DMR region
  if (T) {
    ## for annotation
    if (dmrRegion) {
      aTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
        range = allRanges, name = "DMR")
      tracks <- c(tracks, list(aTrack))
      sizes <- c(sizes, 3)
    }
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
    ucscRefseq <- T
  }
  if (cpgIsland) {
    if (is.null(gr.cpgIsland)) {
      cpgIsland <- get_cpgisland_GRranges(genome)
    } else {
      cpgIsland <- gr.cpgIsland
    }
    if (cpgIsland.overlap) {
      cpgIsland <- IRanges::subsetByOverlaps(cpgIsland, allRanges, type = "any")
    }
    if (mustOverlapCpgIsland && !length(cpgIsland)) {
      message("No overlap of cpgIsland of ", chr)
      return()
    } else {
      message(chr, ": overlap number is:", length(cpgIsland))
    }
    islandTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
      range = cpgIsland, fill = "darkgreen",
      name = if (cpgIsland.overlap) "CpG Island overlap" else "CpG Island"
    )
    if (!is.null(symbol)) {
      message("Use ranges of `symbol` location")
    } else if (!is.null(subset.cpgIsland)) {
      ## just found dmr area
      if (is(subset.cpgIsland, "GRanges")) {
        RangeSubset <- range.IRanges(subset.cpgIsland)
      } else if (is(subset.cpgIsland, "numeric")) {
        RangeSubset <- cpgIsland[subset.cpgIsland, ]
        RangeSubset <- range.IRanges(RangeSubset)
      }
      RangeSubset <- pseudo_ranges(RangeSubset, chr)
      targetRanges <- IRanges::subsetByOverlaps(allRanges, RangeSubset, type = "any")
    }
    tracks <- c(tracks, list(islandTrack))
    sizes <- c(sizes, 1)
  }
  if (dnaseHyper && genome == "hg38") {
    stop("This function of code has not finished.\n",
      "See: https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#obtaining-the-data",
      "\n\t", "And `get_dnase_data`")
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
  if (plot_qvalue) {
    if (qvalueTrans) {
      values <- allRanges$qvalue
      minLimit <- min(values[ values != 0 ]) / 10
      values[ values == 0 ] <- minLimit
      values <- -log10(values)
      allRanges$qvalue <- values
    }
    qvalueTrack <- Gviz::DataTrack(genome = genome, chromosome = chr,
      range = allRanges[, "qvalue"],
      name = if (qvalueTrans) "-Log10(Q-value)" else "Q-value",
      type = "hist", fill.histogram = "darkred", col.histogram = "transparent"
    )
    tracks <- c(tracks, list(qvalueTrack))
    sizes <- c(sizes, if (small) 1 else 3)
  }
  if (plot_delta) {
    deltaTrack <- Gviz::DataTrack(genome = genome, chromosome = chr,
      range = allRanges[, "delta"],
      name = "Delta", type = "gradient",
      gradient = rev(RColorBrewer::brewer.pal(3, "RdBu")),
      ylim = c(-1, 1)
    )
    tracks <- c(tracks, list(deltaTrack))
    sizes <- c(sizes, if (small) 1 else 3)
  }
  ## plot
  ranges <- range.IRanges(targetRanges)
  if (is.null(zo)) {
    if (length(subset.cpgIsland) == 1) {
      zo <- 50
    } else {
      zo <- 1.2
    }
  }
  if (!is.null(subset.cpgIsland) || !is.null(symbol)) {
    ranges <- floor(zoRange(ranges, zo))
  }
  message("DMR in this ranges:")
  print(IRanges::subsetByOverlaps(allRanges, pseudo_ranges(ranges, chr)))
  Gviz::plotTracks(tracks, from = ranges[1], to = ranges[2], sizes = sizes)
  p.tracks <- recordPlot()
  wrap(p.tracks, 7, 6)
}

track_ncbiRefseq <- function(genome, chr, from, to, db_dir = .prefix("ucsc_tracks", "db")) {
  dir.create(db_dir, F)
  file <- file.path(db_dir, paste0("ncbiRefSeq_", genome, "_", chr, "_", from, "_", to, ".rds"))
  if (file.exists(file)) {
    readRDS(file)
  } else {
    track <- e(Gviz::UcscTrack(genome = genome, chromosome = chr, 
        from = from, to = to,
        track = "NCBI RefSeq", table = "ncbiRefSeq", 
        trackType = "GeneRegionTrack", 
        rstarts = "exonStarts", rends = "exonEnds", 
        gene = "name2", symbol = "name2", 
        transcript = "name", strand = "strand", 
        fill = "#8282d2", name = "NCBI RefSeq",
        showId = T, geneSymbol = TRUE
        ))
    saveRDS(track, file)
    track
  }
}

get_cpgisland_GRranges <- function(genome) {
  cpgIsland <- get_cpgIsland_data.ucsc(genome)
  cpgIsland <- dplyr::rename(cpgIsland, chr = chrom, start = chromStart, end = chromEnd)
  cpgIsland <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(cpgIsland$chr),
    ranges = IRanges::IRanges(start = cpgIsland$start, end = cpgIsland$end),
    strand = S4Vectors::Rle(rep("*", nrow(cpgIsland)))
  )
  cpgIsland
}

range.IRanges <- function(x) {
  c(min(IRanges::start(x@ranges)), max(IRanges::end(x@ranges)))
}

as_GRanges <- function(data, chr = NULL,
  col.chr = "chr", col.start = "start", col.end = "end", ...)
{
  if (!is.null(chr)) {
    data <- dplyr::filter(data, !!rlang::sym(col.chr) == !!chr)
  }
  if (!nrow(data)) {
    stop("Nothing to be convert.")
  }
  data <- dplyr::mutate(data, dplyr::across(dplyr::all_of(c(col.start, col.end)), as.integer))
  argsMeta <- list(
    seqnames = S4Vectors::Rle(data[[ col.chr ]]),
    ranges = IRanges::IRanges(start = data[[ col.start ]], end = data[[ col.end ]]),
    strand = S4Vectors::Rle(rep("*", nrow(data)))
  )
  metas <- c(col.chr, col.start, col.end)
  data <- dplyr::select(data, -dplyr::all_of(metas))
  rename <- list(...)
  if (length(rename)) {
    data <- dplyr::rename(data, !!!rlang::syms(rename))
  }
  argsOthers <- as.list(data)
  do.call(GenomicRanges::GRanges, c(argsMeta, argsOthers))
}

plot_multiChr <- function(x, genome, name = NULL, zo = 1.2) {
  stopifnot(is(x, "GRanges"))
  dat <- data.frame(chr = sortCh(unique(as.character(x@seqnames)), F, T))
  lst <- split(x, as.character(x@seqnames))
  ranges <- lapply(lst, function(x) floor(zoRange(range.IRanges(x), zo)))
  scale <- signif(min(vapply(ranges, diff, FUN.VALUE = double(1))) / 2, 1)
  p <- lattice::xyplot(1 ~ chr | chr, dat,
    scales = list(draw = FALSE), xlab = NULL, ylab = NULL,
    panel = function(x) {
      x <- as.character(x)
      itrack <- Gviz::IdeogramTrack(genome = genome, chromosome = x)
      gtrack <- Gviz::GenomeAxisTrack(scale = scale, labelPos = "below", exponent = 3)
      atrack <- Gviz::AnnotationTrack(range = lst[[ x ]], genome = genome, chromosome = x, name = name)
      Gviz::plotTracks(list(itrack, atrack, gtrack),
        chromosome = x, add = TRUE, showId = FALSE,
        from = ranges[[ x ]][1], to = ranges[[ x ]][2]
      )
    })
  print(p)
  p <- recordPlot(p)
  p
}

sortCh <- function(x, decreasing = F, as.factor = F) {
  num <- as.integer(strx(x, "[0-9]+"))
  x <- x[ order(num, decreasing = decreasing) ]
  if ( as.factor ) {
    x <- factor(x, levels = unique(x))
  }
  x
}

pseudo_ranges <- function(range, chr, strand = "*") {
  GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(chr),
    range = IRanges::IRanges(range[1], range[2]),
    strand = S4Vectors::Rle(strand)
  )
}
