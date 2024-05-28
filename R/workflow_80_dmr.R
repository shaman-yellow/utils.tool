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
    method = "R package `Gviz` were used for methylation data visualization"
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
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

plot_dmr <- function(data, chr, genome = c("hg38", "mm10", "rn6"),
  col.chr = "chr", col.start = "start", col.end = "end",
  dmrRegion = F, cpgIsland = T, cpgIsland.overlap = T, subset.cpgIsland = NULL, zo = 1.2,
  dnaseHyper = T, ucscRefseq = F,
  plot_qvalue = T, plot_delta = T,
  use.delta = grpf(colnames(data), "_vs_")[1], use.qvalue = grpf(colnames(data), "qvalue", T)[1],
  qvalueTrans = T)
{
  if (missing(use.delta) || missing(use.qvalue)) {
    message("Be Carefully, the columns of delta and q-value is:\n\t", use.delta, ", ", use.qvalue)
  }
  if (is.numeric(chr)) {
    chr <- paste0("chr", chr)
  }
  genome <- match.arg(genome)
  if (is(data, "GRanges")) {
    dmrRanges <- data[ data@seqnames == chr ]
  } else {
    dmrRanges <- as_GRanges(data, chr, col.chr, col.start, col.end,
      delta = use.delta, qvalue = use.qvalue)
  }
  ## Chromosome annotation Track
  iTrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr, name = paste0("Chromosome ", chr))
  gTrack <- Gviz::GenomeAxisTrack(col = "black", cex = 1, name = "Axis", fontcolor = "black")
  tracks <- list(iTrack, gTrack)
  sizes <- c(1, 2)
  ## DMR region
  if (T) {
    ## for annotation
    if (dmrRegion) {
      aTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
        range = dmrRanges, name = "DMR")
      tracks <- c(tracks, list(aTrack))
      sizes <- c(sizes, 3)
    }
  }
  if (cpgIsland) {
    cpgIsland <- get_cpgisland_GRranges(genome)
    if (cpgIsland.overlap) {
      cpgIsland <- IRanges::subsetByOverlaps(cpgIsland, dmrRanges, type = "any")
    }
    islandTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
      range = cpgIsland, fill = "darkgreen",
      name = if (cpgIsland.overlap) "CpG Island overlap" else "CpG Island"
    )
    if (!is.null(subset.cpgIsland)) {
      if (is(subset.cpgIsland, "GRanges")) {
        RangeSubset <- range.IRanges(subset.cpgIsland)
      } else if (is(subset.cpgIsland, "numeric")) {
        RangeSubset <- cpgIsland[subset.cpgIsland, ]
        RangeSubset <- range.IRanges(RangeSubset)
      }
      RangeSubset <- GenomicRanges::GRanges(
        seqnames = S4Vectors::Rle(chr),
        range = IRanges::IRanges(RangeSubset[1], RangeSubset[2]),
        strand = S4Vectors::Rle("*")
      )
      dmrRanges <- IRanges::subsetByOverlaps(dmrRanges, RangeSubset, type = "any")
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
    ucscRange <- range.IRanges(dmrRanges)
    ucscTrack <- track_ncbiRefseq(genome, chr, ucscRange[1], ucscRange[2])
    tracks <- c(tracks, list(ucscTrack))
    sizes <- c(sizes, 3)
  }
  if (plot_qvalue) {
    if (qvalueTrans) {
      values <- dmrRanges$qvalue
      minLimit <- min(values[ values != 0 ]) / 10
      values[ values == 0 ] <- minLimit
      values <- -log10(values)
      dmrRanges$qvalue <- values
    }
    qvalueTrack <- Gviz::DataTrack(genome = genome, chromosome = chr,
      range = dmrRanges[, "qvalue"],
      name = if (qvalueTrans) "-Log10(Q-value)" else "Q-value",
      type = "hist", fill.histogram = "darkred", col.histogram = "transparent"
    )
    tracks <- c(tracks, list(qvalueTrack))
    sizes <- c(sizes, 3)
  }
  if (plot_delta) {
    deltaTrack <- Gviz::DataTrack(genome = genome, chromosome = chr,
      range = dmrRanges[, "delta"],
      name = "Delta", type = "gradient",
      gradient = RColorBrewer::brewer.pal(3, "RdBu")
    )
    tracks <- c(tracks, list(deltaTrack))
    sizes <- c(sizes, 3)
  }
  ## plot
  ranges <- range.IRanges(dmrRanges)
  if (is.null(zo)) {
    if (length(subset.cpgIsland) == 1) {
      zo <- 3
    } else {
      zo <- 1.2
    }
  }
  if (!is.null(subset.cpgIsland)) {
    ranges <- floor(zoRange(ranges, zo))
  }
  Gviz::plotTracks(tracks, from = ranges[1], to = ranges[2], sizes = sizes)
  p.tracks <- recordPlot()
  namel(p.tracks, cpgIsland)
}

track_ncbiRefseq <- function(genome, chr, from, to) {
  e(Gviz::UcscTrack(genome = "rn6", chromosome = chr, 
    track = "NCBI RefSeq", table = "ncbiRefSeq", 
    from = from, to = to,
    trackType = "GeneRegionTrack", 
    rstarts = "exonStarts", rends = "exonEnds", 
    gene = "name2", symbol = "name2", 
    transcript = "name", strand = "strand", 
    fill = "#8282d2", name = "NCBI RefSeq",
    showId = T, geneSymbol = TRUE
  ))
}

get_cpgisland_GRranges <- function(genome) {
  cpgIsland <- get_cpgIsland_data.ucsc(genome)
  cpgIsland <- dplyr::rename(cpgIsland, chr = chrom, start = chromStart, end = chromEnd)
  cpgIsland <- GenomicRanges::GRanges(
    seqnames = S4Vectors::Rle(cpgIsland$chr),
    ranges = IRanges::IRanges(start = cpgIsland$start, end = cpgIsland$end),
    strand = Rle(rep("*", nrow(cpgIsland)))
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
    strand = Rle(rep("*", nrow(data)))
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
