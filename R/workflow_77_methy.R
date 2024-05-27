# ==========================================================================
# workflow of methy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_methy <- setClass("job_methy", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "methy",
    info = c("https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html",
      "https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html",
      "https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html#43_DataTrack"
      ),
    cite = "[@MissmethylAnPhipso2016; @MinfiAFlexibAryee2014; @VisualizingGenHahne2016; @CallingDifferePeters2021]",
    method = "R package `missMethyl`, `minfi`, `Gviz`, `DMRcate` were used for analysing methylation array data"
    ))

job_methy <- function()
{
  .job_methy()
}

setMethod("step0", signature = c(x = "job_methy"),
  function(x){
    step_message("Prepare your data with function `job_methy`.")
  })

setMethod("step1", signature = c(x = "job_methy"),
  function(x){
    step_message("")
    return(x)
  })

plot_dmr <- function(data, chr, genome = c("hg38", "mm10", "rn6"),
  col.chr = "chr", col.start = "start", col.end = "end",
  dmrRegion = F, cpgIsland = T, cpgIsland.overlap = T,
  dnaseHyper = T, ucscRefseq = T,
  plot_qvalue = T, plot_delta = T,
  use.qvalue = "DMR_Qvalue", use.delta = "DMR_Treatment_vs_Model")
{
  if (is.numeric(chr)) {
    chr <- paste0("chr", chr)
  }
  genome <- match.arg(genome)
  data <- dplyr::filter(data, !!rlang::sym(col.chr) == !!chr)
  if (!nrow(data)) {
    stop("Nothing to be plot.")
  }
  data <- dplyr::mutate(data, dplyr::across(dplyr::all_of(c(col.start, col.end)), as.integer))
  ## Chromosome annotation Track
  iTrack <- Gviz::IdeogramTrack(genome = genome, chromosome = chr,
    name = paste0("Chromosome ", chr))
  gTrack <- Gviz::GenomeAxisTrack(col = "black", cex = 1, name = "Axis", fontcolor = "black")
  tracks <- list(iTrack, gTrack)
  ## DMR region
  if (T) {
    ## for annotation
    dmrRanges <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(data[[ col.chr ]]),
      ranges = IRanges::IRanges(start = data[[ col.start ]], end = data[[ col.end ]]),
      strand = Rle(rep("*", nrow(data)))
    )
    aTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
      range = dmrRanges, name = "DMR")
    if (dmrRegion) {
      tracks <- c(tracks, list(aTrack))
    }
  }
  if (cpgIsland) {
    if (T) {
      cpgIsland <- get_cpgIsland_data.ucsc(genome)
      cpgIsland <- dplyr::rename(cpgIsland, chr = chrom, start = chromStart, end = chromEnd)
    } else {
      # cpgIsland <- dplyr::filter(cpgIsland, chr == !!chr)
    }
    cpgIsland <- GenomicRanges::GRanges(
      seqnames = S4Vectors::Rle(cpgIsland$chr),
      ranges = IRanges::IRanges(start = cpgIsland$start, end = cpgIsland$end),
      strand = Rle(rep("*", nrow(cpgIsland)))
    )
    if (cpgIsland.overlap) {
      cpgIsland <- IRanges::subsetByOverlaps(cpgIsland, dmrRanges, type = "any")
    }
    islandTrack <- Gviz::AnnotationTrack(genome = genome, chromosome = chr,
      range = cpgIsland, fill = "darkgreen", name = "CpG Island"
    )
    tracks <- c(tracks, list(islandTrack))
  }
  if (dnaseHyper && genome == "hg38") {
    stop("This function of code has not finished.\n",
      "See: https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html#obtaining-the-data",
      "\n\t", "And `get_dnase_data`")
  }
  if (ucscRefseq && nrow(data) == 1) {
    message("Show UCSC annotation")
    # Gviz::UcscTrack
    trackName <- switch(strx(genome, "[a-z]+"),
      "hg" = "Human mRNAs", "rn" = "Rat mRNAs", "mm" = "Mouse mRNAs"
    )
    stop_debug(namel(dmrRanges))
    tcscTrack <- e(Gviz::UcscTrack(
      genome = genome, chromosome = chr, track = trackName,
      from = IRanges::start(dmrRanges@ranges), to = IRanges::end(dmrRanges@ranges),
      trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name",
      symbol = "name2", transcript = "name", strand = "strand",
      fill = "darkblue", stacking = "squish", name = "mRNA",
      showId = TRUE, geneSymbol = TRUE
    ))
    tracks <- c(tracks, list(tcscTrack))
  }
  if (plot_qvalue) {
    
  }
  ## plot
  Gviz::plotTracks(tracks,
    from = min(data[[ col.start ]]),
    to = max(data[[ col.end ]])
  )
}
