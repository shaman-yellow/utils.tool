# ==========================================================================
# workflow of mp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mp <- setClass("job_mp", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://bioconductor.org/packages/release/bioc/vignettes/MicrobiotaProcess/inst/doc/MicrobiotaProcess.html")
    ))

setGeneric("asjob_mp", 
  function(x, ...) standardGeneric("asjob_mp"))

setMethod("asjob_mp", signature = c(x = "job_qiime"),
  function(x){
    if (x@step < 6L) {
      stop("x@step < 6L")
    }
    object <- e(MicrobiotaProcess::mp_import_qiime2(
        x@params$table_qza, x@params$taxonomy,
        x@params$meta_file
        ))
    .job_mp(object = object)
  })

setMethod("step0", signature = c(x = "job_mp"),
  function(x){
    step_message("Prepare your data with function `asjob_mp`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_mp"),
  function(x){
    step_message("Alpha diversity analysis.
      "
    )
    if (is.null(x@params$mp_cal_rarecurve)) {
      object(x) <- e(MicrobiotaProcess::mp_rrarefy(object(x)))
      object(x) <- e(MicrobiotaProcess::mp_cal_rarecurve(object(x), .abundance = RareAbundance))
      x@params$mp_cal_rarecurve <- T
    }
    p1 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve)) +
      theme(legend.position = "none")
    p2 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve,
        .group = group, plot.group = T))
    p.rarefaction <- wrap(p1 + p2, 15, 5)
    if (is.null(x@params$mp_cal_alpha)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_alpha(object(x), .abundance = RareAbundance))
      x@params$mp_cal_alpha <- T
    }
    p.alpha_index <- e(MicrobiotaProcess::mp_plot_alpha(object(x), .group = group,
        .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
        ))
    p.alpha_index <- wrap(p.alpha_index, 15, 5)
    x@plots[[ 1 ]] <- namel(p.rarefaction, p.alpha_index)
    return(x)
  })

setMethod("step2", signature = c(x = "job_mp"),
  function(x){
    step_message("The visualization of taxonomy abundance")
    if (is.null(x@params$mp_cal_abundance)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = RareAbundance), text = "For samples")
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = RareAbundance, .group = group), text = "For groups")
      x@params$mp_cal_abundance <- T
    }
    p.rel_abundance_genus <- e(MicrobiotaProcess::mp_plot_abundance(
        object(x), .abundance = RareAbundance, .group = group,
        taxa.class = Genus, topn = 19, plot.group = T
        ), text = "For alluvium")
    p.rel_abundance_genus <- wrap(p.rel_abundance_genus, 12)
    p.heatmap_abundance_genus <- e(MicrobiotaProcess::mp_plot_abundance(
        object(x), .abundance = RareAbundance, .group = group,
        taxa.class = Genus, topn = 50, geom = 'heatmap',
        sample.dist = 'bray', sample.hclust = 'average'
        ), text = "For heatmap")
    p.heatmap_abundance_genus <- wrap(p.heatmap_abundance_genus,
      nrow(object(x)@colData) * .25, 15)
    x@plots[[ 2 ]] <- namel(p.rel_abundance_genus, p.heatmap_abundance_genus)
    return(x)
  })

setMethod("step3", signature = c(x = "job_mp"),
  function(x){
  })

setMethod("step4", signature = c(x = "job_mp"),
  function(x){})
