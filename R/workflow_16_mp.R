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
  function(x, group = "group"){
    step_message("Alpha diversity analysis.
      "
    )
    x@params$group <- group
    if (is.null(x@params$mp_cal_rarecurve)) {
      object(x) <- e(MicrobiotaProcess::mp_rrarefy(object(x)))
      object(x) <- e(MicrobiotaProcess::mp_cal_rarecurve(object(x), .abundance = RareAbundance))
      x@params$mp_cal_rarecurve <- T
    }
    p1 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve)) +
      theme(legend.position = "none")
    p2 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve,
        .group = !!rlang::sym(x@params$group), plot.group = T))
    p.rarefaction <- wrap(p1 + p2, 15, 5)
    if (is.null(x@params$mp_cal_alpha)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_alpha(object(x), .abundance = RareAbundance))
      x@params$mp_cal_alpha <- T
    }
    p.alpha_index <- e(MicrobiotaProcess::mp_plot_alpha(object(x),
        .group = !!rlang::sym(x@params$group),
        .alpha = c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
        ))
    p.alpha_index <- wrap(p.alpha_index, 15, 5)
    x@plots[[ 1 ]] <- namel(p.rarefaction, p.alpha_index)
    return(x)
  })

setMethod("step2", signature = c(x = "job_mp"),
  function(x){
    step_message("The visualization of taxonomy abundance")
    if (is.null(x@params$ontology))
      x@params$ontology <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (is.null(x@params$mp_cal_abundance)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = RareAbundance), text = "For samples")
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = RareAbundance, .group = !!rlang::sym(x@params$group)),
        text = "For groups")
      x@params$mp_cal_abundance <- T
    }
    p.rel_abundance <- e(sapply(x@params$ontology, simplify = F,
        function(class) {
          p <- MicrobiotaProcess::mp_plot_abundance(
            object(x), .abundance = RareAbundance,
            .group = !!rlang::sym(x@params$group),
            taxa.class = !!rlang::sym(class), topn = 19, plot.group = T
          )
          if (class == "Species") width <- 16
          else width <- 14
          wrap(p, width)
        }), text = "For alluvium")
    if (FALSE) {
      p.heatmap_abundance <- e(sapply(x@params$ontology, simplify = F,
          function(class) {
            p <- MicrobiotaProcess::mp_plot_abundance(
              object(x), .abundance = RareAbundance,
              .group = !!rlang::sym(x@params$group),
              taxa.class = !!rlang::sym(class), topn = 50, geom = 'heatmap',
              sample.dist = 'bray', sample.hclust = 'average'
            )
            wrap(p, nrow(object(x)@colData) * .23 + 4, 12)
          }), text = "For heatmap")
    }
    x@plots[[ 2 ]] <- namel(p.rel_abundance)
    return(x)
  })

setMethod("step3", signature = c(x = "job_mp"),
  function(x){
    step_message("Beta diversity analysis.")
    if (is.null(x@params$mp_decostand)) {
      object(x) <- e(MicrobiotaProcess::mp_decostand(object(x), .abundance = Abundance))
      object(x) <- e(MicrobiotaProcess::mp_cal_dist(object(x), .abundance = hellinger,
          distmethod = "bray"))
      x@params$mp_decostand <- T
    }
    if (is.null(x@params$mp_cal_pcoa)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_pcoa(object(x), .abundance = hellinger,
          distmethod = "bray"))
      object(x) <- e(MicrobiotaProcess::mp_adonis(object(x),
          .abundance = hellinger, .formula = as.formula(paste0("~ ", x@params$group)),
          distmethod = "bray", action = "add"))
      x@params$mp_cal_pcoa <- T
    }
    if (is.null(x@params$mp_cal_clust)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_clust(
          object(x), .abundance = hellinger, distmethod = "bray",
          hclustmethod = "average", action = "add"
          ))
      x@params$mp_cal_clust <- T
    }
    p.sample_dist <- e(MicrobiotaProcess::mp_plot_dist(object(x),
        .distmethod = bray, .group = !!rlang::sym(x@params$group)),
      text = "Heatmap")
    p.group_test <- e(MicrobiotaProcess::mp_plot_dist(object(x), .distmethod = bray,
        .group = !!rlang::sym(x@params$group),
        group.test = TRUE, textsize = 3), text = "Boxplot")
    p.group_test <- set_palette(p.group_test)
    p.group_test <- wrap(p.group_test, 4)
    p.pcoa <- e(MicrobiotaProcess::mp_plot_ord(
        object(x), .ord = pcoa, .group = !!rlang::sym(x@params$group),
        .color = !!rlang::sym(x@params$group),
        .size = Observe, .alpha = Shannon, ellipse = TRUE))
    p.pcoa <- set_palette(p.pcoa)
    ## hierarchy
    data_hier <- e(MicrobiotaProcess::mp_extract_internal_attr(object(x), name = 'SampleClust'))
    require(ggtree)
    require(ggtreeExtra)
    p.hier <- ggtree(data_hier) + 
      geom_tippoint(aes(color = !!rlang::sym(x@params$group))) +
      geom_tiplab(as_ylab = TRUE) +
      scale_x_continuous(expand = c(0, 0.01))
    p.hier <- set_palette(p.hier)
    ## class
    p.hier <- e(sapply(x@params$ontology, simplify = F,
        function(class) {
          data <- e(MicrobiotaProcess::mp_extract_abundance(object(x),
              taxa.class = !!rlang::sym(class), topn = 30))
          data <- tidyr::unnest(data, cols = RareAbundanceBySample)
          data <- dplyr::rename(data, Ontology = label)
          p <- p.hier + geom_fruit(data = data, geom = geom_col,
            mapping = aes(x = RelRareAbundanceBySample, y = Sample, fill = Ontology),
            orientation = "y", pwidth = 3
          )
          p <- p + guides(fill = guide_legend(ncol = 1))
          wrap(p, 12, nrow(object(x)@colData) * .13 + 2)
        }))
    x@plots[[ 3 ]] <- namel(p.sample_dist, p.group_test, p.pcoa, p.hier)
    return(x)
  })

setMethod("step4", signature = c(x = "job_mp"),
  function(x){
    step_message("Biomarker discovery.")
    if (is.null(x@params$mp_diff_analysis)) {
      object(x) <- e(MicrobiotaProcess::mp_diff_analysis(
          object(x), .abundance = RelRareAbundanceBySample,
          .group = !!rlang::sym(x@params$group)))
      x@params$mp_diff_analysis <- T
    }
    ## plot boxplot
    p.box <- e(MicrobiotaProcess::mp_plot_diff_boxplot(object(x),
        taxa.class = c(Phylum, Class, Order, Family, Genus, Species),
        .group = !!rlang::sym(x@params$group)))
    p.box <- MicrobiotaProcess::set_diff_boxplot_color(p.box,
      values = ggsci::pal_npg()(10),
      guide = guide_legend(title = NULL)
    )
    ## plot radical tree
    p.tree <- e(MicrobiotaProcess::mp_plot_diff_cladogram(object(x),
        label.size = 2.5, hilight.alpha = .3, bg.tree.size = .5,
        bg.point.size = 2, bg.point.stroke = .25))
    p.tree <- p.tree + 
      MicrobiotaProcess::scale_fill_diff_cladogram(values = ggsci::pal_npg()(10)) +
      scale_size_continuous(range = c(1, 4))
    x@plots[[ 4 ]] <- namel(p.box, p.tree)
    ## extract results
    tree <- e(MicrobiotaProcess::mp_extract_tree(object(x), type = "taxatree"))
    sign <- paste0("Sign_", x@params$group)
    top_table <- e(MicrobiotaProcess::select(tree, label, nodeClass,
        pvalue, fdr, LDAupper, LDAmean, LDAlower, !!rlang::sym(sign)))
    top_table <- dplyr::filter(top_table, !is.na(fdr))
    top_table <- dplyr::arrange(top_table, fdr)
    x@tables[[ 4 ]] <- namel(top_table)
    return(x)
  })