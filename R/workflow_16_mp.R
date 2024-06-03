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
    info = c("Tutorial: https://bioconductor.org/packages/release/bioc/vignettes/MicrobiotaProcess/inst/doc/MicrobiotaProcess.html"),
    cite = "[@MicrobiotaproceXuSh2023]",
    method = "R package `MicrobiotaProcess` used for microbiome data visualization",
    tag = "16s:qiime2+mp"
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
    error <- F
    plots <- list()
    if (relative(x)) {
      x$use <- "mpa_relativeAbundance"
      x$force <- T
    } else {
      x$use <- "RareAbundance"
      x$force <- F
      if (is.null(x@params$mp_cal_rarecurve)) {
        object(x) <- e(MicrobiotaProcess::mp_rrarefy(object(x)))
        res <- try(e(MicrobiotaProcess::mp_cal_rarecurve(object(x), .abundance = RareAbundance)), F)
        if (inherits(res, "try-error")) {
          error <- T
        } else {
          object(x) <- res
        }
        x@params$mp_cal_rarecurve <- T
      }
      if (!error) {
        p1 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve)) +
          theme(legend.position = "none")
        p2 <- e(MicrobiotaProcess::mp_plot_rarecurve(object(x), .rare = RareAbundanceRarecurve,
            .group = !!rlang::sym(x@params$group), plot.group = T))
        p.rarefaction <- wrap(p1 + p2, 15, 5)
        p.rarefaction <- .set_lab(p.rarefaction, sig(x), "alpha rarefaction")
      } else {
        p.rarefaction <- NULL
      }
      plots <- c(plots, namel(p.rarefaction))
    }
    if (is.null(x@params$mp_cal_alpha)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_alpha(object(x),
          .abundance = !!rlang::sym(x$use), force = x$force))
      x@params$mp_cal_alpha <- T
    }
    if (relative(x)) {
      use.alpha <- c("Observe", "Shannon", "Simpson", "Pielou")
      width <- 10
    } else {
      use.alpha <- c("Observe", "Chao1", "ACE", "Shannon", "Simpson", "Pielou")
      width <- 15
    }
    p.alpha_index <- e(MicrobiotaProcess::mp_plot_alpha(object(x),
        .group = !!rlang::sym(x@params$group),
        .alpha = !!use.alpha
        ))
    p.alpha_index <- wrap(p.alpha_index, width, 5)
    p.alpha_index <- .set_lab(p.alpha_index, sig(x), "alpha diversity")
    plots <- c(plots, namel(p.alpha_index))
    x@plots[[ 1 ]] <- plots
    return(x)
  })

setMethod("step2", signature = c(x = "job_mp"),
  function(x){
    step_message("The visualization of taxonomy abundance")
    if (is.null(x@params$ontology))
      x@params$ontology <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (is.null(x@params$mp_cal_abundance)) {
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = !!rlang::sym(x$use), force = x$force), text = "For samples")
      object(x) <- e(MicrobiotaProcess::mp_cal_abundance(object(x),
          .abundance = !!rlang::sym(x$use), force = x$force,
          .group = !!rlang::sym(x@params$group)),
        text = "For groups")
      x@params$mp_cal_abundance <- T
    }
    p.rel_abundance <- e(sapply(x@params$ontology, simplify = F,
        function(class) {
          p <- MicrobiotaProcess::mp_plot_abundance(
            object(x), .abundance = !!rlang::sym(x$use), force = x$force,
            .group = !!rlang::sym(x@params$group),
            taxa.class = !!rlang::sym(class), topn = 19, plot.group = T
          )
          if (class == "Species") width <- 16
          else width <- 14
          wrap(p, width)
        }), text = "For alluvium")
    p.rel_abundance <- .set_lab(p.rel_abundance, sig(x), names(p.rel_abundance), "abundance")
    if (FALSE) {
      p.heatmap_abundance <- e(sapply(x@params$ontology, simplify = F,
          function(class) {
            p <- MicrobiotaProcess::mp_plot_abundance(
              object(x), .abundance = !!rlang::sym(x$use), force = x$force,
              .group = !!rlang::sym(x@params$group),
              taxa.class = !!rlang::sym(class), topn = 50, geom = 'heatmap',
              sample.dist = 'bray', sample.hclust = 'average'
            )
            wrap(p, nrow(object(x)@colData) * .23 + 4, 12)
          }), text = "For heatmap")
    }
    x@plots[[ 2 ]] <- namel(p.rel_abundance)
    if (relative(x)) {
      x$abun.by.sample <- "mpa_relativeAbundanceBySample"
      x$rel.abun.by.sample <- "Relmpa_relativeAbundanceBySample"
    } else {
      x$abun.by.sample <- "RareAbundanceBySample"
      x$rel.abun.by.sample <- "RelRareAbundanceBySample"
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_mp"),
  function(x){
    step_message("Beta diversity analysis.")
    if (relative(x)) {
      use <- "mpa_relativeAbundance"
    } else {
      use <- "Abundance"
    }
    if (is.null(x@params$mp_decostand)) {
      object(x) <- e(MicrobiotaProcess::mp_decostand(object(x), .abundance = !!rlang::sym(use)))
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
    p.sample_dist <- try(wrap(e(MicrobiotaProcess::mp_plot_dist(object(x),
        .distmethod = bray, .group = !!rlang::sym(x@params$group)),
      text = "Heatmap")), T)
    p.sample_dist <- .set_lab(p.sample_dist, sig(x), "Sampe distance")
    p.group_test <- e(MicrobiotaProcess::mp_plot_dist(object(x), .distmethod = bray,
        .group = !!rlang::sym(x@params$group),
        group.test = TRUE, textsize = 3), text = "Boxplot")
    p.group_test <- set_palette(p.group_test)
    p.group_test <- wrap(p.group_test, 4)
    p.group_test <- .set_lab(p.group_test, sig(x), "Beta diversity group test")
    p.pcoa <- e(MicrobiotaProcess::mp_plot_ord(
        object(x), .ord = pcoa, .group = !!rlang::sym(x@params$group),
        .color = !!rlang::sym(x@params$group),
        .size = Observe, .alpha = Shannon, ellipse = TRUE))
    p.pcoa <- p.pcoa + ggrepel::geom_label_repel(aes(label = Sample), fill = "white", alpha = 1, size = 5)
    p.pcoa <- set_palette(p.pcoa)
    p.pcoa <- .set_lab(wrap(p.pcoa), sig(x), "PCoA")
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
          data <- tidyr::unnest(data, cols = !!rlang::sym(x$abun.by.sample))
          data <- dplyr::rename(data, Ontology = label)
          p <- p.hier + geom_fruit(data = data, geom = geom_col,
            mapping = aes(x = !!rlang::sym(x$rel.abun.by.sample), y = Sample, fill = Ontology),
            orientation = "y", pwidth = 3
          )
          p <- p + guides(fill = guide_legend(ncol = 1))
          width <- nrow(object(x)@colData) * .13 + 2
          wrap(as_grob(p), 12, if (width < 8) 8 else width)
        }))
    p.hier <- .set_lab(p.hier, sig(x), names(p.hier), "hierarchy")
    lab(p.hier) <- paste0(sig(x), " all hierarchy data")
    x@plots[[ 3 ]] <- namel(p.sample_dist, p.group_test, p.pcoa, p.hier)
    return(x)
  })

setMethod("step4", signature = c(x = "job_mp"),
  function(x, legend.color = F)
  {
    step_message("Biomarker discovery.")
    if (is.null(x@params$mp_diff_analysis)) {
      object(x) <- e(MicrobiotaProcess::mp_diff_analysis(
          object(x), .abundance = !!rlang::sym(x$use),
          .group = !!rlang::sym(x@params$group), force = x$force, filter.p = "fdr"))
      x@params$mp_diff_analysis <- T
    }
    ## plot radical tree
    p.tree <- try(e(MicrobiotaProcess::mp_plot_diff_cladogram(object(x),
        .group = !!rlang::sym(x@params$group),
        label.size = 2.5, hilight.alpha = .3, bg.tree.size = .5,
        bg.point.size = 2, bg.point.stroke = .25
        )), T)
    if (!inherits(p.tree, "try-error")) {
      p.tree <- p.tree + 
        MicrobiotaProcess::scale_fill_diff_cladogram(values = ggsci::pal_npg()(10)) +
        scale_size_continuous(range = c(1, 4))
      get_tree <- T 
    } else {
      get_tree <- F
    }
    ## extract results
    tree <- e(MicrobiotaProcess::mp_extract_tree(object(x), type = "taxatree"))
    sign <- paste0("Sign_", x@params$group)
    tryPvalue <- try(MicrobiotaProcess::select(tree, pvalue), T)
    if (inherits(tryPvalue, "try-error")) {
      message("\n\tNo significant features found.")
      Sys.sleep(3)
      tryPvalue <- F
    } else {
      top_table <- e(MicrobiotaProcess::select(tree, label, nodeClass,
          pvalue, fdr, LDAupper, LDAmean, LDAlower, !!rlang::sym(sign)))
      top_table <- dplyr::filter(top_table, !is.na(fdr))
      top_table <- dplyr::arrange(top_table, fdr)
      top_table <- .set_lab(top_table, sig(x), "statistic of all difference microbiota")
      x@tables[[ 4 ]] <- namel(top_table)
      if (get_tree) {
        if (legend.color) {
          tol <- length(unique(dplyr::filter(top_table, pvalue < .01)$label))
          if (tol) {
            h <- ceiling(sqrt(tol * 15))
            p.tree <- p.tree + guides(color = guide_legend(nrow = h))
          } else {
            h <- 0
          }
          print(p.tree)
          p.tree <- recordPlot()
          p.tree <- wrap(p.tree, 12 + h * .12, 7 + h * .06)
        } else {
          p.tree <- p.tree + guides(color = "none")
          print(p.tree)
          p.tree <- recordPlot()
          p.tree <- wrap(p.tree, 12, 11)
        }
        p.tree <- .set_lab(p.tree, sig(x), "The cladogram of differential species")
      }
      tryPvalue <- T
    }
    ## plot boxplot
    if (tryPvalue) {
      classes <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
      fun_test <- function(classes) {
        max <- 50
        data <- dplyr::filter(top_table, pvalue < .01)
        stats <- table(data$nodeClass)
        stats <- stats[ match(classes, names(stats)) ]
        classes[ cumsum(stats) <= max ]
      }
      classes <- fun_test(classes)
      p.box <- e(MicrobiotaProcess::mp_plot_diff_boxplot(object(x),
          taxa.class = !!classes,
          .group = !!rlang::sym(x@params$group)))
      p.box <- MicrobiotaProcess::set_diff_boxplot_color(p.box,
        values = ggsci::pal_npg()(10),
        guide = guide_legend(title = NULL)
      )
      p.box <- wrap(p.box, 10, 10)
      p.box <- .set_lab(p.box, sig(x), "The abundance and LDA ",
        paste0("from ", classes[1], " to ", tail(classes, n = 1)))
    } else {
      p.box <- NULL
    }
    x@plots[[ 4 ]] <- namel(p.box, p.tree)
    return(x)
  })

setMethod("step5", signature = c(x = "job_mp"),
  function(x, match,
    use = c("pvalue", "fdr"), cutoff = .05, classes = c("Species", "Genus", "Family"))
  {
    step_message("From Microbiota to metabolites.")
    ## prepare patterns to match
    use <- match.arg(use)
    mic_tops <- filter(mp@tables$step4$top_table, !!rlang::sym(use) < cutoff,
      nodeClass %in% dplyr::all_of(classes))
    mt.pattern <- .create_mt_pattern(mic_tops$label)
    gm <- job_gutmd()
    .add_internal_job(gm)
    db_data <- object(gm)
    ## find related metabolites in database
    matched <- matchThats(db_data[[1]], mt.pattern)
    db_data_matched <- filter(db_data,
      Gut.Microbiota %in% dplyr::all_of(unique(unlist(matched))),
      Substrate != "")
    require(ggalluvial)
    axes <- c("Gut.Microbiota", "Substrate", "Metabolite")
    data <- dplyr::select(db_data_matched, dplyr::all_of(axes))
    data <- dplyr::mutate(data,
      match = ifelse(Substrate %in% !!match | Metabolite %in% !!match, T, F),
      Match = ifelse(match, "Sig.", "unSig."),
    )
    data <- to_lodes_form(data, key = "Types", axes = 1:3)
    freq <- table(data$stratum)
    label.notshow <- names(freq)[ as.integer(freq) <= fivenum(as.integer(freq))[4] ]
    fun <- function(x) {
      ifelse(x %in% label.notshow, "", as.character(x))
    }
    data <- dplyr::mutate(data, label = fun(stratum),
      label = ifelse(match, as.character(stratum), label)
    )
    aes <- aes(x = Types, y = 1, label = label,
      stratum = stratum, alluvium = alluvium)
    p.alluvial <- ggplot(data, aes) +
      geom_alluvium(aes(fill = Match)) +
      geom_stratum(fill = "lightyellow") +
      geom_text(stat = "stratum") +
      labs(fill = "Metabolites", y = "") +
      scale_fill_manual(values = ggsci::pal_npg()(2)) +
      theme_minimal() +
      theme(axis.text.y = element_blank())
    x$matched <- matched
    x@tables[[ 5 ]] <- namel(db_data, db_data_matched)
    x@plots[[ 5 ]] <- namel(p.alluvial)
    return(x)
  })

setMethod("pattern", signature = c(x = "job_mp"),
  function(x, use = c("pvalue", "fdr"), cutoff = .05, not = c("GB", "^un$", "^bacterium$"))
  {
    if (x@step < 4L) {
      stop("x@step < 4L")
    }
    data <- x@tables$step4$top_table
    if (is.null(data) || !nrow(data)) {
      stop("No significant microbiota found.")
    }
    use <- match.arg(use)
    data <- dplyr::filter(data, !!rlang::sym(use) < !!cutoff)
    lst <- data$label
    # onto <- strx(lst, "^[a-z]")
    lst <- gs(lst, "^.__([^_]+).*", "\\1")
    # names(lst) <- onto
    lst <- lst[ !grpl(lst, paste0(not, collapse = "|")) ]
    lst[ !duplicated(lst) ]
  })

.create_mt_pattern <- function(strings) {
  mt.pattern <- stringr::str_extract_all(strings, "(?<=__|^).*?(?=__|$)")
  mt.pattern <- unique(unlist(mt.pattern))
  mt.pattern <- mt.pattern[nchar(mt.pattern) > 1]
  mt.pattern <- mt.pattern[ mt.pattern != "bacterium" ]
  mt.pattern
}

matchThats <- function(x, patterns, ignore.case = T) {
  patterns <- patterns[ !is.na(patterns) ]
  if (length(patterns)) {
    matched <- .find_and_sort_strings(unique(x), patterns, ignore.case = ignore.case)
    matched <- nl(patterns, matched)
    matched <- lst_clear0(matched)
    matched
  } else {
    list()
  }
}

setMethod("relative", signature = c(x = "job_mp"),
  function(x){
    if (is.null(x$relative)) {
      F
    } else T
  })
