# ==========================================================================
# workflow of herb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClassUnion("JOB_herb", c("job_herb", "job_tcmsp", "job_batman", "job_tcmsp2"))

setMethod("slice", signature = c(x = "JOB_herb"),
  function(x, ...){
    x@params$herbs_info <- dplyr::slice(x@params$herbs_info, ...)
    return(x)
  })

setMethod("vis", signature = c(x = "JOB_herb"),
  function(x, col.fill = 2, axes = 1:4, label.auto = F, label.freq = NULL, label.factor = 1)
  {
    symDisease <- x@tables$step3$disease_targets_annotation$hgnc_symbol
    data <- x@params$data.allu
    data <- dplyr::filter(data, Target.name %in% !!symDisease)
    data <- dplyr::mutate(data, Disease =  x$disease)
    p <- new_allu(data, col.fill, axes = axes,
      label.auto = label.auto, label.freq = label.freq,
      label.factor = label.factor)
    p <- .set_lab(p, sig(x), "Targets of compounds and related disease" )
    p
  })

setMethod("intersect", signature = c(x = "JOB_herb", y = "JOB_herb"),
  function(x, y, names)
  {
    lst <- list(unique(x$data.allu$Target.name),
      unique(y$data.allu$Target.name))
    names(lst) <- names
    venn <- new_venn(lst = lst)
    venn <- .set_lab(venn, "intersected targets of ", paste0(names, collapse = " and "))
    if (length(venn$ins)) {
      lst <- list(x, y)
      names(lst) <- names
      lst <- lapply(lst,
        function(obj) {
          dplyr::filter(obj$data.allu, Target.name %in% !!venn$ins)
        })
      lst <- frbind(lst, idcol = "From")
      p.pharm <- plot_network.pharm(lst[, -1])
      list(p.venn = venn, p.pharm = p.pharm, data = lst)
    } else {
      message("No intersection.")
      venn
    }
  })

setMethod("map", signature = c(x = "JOB_herb", ref = "list"),
  function(x, ref, HLs = NULL, levels = NULL, lab.level = "Level", name = "dis", compounds = NULL,
    syns = NULL, ...)
  {
    message("Filter compounds targets with disease targets.")
    data <- x$data.allu
    if (!is.null(compounds)) {
      data <- dplyr::filter(data, Ingredient.name %in% !!compounds)
    }
    if (!is.null(syns)) {
      data <- map(data, colnames(data)[2], syns, colnames(syns)[1], colnames(syns)[2], rename = F)
    }
    if (length(ref)) {
      data <- dplyr::filter(data, Target.name %in% !!unlist(ref, use.names = F))
      p.venn2dis <- new_venn(Diseases = unlist(ref, use.names = F), Targets = rm.no(x$data.allu$Target.name))
      x[[ paste0("p.venn2", name) ]] <- .set_lab(p.venn2dis, sig(x), "Targets intersect with targets of diseases")
    }
    herbs <- unique(x$data.allu[[1]])
    p.pharm <- plot_network.pharm(data, HLs = HLs, ax2.level = levels,
      lab.fill = lab.level, force.ax1 = herbs, ...)
    p.pharm$.data <- data
    if (length(ref)) {
      x[[ paste0("p.pharm2", name) ]] <- .set_lab(p.pharm, sig(x), "network pharmacology with disease")
    } else {
      x$p.pharmMap <-  .set_lab(p.pharm, sig(x), "network pharmacology")
    }
    return(x)
  })



