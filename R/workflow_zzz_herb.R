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
    compounds.keep.intersection = F,
    syns = NULL, enrichment = NULL, en.top = 10, use.enrich = c("kegg", "go"), ...)
  {
    message("Filter compounds targets with disease targets.")
    data <- x$data.allu
    if (!is.null(compounds)) {
        dat.compounds <- dplyr::filter(data, Ingredient.name %in% !!compounds)
      if (compounds.keep.intersection) {
        data <- dplyr::filter(data, Target.name %in% !!dat.compounds$Target.name)
      } else {
        data <- dat.compounds
      }
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
    if (!is.null(enrichment)) {
      res <- plot_network.enrich(data, enrichment, use.enrich, en.top, HLs = HLs, ...)
      p.pharm <- res$p.pharm
      p.pharm$.path <- res$dataPath
    } else {
      p.pharm <- plot_network.pharm(data, HLs = HLs, ax2.level = levels,
        lab.fill = lab.level, force.ax1 = herbs, ...)
    }
    p.pharm$.data <- data
    if (length(ref)) {
      x[[ paste0("p.pharm2", name) ]] <- .set_lab(p.pharm, sig(x), "network pharmacology with disease")
    } else {
      x$p.pharmMap <-  .set_lab(p.pharm, sig(x), "network pharmacology")
    }
    return(x)
  })


plot_network.enrich <- function(data, enrichment, use.enrich = c("kegg", "go"), en.top = 10,
  use.sub = c("BP", "CC", "MF"), ...)
{
  if (is(enrichment, "job_enrich")) {
    use.enrich <- match.arg(use.enrich)
    enrichment <- enrichment@tables$step1[[ paste0("res.", use.enrich) ]][[ 1 ]]
    if (use.enrich == "go") {
      enrichment <- dplyr::filter(enrichment, ont == !!use.sub)
    }
    message("Use fisrt enrichment results of ", use.enrich, ": ", use.sub)
  }
  en <- en.sig <- head(enrichment, en.top)
  en <- dplyr::select(en, Description, geneName_list)
  en <- reframe_col(en, "geneName_list", unlist)
  ax.names <- colnames(data)
  data <- tbmerge(data, en, by.x = ax.names[3], by.y = "geneName_list", all.x = T, allow.cartesian = T)
  data <- dplyr::relocate(data, !!!rlang::syms(ax.names))
  en.sig <- dplyr::select(en.sig, Description, p.adjust)
  ## plot
  p.pharm <- plot_network.pharm(
    data, lab.fill = "P.adjust",
    ax4 = "Pathway", ax4.level = en.sig,
    decImport.ax4Level = T, ...
  )
  dataPath <- reframe_col(data, "Description",
    function(x) paste0(length(x), "###", paste0(x, collapse = "; ")))
  dataPath <- tidyr::separate(dataPath, Description, c("Hit_pathway_number", "Enriched_pathways"), "###")
  dataPath <- dplyr::mutate(dataPath, Hit_pathway_number = as.integer(Hit_pathway_number))
  dataPath <- dplyr::arrange(dataPath, dplyr::desc(Hit_pathway_number))
  namel(p.pharm, dataPath)
}

filter_hob_dl <- function(ids, filter.hob, filter.dl, test = F, use.cids = T) {
  if (filter.hob || filter.dl) {
    if (test) {
      ids <- head(ids, n = 100)
    }
    if (use.cids) {
      if (!is.numeric(ids)) {
        stop("is.numeric(ids) == F")
      }
      smiles <- get_smiles_batch(ids, n = 50)
      smiles <- dplyr::rename(smiles, ID = CID, SMILES = IsomericSMILES)
    } else {
      message("use.cids == F, so `ids` should be SMILES with names.")
      smiles <- ids
      message("De duplicated smiles.")
      smiles <- smiles[ !duplicated(smiles) ]
      if (any(duplicated(names(smiles)))) {
        stop("any(duplicated(names(smiles)))")
      }
      smiles <- tibble::tibble(ID = names(smiles), SMILES = unname(smiles))
    }
    message("Number of unique compounds: ", nrow(smiles))
    setSmiles <- list(All_compounds = smiles$SMILES)
  }
  if (filter.hob) {
    hob <- job_hob(smiles)
    hob <- step1(hob)
    hobSmiles <- dplyr::filter(res(hob), isOK)$smiles
    .add_internal_job(hob, T)
    setSmiles <- c(setSmiles, list(HOB_is_OK = hobSmiles))
  }
  if (filter.dl) {
    dl <- job_dl(smiles)
    dl <- step1(dl)
    dl <- step2(dl)
    dlSmiles <- dplyr::filter(res(dl), isOK)$smiles
    .add_internal_job(dl, T)
    setSmiles <- c(setSmiles, list(DL_is_OK = dlSmiles))
  }
  if (filter.hob || filter.dl) {
    p.upset <- new_upset(lst = setSmiles)
    ids <- dplyr::filter(smiles, SMILES %in% !!p.upset$ins)$ID
    p.upset$ins <- smiles$ID[ match(p.upset$ins, smiles$SMILES) ]
    namel(ids, p.upset)
  } else {
    return()
  }
}

plot_network.pharm <- function(data, f.f = 2.5, f.f.mul = .7, f.f.sin = .2, f.ax4 = 2, seed = sample(1:10, 1), HLs = NULL,
  ax1 = "Herb", ax2 = "Compound", ax3 = "Target", ax4 = NULL, less.label = T,
  ax2.level = NULL, ax4.level = NULL, decImport.ax4Level = F, lab.fill = "",
  edge_width = .1, force.ax1 = NULL, ax3.layout = c("spiral", "grid"), ax4.layout = c("circle", "linear"),
  HLs.label = T, ...)
{
  if (length(unique(data[[1]])) == 1) {
    sherb <- 1L
  } else {
    sherb <- 0L
  }
  spiral <- F
  spiral_order <- NULL
  if (!is.null(ax2.level)) {
    if (ncol(ax2.level) > 2) {
      message("Too many columns found in `ax2.level`, try select `name` and `level`")
      ax2.level <- dplyr::select(ax2.level, name, level)
    }
    if (!identical(colnames(ax2.level), c("name", "level"))) {
      message("Rename colnames of `ax2.level` with 'name', 'level'")
      colnames(ax2.level) <- c("name", "level")
    }
  }
  ax3.layout <- match.arg(ax3.layout)
  ax4.layout <- match.arg(ax4.layout)
  prepare_data <- function(data) {
    colnames(data) <- c(ax1, ax2, ax3, ax4)
    nodes <- tidyr::gather(data, type, value)
    nodes <- dplyr::distinct(nodes)
    nodes <- dplyr::relocate(nodes, name = value, type)
    nodes <- dplyr::filter(nodes, !is.na(name))
    nodes <- dplyr::mutate(nodes, name = ifelse(duplicated(name), paste0(type, ":", name), name))
    ed.12 <- dplyr::distinct(data, dplyr::pick(1:2))
    if (sherb) {
      ComMul <- ed.12[[ ax2 ]]
    } else {
      freq12 <- table(ed.12[[ ax2 ]])
      ComSin <- names(freq12[ freq12 == 1 ])
      ComMul <- names(freq12[ freq12 > 1 ])
      ed.12sin <- dplyr::filter(ed.12, !!rlang::sym(ax2) %in% !!ComSin)
    }
    resize <- function(data, f) {
      data$x <- data$x * f
      data$y <- data$y * f
      data
    }
    shift <- function(data, x = 0, y = 0) {
      data$x <- data$x + x
      data$y <- data$y + y
      data
    }
    ## axis 3 coords
    nodesTgt <- dplyr::filter(nodes, type == !!ax3)
    # if (ax3.layout == "spiral" && nrow(nodesTgt) < 30) {
      # if (ax4.layout != "circle") {
      #   message("Nodes of 'ax3' too small, set layout as: grid")
      #   ax3.layout <- "grid"
      # }
    # }
    if (ax3.layout == "grid") {
      crds.Tgt <- get_layout(NULL, "grid", nodes = nodesTgt)
      crds.Tgt <- shift(crds.Tgt, -max(crds.Tgt$x) / 2, -max(crds.Tgt$y) / 2)
    } else if (ax3.layout == "spiral") {
      crds.Tgt <- get_coords.spiral(nrow(nodesTgt), nCir = 2, minRad = 2)
      crds.Tgt <- get_layout(NULL, crds.Tgt, nodes = nodesTgt)
      f.f <- f.f * .7
    }
    f.rsz <- max(crds.Tgt$x) * f.f
    if (!is.null(ax4)) {
      crds.Tgt$y <- crds.Tgt$y * 0.6
      if (is.null(ax4.level)) {
        crds.ax4 <- get_layout(NULL, ax4.layout, nodes = dplyr::filter(nodes, type == !!ax4))
      } else {
        colnames(ax4.level) <- c("name", "level")
        if (any(duplicated(ax4.level))) {
          message("`ax4.level` is duplicated, de duplicated herein.")
          ax4.level <- dplyr::distinct(ax4.level, name, .keep_all = T)
        }
        ## change for parent.frame
        ax4.level <<- ax4.level
        nodesAx4 <- dplyr::arrange(ax4.level, if (decImport.ax4Level) dplyr::desc(level) else level)
        crds.ax4 <- get_layout(NULL, ax4.layout, nodes = nodesAx4)
      }
      if (ax4.layout == "linear") {
        crds.ax4 <- shift(crds.ax4, -max(crds.ax4$x) / 2, -max(crds.Tgt$y) * f.ax4)
        f.rsz <- f.rsz * f.ax4 / 2
      } else if (ax4.layout == "circle") {
        if (!f.rsz) {
          f.rsz <- 3 * f.f
        }
        # f.rsz / f.f, e.g., max(crds.Tgt$x)
        crds.ax4 <- resize(crds.ax4, f.rsz / f.f * .3)
      }
    }
    if (!f.rsz) {
      f.rsz <- 3 * f.f
    } 
    if (!sherb) {
      if (!is.null(force.ax1)) {
        force.ax1 <- force.ax1[ !force.ax1 %in% dplyr::filter(nodes, type == !!ax1)$name ]
        nodes <- tibble::add_row(nodes, type = ax1, name = force.ax1)
      }
      crds.Hrb <- get_layout(NULL, "circle", nodes = dplyr::filter(nodes, type == !!ax1))
      crds.Hrb <- resize(crds.Hrb, f.rsz)
      lst <- split(ed.12sin, ed.12sin[[ ax1 ]])
      # coords for compounds from separate herb
      crds.ComSin <- mapply(lst, names(lst), SIMPLIFY = F,
        FUN = function(ed, nm) {
          lay <- resize(get_layout(ed, "star"), f.rsz * f.f.sin)
          crd.hrb <- dplyr::filter(crds.Hrb, name == !!nm)
          dplyr::filter(shift(lay, crd.hrb$x, crd.hrb$y), name != !!nm)
        })
      crds.ComSin <- do.call(dplyr::bind_rows, crds.ComSin)
    } else {
      crds.Hrb <- NULL
      crds.ComSin <- NULL
    }
    ## deside spiral
    crds.ComMul <- get_layout(NULL, "circle", nodes = data.frame(name = ComMul))
    if (sherb & length(ComMul) > 50) {
      ## coords for spiral nodes
      if (!is.null(ax2.level)) {
        ax2.level <- dplyr::arrange(ax2.level, dplyr::desc(level))
        if (!all(isIn <- crds.ComMul$name %in% ax2.level$name)) {
          message("!all(crds.ComMul %in% ax2.level$name), fill with NA")
          ax2.level <- tibble::add_row(ax2.level,
            name = crds.ComMul$name[ !isIn ])
        }
        ax2.level <- dplyr::filter(ax2.level, name %in% !!crds.ComMul$name)
        ## reorder by level
        crds.ComMul <- crds.ComMul[ match(ax2.level$name, crds.ComMul$name), ]
      }
      crds.ComMul[, c("x", "y")] <- get_coords.spiral(nrow(crds.ComMul))
      spiral <<- T
      spiral_order <<- crds.ComMul$name
    }
    crds.ComMul <- resize(crds.ComMul, f.rsz * f.f.mul)
    if (!is.null(ax4)) {
      crds <- lst_clear0(list(crds.Hrb, crds.ComSin, crds.ComMul, crds.Tgt, crds.ax4))
    } else {
      crds <- lst_clear0(list(crds.Hrb, crds.ComSin, crds.ComMul, crds.Tgt))
    }
    crds <- do.call(dplyr::bind_rows, crds)
    crds <- dplyr::distinct(crds, name, .keep_all = T)
    if (!is.null(ax4)) {
      if (sherb) {
        use.cols <- 2:3
      } else {
        use.cols <- 1:3
      }
    } else {
      if (sherb) {
        use.cols <- 2
      } else {
        use.cols <- 1:2
      }
    }
    edges <- lapply(use.cols,
      function(n) {
        edges <- dplyr::distinct(data, dplyr::pick(n:(n+1)))
        colnames(edges) <- c("source", "target")
        edges
      })
    edges <- do.call(dplyr::bind_rows, edges)
    edges <- dplyr::filter(edges, !is.na(target))
    nodes <- dplyr::mutate(nodes, .type = type)
    if (!sherb) {
      nodes <- dplyr::mutate(nodes,
        type = ifelse(name %in% !!ComSin, paste0(!!ax2, " Unique"),
          ifelse(name %in% !!ComMul, paste0(!!ax2, " Sharing"), type)))
    }
    nodes <- dplyr::filter(nodes, name %in% !!crds$name)
    crds <- crds[ match(nodes$name, crds$name), ]
    layout <- dplyr::select(crds, x, y)
    if (!is.null(HLs)) {
      edges <- dplyr::mutate(edges,
        highlight = ifelse(source %in% HLs & target %in% HLs, "Highlight", "Non-highlight")
      )
    }
    message("Generating network layout.")
    fast_layout(edges, layout = layout, nodes = nodes)
  }
  x <- prepare_data(data)
  data <- as_tibble(data.frame(x))
  data <- dplyr::mutate(
    data, cent = ifelse(type == !!ax1, cent * 1000 + 2000,
      ifelse(grpl(type, !!ax2, ignore.case = T), cent * 20, cent)))
  data.tgt <- dplyr::filter(data, type == !!ax3)
  set.seed(seed)
  minSize <- .5
  edge_data <- get_edges("short")(x)
  if (!is.null(ax4)) {
    edge_ax4 <- dplyr::filter(edge_data, node2.type == !!ax4)
    edge_data <- dplyr::filter(edge_data, node2.type != !!ax4)
    geom_edge_ax4 <- ggraph::geom_edge_link(data = edge_ax4,
      # aes(edge_color = highlight, edge_width = highlight),
      edge_color = sample(color_set()[1:14], 1), alpha = .5)
  } else {
    geom_edge_ax4 <- geom_blank()
  }
  if (is.null(HLs)) {
    if (missing(edge_width)) {
      if (nrow(dplyr::filter(data, type == !!ax2)) > 10L) {
        width <- .1
      } else {
        width <- 1
        minSize <- 3
      }
    } else {
      width <- edge_width
    }
    geom_edge <- ggraph::geom_edge_link(data = edge_data,
      color = sample(color_set()[1:6], 1), alpha = .2, width = width)
  } else {
    geom_edge <- geom_edge_link(data = edge_data,
      aes(edge_color = highlight, edge_width = highlight), alpha = .2)
  }
  if (is.null(ax4)) {
    dataAx123 <- data
    geom_label_3 <- geom_blank()
  } else {
    dataAx123 <- dplyr::filter(data, type != !!ax4)
    if (ax4.layout == "linear") {
      geom_label_3 <- ggplot2::geom_text(
        data = dplyr::filter(data, type == !!ax4),
        aes(x = x, y = y * 1.1, label = name, color = type),
        angle = 45, hjust = 1
      )
    } else if (ax4.layout == "circle") {
      geom_label_3 <- ggrepel::geom_label_repel(
        data = dplyr::filter(data, type == !!ax4),
        aes(x = x, y = y, label = name, color = type)
      )
    }
  }
  if (less.label) {
    geom_label_1 <- ggrepel::geom_label_repel(
      data = dplyr::filter(dataAx123, cent > if (sherb) 100 else 1000),
      aes(x = x, y = y, label = name, color = type))
    geom_label_2 <- ggrepel::geom_label_repel(
      data = dplyr::filter(data.tgt, cent >= sort(cent, T)[10]),
      aes(x = x, y = y, label = name, color = type))
  } else {
    geom_label_1 <- ggrepel::geom_label_repel(
      data = dataAx123, aes(x = x, y = y, label = name, color = type))
    geom_label_2 <- geom_blank()
  }
  if (is.null(ax2.level)) {
    geom_node_1 <- geom_node_point(data = dataAx123,
      aes(x = x, y = y, color = type, size = cent, shape = type))
    geom_node_2 <- geom_blank()
    scale_fill_gradient <- geom_blank()
  } else {
    data.level <- dplyr::filter(dataAx123, .type == !!ax2)
    data.level <- map(data.level, "name", ax2.level, "name", "level", col = "level")
    geom_node_1 <- geom_node_point(data = data.level,
      aes(x = x, y = y, fill = level, size = cent), shape = 21,
      stroke = 0, color = "transparent"
    )
    data.nonlevel <- dplyr::filter(dataAx123, .type != !!ax2)
    geom_node_2 <- geom_node_point(data = data.nonlevel,
      aes(x = x, y = y, color = type, size = cent, shape = type))
    scale_fill_gradient <- scale_fill_gradientn(colors = rev(color_set2()))
  }
  sizeRange <- c(minSize, 15)
  if (is.null(ax4.level)) {
    geom_node_3 <- geom_blank()
  } else {
    dataNode3 <- dplyr::filter(data, type == !!ax4)
    dataNode3 <- map(dataNode3, "name", ax4.level, "name", "level", col = "level")
    ## must sort the data, else the `size` were mapped error
    dataNode3 <- dplyr::arrange(dataNode3, x)
    sizeNode3 <- seq(sizeRange[1], sizeRange[2] / 2, length.out = length(unique(dataNode3$name)))
    ## the order of the order
    dataNode3 <- dplyr::mutate(dataNode3, size = sizeNode3[ order(order(level, decreasing = decImport.ax4Level)) ])
    geom_node_3 <- geom_node_point(
      data = dataNode3, aes(fill = level),
      size = dataNode3$size, shape = 21,
      stroke = 0, color = "transparent"
    )
  }
  p <- ggraph(x) + geom_edge + geom_edge_ax4 +
    geom_node_1 + geom_node_2 + geom_node_3 +
    geom_label_1 + geom_label_2 + geom_label_3 +
    scale_size(range = sizeRange) +
    guides(size = "none", shape = "none") +
    scale_color_manual(values = rstyle("pal", seed)) +
    scale_edge_color_manual(values = c("Highlight" = "red", "Non-highlight" = "lightblue")) +
    scale_edge_width_manual(values = c("Highlight" = edge_width * 5, "Non-highlight" = edge_width)) +
    scale_fill_gradient +
    theme_minimal() +
    labs(color = "Type", edge_color = "Highlight", edge_width = "Highlight", fill = lab.fill) +
    theme(axis.text = element_blank(),
      axis.title = element_blank()) +
    geom_blank()
  if (!is.null(HLs)) {
    data <- dplyr::filter(data, name %in% HLs)
    p <- p + geom_point(data = data, aes(x = x, y = y), shape = 21, color = "red", size = 10)
    if (HLs.label) {
      p <- p + ggrepel::geom_label_repel(data = data, aes(x = x, y = y, label = name),
        size = 7, color = "black")
    }
  }
  if (sherb) {
    if (spiral) {
      p <- wrap(p, 15, 12)
      p$spiral_order <- spiral_order
      p
    } else {
      wrap(p, 10, 8)
    }
  } else {
    wrap(p, 15, 12)
  }
}
