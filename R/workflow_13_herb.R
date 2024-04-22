# ==========================================================================
# workflow of herb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_herb <- setClass("job_herb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@HerbAHighThFang2021]",
    method = "Website `HERB` <http://herb.ac.cn/> used for TCM data source"
    ))

job_herb <- function(herbs, db = get_herb_data())
{
  x <- .job_herb(object = db)
  x@params$herbs <- herbs
  herbs_info <- filter(object(x)$herb, Herb_cn_name %in% !!params(x)$herbs)
  herbs_info <- .set_lab(herbs_info, "herbs information")
  x@params$herbs_info <- herbs_info
  print(herbs_info)
  message("Got the herbs:\n\n\t", paste0(herbs_info$Herb_cn_name, collapse = ", "),
      "\n\n", "\tTotal: ", colSum(herbs_info$Herb_cn_name))
  message("By modify `x@params$herbs_info` to filter the results.")
  x
}

setMethod("step0", signature = c(x = "job_herb"),
  function(x){
    step_message("Prepare your data with function `job_herb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_herb"),
  function(x, filter.hob = T, filter.dl = T, ...,
    tempdir = "download", db = .prefix("herb/herbs_ingredient.rds", "db"))
  {
    step_message("Dowload compounds of herbs.")
    ids <- params(x)$herbs_info$Herb_
    db <- new_db(db, "herb_id")
    db <- not(db, ids)
    if (length(db@query)) {
      x$tempdir <- tempdir
      link <- start_drive(browser = "firefox", download.dir = x$tempdir)
      Sys.sleep(3)
      tryCatch(
        download_herbCompounds(link, db@query, tempdir = x$tempdir),
        finally = end_drive()
      )
      fun_clear <- function() {
        if (file.exists("herbs_ingredient")) {
          unlink("herbs_ingredient", T, T)
        }
      }
      fun_clear()
      herbs_compounds <- collateFiles_herbs(db@query, from = x$tempdir, to = "herbs_ingredient")
      db <- upd(db, herbs_compounds, db@query)
      fun_clear()
      unlink(tempdir, T, T)
    }
    herbs_compounds <- dplyr::filter(db@db, herb_id %in% ids)
    herbs_compounds <- .set_lab(herbs_compounds, sig(x), "Components of Herbs")
    if (filter.hob || filter.dl) {
      smiles <- dplyr::filter(object(x)$component,
        Ingredient_id %in% !!unique(herbs_compounds$Ingredient.id),
        Ingredient_Smile != "Not Available")
      smiles <- dplyr::select(smiles, Ingredient_id, Ingredient_Smile)
      smiles <- nl(smiles$Ingredient_id, smiles$Ingredient_Smile, F)
      ## Check smiles available
      message("Use rcdk to check smiles validity.")
      notValide <- vapply(rcdk::parse.smiles(unname(smiles)), is.null, logical(1))
      if (any(notValide)) {
        message("Some smiles was not valid: ", length(which(notValide)))
      }
      smiles <- smiles[ !notValide ]
      lst <- filter_hob_dl(smiles, filter.hob, filter.dl, use.cids = F, ...)
      if (length(lst)) {
        p.upset <- .set_lab(lst$p.upset, sig(x), "Prediction of HOB and Drug-Likeness")
        x@plots[[ 1 ]] <- namel(p.upset)
        herbs_compounds <- dplyr::filter(herbs_compounds, Ingredient.id %in% !!lst$ids)
      }
    }
    x@tables[[ 1 ]] <- namel(herbs_compounds)
    return(x)
  })

setMethod("step2", signature = c(x = "job_herb"),
  function(x, group_number = 50, group_sleep = 3,
    searchDir = .prefix(), db = .prefix("herb/compounds_target.rds", "db"),
    removeExistsTemp = F)
  {
    step_message("Dowload targets of compounds")
    all_ids <- unique(x@tables$step1$herbs_compounds$Ingredient.id)
    db <- new_db(db, "Ingredient_id")
    db <- not(db, all_ids)
    if (length(db@query)) {
      queries <- db@query
      merge_data <- F
      if (!is.null(searchDir)) {
        message("The feature of `searchDir` was deprecated.")
        dbfile <- list.files(searchDir, "^HBIN.*.xlsx$", full.names = T, recursive = T)
        names(dbfile) <- get_realname(dbfile)
        dbfile <- dbfile[ !duplicated(dbfile) ]
        has_ids <- queries[ queries %in% names(dbfile) ]
        if (length(has_ids)) {
          files <- dbfile[ names(dbfile) %in% has_ids ]
          data <- collateFiles_herbs(readFrom = files, .id = "Ingredient_id", from = x$tempdir)
          queries <- queries[ !queries %in% has_ids ]
          merge_data <- T
        }
      }
      if (length(queries)) {
        queries <- grouping_vec2list(queries, group_number, T)
        link <- start_drive(browser = "firefox")
        Sys.sleep(3)
        lst <- pbapply::pblapply(queries,
          function(ids) {
            .dir <- "compounds_target"
            to <- paste0(.dir, "/", attr(ids, "name"))
            if (file.exists(to)) {
              .dir <- timeName("compounds_target")
              to <- paste0(.dir, "/", attr(ids, "name"))
            }
            download_compoundTargets(link, ids, tempdir = x$tempdir)
            data <- collateFiles_herbs(ids, to = to, .id = "Ingredient_id", from = x$tempdir)
            cli::cli_alert_info("Sleep between groups ...")
            Sys.sleep(group_sleep)
            data
          })
        lst <- data.table::rbindlist(lst, fill = T)
        end_drive()
      } else {
        lst <- data.frame()
      }
      if (merge_data) {
        if (!is(data, "list")) {
          data <- list(data)
        }
        if (nrow(lst)) {
          lst <- c(list(lst), data)
        } else {
          lst <- data
        }
      }
      if (is(lst, "list")) {
        data <- data.table::rbindlist(lst, fill = T)
      } else {
        data <- lst
      }
      compounds_targets <- as_tibble(data)
      db <- upd(db, compounds_targets, db@query)
    }
    compounds_targets <- dplyr::filter(db@db, Ingredient_id %in% !!all_ids)
    x@tables[[ 2 ]] <- namel(compounds_targets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_herb"),
  function(x, disease = NULL, disease.score = 5,
    mart_dataset = "hsapiens_gene_ensembl", HLs = NULL)
  {
    step_message("Get disease targets from genecards and annotate targets with biomaRt.
      As well, plot the intersection of herbs targets.
      "
    )
    if (is.null(x@params$mart)) {
      mart <- new_biomart(mart_dataset)
      x@params$mart <- mart
    } else {
      mart <- x@params$mart
    }
    ## annotate compounds targets
    message("Annotate compounds targets")
    compounds_targets <- x@tables$step2$compounds_targets
    targets <- unique(compounds_targets$Target.name)
    targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol", targets)
    x@params$targets_annotation <- targets_annotation
    compounds_targets_annotation <- tbmerge(
      compounds_targets, targets_annotation,
      by.x = "Target.name", by.y = "hgnc_symbol", all.x = T,
      allow.cartesian = T
    )
    ## create data: herbs_target
    message("Create data: herbs_target")
    herbs_compounds <- x@tables$step1$herbs_compounds
    if (!identical(class(herbs_compounds$Ingredient.id), class(compounds_targets$Ingredient_id))) {
      message("'ID' not identical in class, convert as character.")
      herbs_compounds$Ingredient.id %<>% as.character()
      compounds_targets$Ingredient_id %<>% as.character()
    }
    herbs_targets <- tbmerge(
      herbs_compounds, compounds_targets,
      by.x = "Ingredient.id", by.y = "Ingredient_id", all.x = T,
      allow.cartesian = T
    )
    x@params$herbs_targets <- herbs_targets
    easyRead <- map(x@params$herbs_targets, "herb_id", object(x)$herb, "Herb_", "Herb_pinyin_name")
    easyRead <- .set_lab(easyRead, sig(x), "Herbs compounds and targets")
    x@params$easyRead <- easyRead
    if (length(unique(herbs_targets$herb_id)) > 1) {
      ## plot upset of herbs targets
      sets <- dplyr::distinct(herbs_targets, herb_id, Target.name)
      sets <- split(sets$Target.name, sets$herb_id)
      fun <- function(sets) x@params$herbs_info$Herb_pinyin_name[match(names(sets), x@params$herbs_info$Herb_)]
      names(sets) <- fun(sets)
      p.herbs_targets <- new_upset(lst = sets)
      p.herbs_targets <- .set_lab(p.herbs_targets, sig(x), "Intersection of herbs all targets")
      ## plot upset of herbs ingredient
      sets <- split(herbs_compounds$Ingredient.id, herbs_compounds$herb_id)
      names(sets) <- fun(sets)
      p.herbs_compounds <- new_upset(lst = sets)
      p.herbs_compounds <- .set_lab(p.herbs_compounds, sig(x), "Intersection of herbs compounds")
    } else {
      p.herbs_targets <- NULL
      p.herbs_compounds <- NULL
    }
    ## set plots and tables
    plots <- namel(p.herbs_targets, p.herbs_compounds)
    tables <- namel(compounds_targets_annotation)
    if (!is.null(disease)) {
      x$disease <- disease
      disease_targets <- get_from_genecards(disease, score = disease.score)
      .add_internal_job(.job_genecard())
      x@params$disease_targets <- disease_targets
      disease_targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol",
        disease_targets$Symbol)
      .add_internal_job(.job_biomart())
      tables <- c(tables, namel(disease_targets_annotation))
      lst <- list(
        targets = targets,
        targets_annotation = targets_annotation$hgnc_symbol,
        disease_targets = disease_targets$Symbol,
        disease_targets_annotation = disease_targets_annotation$hgnc_symbol
      )
      p.targets_disease <- new_upset(lst = lst)
      plots <- c(plots, namel(p.targets_disease))
      x@params$ppi_used <- dplyr::filter(
        targets_annotation,
        hgnc_symbol %in% dplyr::all_of(disease_targets_annotation$hgnc_symbol)
      )
    }
    data.allu <- dplyr::select(easyRead, Herb_pinyin_name, Ingredient.name, Target.name)
    data.allu <- dplyr::mutate(data.allu, Target.name = ifelse(is.na(Target.name), "Unkown", Target.name))
    data.allu <- dplyr::distinct(data.allu)
    x$data.allu <- data.allu
    p.pharm <- plot_network.pharm(data.allu, seed = x$seed, HLs = HLs)
    p.pharm <- .set_lab(p.pharm, sig(x), "network pharmacology visualization")
    plots <- c(plots, namel(p.pharm))
    x@plots[[ 3 ]] <- plots
    x@tables[[ 3 ]] <- tables
    return(x)
  })

# setMethod("merge", signature = c(x = "job_herb", y = "list"),
#   function(x, y){
#   })


rstyle <- function(get = "pal", seed = NULL, n = 1L) {
  if (get == "pal") {
    set <- lapply(ggsci:::ggsci_db[1:19], function(x) unname(x[[1]]))
  } else if (get == "shape.c") {
    set <- as.list(15:19)
  } else if (get == "shape.f") {
    set <- as.list(21:25)
  } else if (get == "theme") {
    set <- list(theme_grey(), theme_minimal(), theme_bw(), theme_light())
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  res <- sample(set, n)
  if (n == 1) {
    return(res[[ 1 ]])
  }
  res
}

plot_network.pharm <- function(data, f.f = 2.5, f.f.mul = .7, f.f.sin = .2, f.ax4 = 2, seed = sample(1:10, 1), HLs = NULL,
  ax1 = "Herb", ax2 = "Compound", ax3 = "Target", ax4 = NULL, less.label = T,
  ax2.level = NULL, ax4.level = NULL, decImport.ax4Level = F, lab.fill = "",
  edge_width = .1, force.ax1 = NULL, ax3.layout = c("spiral", "grid"), ax4.layout = c("circle", "linear"),
  HLs.label = T)
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
    if (ax3.layout == "spiral" && nrow(nodesTgt) < 30) {
      message("Nodes of 'ax3' too small, set layout as: grid")
      ax3.layout <- "grid"
    }
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
      edge_color = "lightblue")
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

get_herb_data <- function(herb = .prefix("HERB_herb_info.txt", "db"),
  component = .prefix("HERB_ingredient_info.txt", "db"),
  target = .prefix("HERB_target_info.txt", "db"))
{
  db <- lapply(namel(herb, component, target),
    function(file) {
      tibble::as_tibble(suppressWarnings(data.table::fread(file)))
    })
  names <- colnames(db$component)
  colnames(db$component) <- c(names[-1], "extra_id")
  db
}

download_herbTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Gene Targets", ...)
{
  download_herbCompounds(link, ids, heading = heading, ...)
}

download_compoundTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Ingredient"),
  heading = "Related Gene Targets", ...)
{
  download_herbCompounds(link, ids, urls, heading = heading, ...)
}

download_herbCompounds <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Ingredients", tempdir = "download")
{
  get_all <- function() list.files(tempdir, "\\.xlsx$", full.names = T)
  checks <- get_all()
  if (length(checks) > 0)
    file.remove(checks)
  link$open()
  jobn <- 0
  res <- lapply(urls,
    function(url) {
      link$navigate(url)
      jobn <<- jobn + 1
      if (jobn > 1) {
        link$refresh()
      }
      Sys.sleep(.5)
      ## check whether all the web page is all loaded
      isOK <- F
      if (grpl(url, "Herb$")) {
        checkWhat <- "Herb id"
      } else if (grpl(url, "Ingredient$")) {
        checkWhat <- "Ingredient id"
      }
      while (!isOK) {
        status <- try(link$findElement(
            using = "xpath",
            value = paste0("//main//div//ul//span//b[text()='", checkWhat, "']")
            ), T)
        if (!inherits(status, "try-error")) {
          isOK <- T
        } else {
          message("Guess the page is loading...")
          Sys.sleep(1)
        }
      }
      res <- try(silent = T, {
        ele <- link$findElement(
          using = "xpath",
          value = paste0("//main//div//h4[text()='", heading, "']/../div//button")
        )
        Sys.sleep(.5)
        ele$sendKeysToElement(list("Download", key = "enter"))
        })
      Sys.sleep(.5)
      if (inherits(res, "try-error")) {
        writeLines("", paste0(tempdir, "/empty_", jobn, ".xlsx"))
      }
      checks <- get_all()
      cli::cli_alert_info(paste0("Job number ", jobn))
      ## sometimes these is no response
      while (length(checks) < jobn) {
        cli::cli_alert_info("Retry remotes response ...")
        ele$sendKeysToElement(list("Download", key = "enter"))
        Sys.sleep(1)
        checks <- get_all()
      }
    })
  link$close()
}

start_drive <- function(command = paste0("java -jar ", .prefix("/selenium.jar", "op")),
  port = 4444, extra = NULL, browser = c("firefox", "chrome"), download.dir = "download", ...)
{
  system(paste(command, "-port", port, extra), wait = F)
  new_link(port = port, browser = match.arg(browser), download.dir = download.dir, ...)
}

end_drive <- function(pattern = "[0-9]\\s*java.*-jar") {
  if (.Platform$OS.type == "unix") {
    if (getOption("end_drive", F)) {
      tmp <- tempfile(fileext = ".txt")
      system(paste0("ps aux > ", tmp))
      text <- readLines(tmp)
      text <- text[ grepl(pattern, text) ]
      pid <- stringr::str_extract(text, "(?<=\\s)[0-9]{1,}(?=\\s)")
      system(paste("kill", paste0(pid, collapse = " ")))
    }
  }
}

get_table.html <- function(x, ...) {
  if (file.exists(x)) {
    str <- readLines(x)
  } else {
    str <- x
  }
  ht <- e(XML::htmlParse(str))
  XML::readHTMLTable(ht, ...)
}

collateFiles_herbs <- function(ids,
  file.pattern = "\\.xlsx$", index.pfun = file_seq.by_time,
  from = "download", to = "herbs_ingredient", suffix = ".xlsx", .id = "herb_id", readFrom = NULL)
{
  if (is.null(readFrom)) {
    if (!file.exists(to)) {
      args <- as.list(environment())
      args$readFrom <- NULL
      files <- do.call(collateFiles, args)
      isOrdered <- T
    } else {
      files <- list.files(to, file.pattern, full.names = T)
      isOrdered <- F
    }
  } else {
    files <- readFrom
    isOrdered <- F
  }
  data <- lapply(files,
    function(file) {
      res <- try(openxlsx::read.xlsx(file), T)
      if (inherits(res, "try-error")) {
        tibble::tibble()
      } else res
    })
  if (isOrdered) {
    names(data) <- ids
  } else {
    names(data) <- get_realname(files)
  }
  data <- frbind(data, idcol = .id, fill = T)
  if (nrow(data)) {
    data <- dplyr::relocate(data, !!rlang::sym(.id))
  }
  data
}

collateFiles <- function(ids,
  file.pattern, index.pfun = file_seq.by_time,
  from, to, suffix, ...)
{
  files <- list.files(from, file.pattern, full.names = T)
  index <- index.pfun(files)
  files <- files[order(index)][1:length(ids)]
  files <- rev(files)
  dir.create(to, F, T)
  files <- lapply(1:length(ids),
    function(n) {
      file.copy(files[ n ], file <- paste0(to, "/", ids[n], suffix), T)
      file
    })
  files
}

file_seq.by_time <- function(files) {
  info <- fs::file_info(files)
  info <- dplyr::mutate(info, .INDEX = 1:nrow(info))
  info <- dplyr::arrange(info, dplyr::desc(modification_time))
  info <- dplyr::mutate(info, .SEQ = 1:nrow(info))
  info <- dplyr::arrange(info, .INDEX)
  info$.SEQ
}

format_index.by_num <- function(index) {
  index <- gsub("\\(|\\)", "", index)
  index <- ifelse(is.na(index), "0", index)
  as.integer(index)
}

new_link <- function(port = 4444L, browser = c("firefox", "chrome"), addr = "localhost", download.dir = "download")
{
  if (!dir.exists(download.dir)) {
    dir.create(download.dir)
  }
  browser <- match.arg(browser)
  if (browser == "firefox") {
    download.dir <- normalizePath(download.dir)
    # about:config
    if (!is.null(path_prof <- getOption("FirefoxProfile"))) {
      prof <- RSelenium::getFirefoxProfile(path_prof)
    } else {
      prof <- RSelenium::makeFirefoxProfile(
        list('permissions.default.image' = 2L,
          'browser.download.folderList' = 2L,
          'browser.download.manager.showWhenStarting' = F,
          'browser.download.lastDir' = download.dir,
          'browser.download.dir' = download.dir
        )
      )
    }
    link <- RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = prof
    )
    link
  } else if (browser == "chrome") {
    RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = list()
    )
  }
}

merge.componentsGenes <- function(data, genes.lst) {
  names(genes.lst) <- data[[1]]
  genes.lst <- lapply(genes.lst,
    function(lst) {
      if (nrow(lst$ids) == 0) {
        return(NULL)
      } else {
        lst$ids$genes <- lst$data
        return(lst$ids)
      }
    })
}

get_tcm.components <- function(id, key = "Components",
  link_prefix = "http://www.tcmip.cn/ETCM/index.php/Home/Index/yc_details.html?id=",
  src = NULL)
{
  if (is.null(src)) {
    src <- get_tcm.base(id, link_prefix)
  }
  pos <- grep(paste0("^\\s*\\{\"ID*.*", key), src)
  data <- src[pos]
  ids.data <- extract_id(data)
  data <- gsub("<[^<^>]*>", "##", data)
  data <- gsub("^.*?####|####[^#]*$", "", data)
  data <- strsplit(data, "##, ##")[[1]]
  list(data = data, ids = ids.data, src = src)
}

extract_id <- function(ch) {
  links <- stringr::str_extract_all(ch, "(?<=href=\').*?(?=\')")[[1]]
  ids <- stringr::str_extract(links, "[^=]*$")
  data.frame(links = links, ids = ids)
}

get_tcm.componentToGenes <- function(ids, key = "Candidate Target Genes", dep = T) 
{
  link_prefix <- "http://www.tcmip.cn/ETCM/index.php/Home/Index/cf_details.html?id="
  lst <- pbapply::pblapply(ids,
    function(id) {
      data <- get_tcm.components(id, key, link_prefix)
      if (dep) {
        data$src <- NULL
      }
      data
    })
  lst
}

get_tcm.base <- function(id, link_prefix) {
  link <- paste0(link_prefix, id)
  src <- xml2::read_html(link)
  data <- rvest::html_text(src)
  data <- strsplit(data, "\r\n")[[1]]
  data
}

get_coords.sin <- function(num = 100, nWave = 4) {
  x <- seq(0, nWave * 2 * pi, length.out = num)
  tibble::tibble(x = x, y = sin(x))
}

get_coords.spiral <- function(num = 100, nCir = 3, minRad = 1, maxRad = 3) {
  ang <- seq(0, by = 2 * pi * nCir / num, length.out = num)
  rad <- seq(minRad, maxRad, length.out = num)
  tibble::tibble(x = sin(ang) * rad, y = cos(ang) * rad)
}

