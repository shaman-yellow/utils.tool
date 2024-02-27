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
    method = "Website `HERB` <http://herb.ac.cn/> used for data source"
    ))

job_herb <- function(herbs, db = get_herb_data())
{
  x <- .job_herb(object = db)
  x@params$herbs <- herbs
  herbs_info <- filter(object(x)$herb, Herb_cn_name %in% !!params(x)$herbs)
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
  function(x, tempdir = "download"){
    step_message("Dowload compounds of herbs.")
    x$tempdir <- tempdir
    link <- start_drive(browser = "firefox", download.dir = x$tempdir)
    Sys.sleep(3)
    tryCatch(
      download_herbCompounds(link, params(x)$herbs_info$Herb_, tempdir = x$tempdir),
      finally = end_drive()
    )
    if (file.exists("herbs_ingredient")) {
      unlink("herbs_ingredient", T, T)
    }
    herbs_compounds <- moveToDir_herbs(params(x)$herbs_info$Herb_, from = x$tempdir)
    x@tables[[ 1 ]] <- namel(herbs_compounds)
    return(x)
  })

setMethod("step2", signature = c(x = "job_herb"),
  function(x, group_number = 50, group_sleep = 3, db = .prefix(), removeExistsTemp = F)
  {
    step_message("Dowload targets of compounds")
    all_ids <- unique(x@tables$step1$herbs_compounds$Ingredient.id)
    merge_data <- F
    if (!is.null(db)) {
      dbfile <- list.files(db, "^HBIN.*.xlsx$", full.names = T, recursive = T)
      names(dbfile) <- get_realname(dbfile)
      dbfile <- dbfile[ !duplicated(dbfile) ]
      has_ids <- all_ids[ all_ids %in% names(dbfile) ]
      if (length(has_ids)) {
        files <- dbfile[ names(dbfile) %in% has_ids ]
        data <- moveToDir_herbs(readFrom = files, .id = "Ingredient_id", from = x$tempdir)
        all_ids <- all_ids[ !all_ids %in% has_ids ]
        merge_data <- T
      }
    }
    if (length(all_ids)) {
      all_ids <- grouping_vec2list(all_ids, group_number, T)
      link <- start_drive(browser = "firefox")
      Sys.sleep(3)
      lst <- pbapply::pblapply(all_ids,
        function(ids) {
          if (file.exists(.dir <- "compounds_target")) {
            .dir <- timeName("compounds_target")
          }
          to <- paste0(.dir, "/", attr(ids, "name"))
          if (!file.exists(to)) {
            download_compoundTargets(link, ids, tempdir = x$tempdir)
          } else {
            if (!removeExistsTemp) {
              isThat <- usethis::ui_yeah(paste0("Remove the exists tempdir:", to, " ?"))
            } else {
              isThat <- T
            }
            if (isThat) {
              unlink(to, recursive = T, force = T)
              download_compoundTargets(link, ids, tempdir = x$tempdir)
            }
          }
          data <- moveToDir_herbs(ids, to = to, .id = "Ingredient_id", from = x$tempdir)
          cli::cli_alert_info("Sleep between groups ...")
          Sys.sleep(group_sleep)
          data
        })
      lst <- data.table::rbindlist(lst, fill = T)
      end_drive()
    } else {
      lst <- list()
    }
    if (merge_data) {
      if (length(lst))
        lst <- list(lst, list(data))
      else
        lst <- list(data)
    }
    if (is(lst, "list")) {
      data <- data.table::rbindlist(lst, fill = T)
    } else {
      data <- lst
    }
    compounds_targets <- as_tibble(data)
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
      by.x = "Target.name", by.y = "hgnc_symbol", all.x = T
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
    data.allu <- dplyr::filter(data.allu, !is.na(Target.name))
    x$data.allu <- data.allu
    p.pharm <- plot_network.pharm(data.allu, seed = x$seed, HLs = HLs)
    p.pharm <- .set_lab(p.pharm, sig(x), "network pharmacology visualization")
    plots <- c(plots, namel(p.pharm))
    x@plots[[ 3 ]] <- plots
    x@tables[[ 3 ]] <- tables
    return(x)
  })

setMethod("intersect", signature = c(x = "job_herb", y = "job_herb"),
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

setMethod("map", signature = c(x = "job_herb", ref = "list"),
  function(x, ref, HLs = NULL, levels = NULL, lab.level = "Level", name = "dis", compounds = NULL, ...)
  {
    data <- x$data.allu
    if (!is.null(compounds)) {
      data <- dplyr::filter(data, Ingredient.name %in% !!compounds)
    }
    if (length(ref)) {
      data <- dplyr::filter(data, Target.name %in% !!unlist(ref, use.names = F))
      p.venn2dis <- new_venn(Diseases = unlist(ref, use.names = F), Targets = data$Target.name)
      x[[ paste0("p.venn2", name) ]] <- .set_lab(p.venn2dis, sig(x), "Targets intersect with targets of diseases")
    }
    p.pharm <- plot_network.pharm(data, HLs = HLs, ax2.level = levels, lab.fill = lab.level, ...)
    if (length(ref)) {
      x[[ paste0("p.pharm2", name) ]] <- .set_lab(p.pharm, sig(x), "network pharmacology with disease")
    } else {
      x$p.pharmMap <-  .set_lab(p.pharm, sig(x), "network pharmacology")
    }
    return(x)
  })

setMethod("slice", signature = c(x = "job_herb"),
  function(x, ...){
    x@params$herbs_info <- dplyr::slice(x@params$herbs_info, ...)
    return(x)
  })

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

plot_network.pharm <- function(data, f.f = 2.5, f.f.mul = .7, f.f.sin = .2, seed = 1, HLs = NULL,
  ax1 = "Herb", ax2 = "Compound", ax3 = "Target", less.label = T,
  ax2.level = NULL, lab.fill = "", edge_width = .1)
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
  prepare_data <- function(data) {
    colnames(data) <- c(ax1, ax2, ax3)
    nodes <- tidyr::gather(data, type, value)
    nodes <- dplyr::distinct(nodes)
    nodes <- dplyr::relocate(nodes, name = value, type)
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
    crds.Tgt <- get_layout(NULL, "grid", nodes = dplyr::filter(nodes, type == !!ax3))
    crds.Tgt <- shift(crds.Tgt, -max(crds.Tgt$x) / 2, -max(crds.Tgt$y) / 2)
    f.rsz <- max(crds.Tgt$x) * f.f
    if (!sherb) {
      crds.Hrb <- get_layout(NULL, "circle", nodes = dplyr::filter(nodes, type == !!ax1))
      crds.Hrb <- resize(crds.Hrb, f.rsz)
      lst <- split(ed.12sin, ed.12sin[[ ax1 ]])
      # coords for compounds from sigle herb
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
    crds <- lst_clear0(list(crds.Hrb, crds.ComSin, crds.ComMul, crds.Tgt))
    crds <- do.call(dplyr::bind_rows, crds)
    crds <- dplyr::distinct(crds, name, .keep_all = T)
    if (sherb) {
      use.cols <- 2
    } else {
      use.cols <- 1:2
    }
    edges <- lapply(use.cols,
      function(n) {
        edges <- dplyr::distinct(data, dplyr::pick(n:(n+1)))
        colnames(edges) <- c("source", "target")
        edges
      })
    edges <- do.call(dplyr::bind_rows, edges)
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
    fast_layout(edges, layout = layout, nodes = nodes)
  }
  x <- prepare_data(data)
  data <- as_tibble(data.frame(x))
  data <- dplyr::mutate(
    data, cent = ifelse(type == !!ax1, cent * 150,
      ifelse(grpl(type, !!ax2, ignore.case = T), cent * 20, cent)))
  data.tgt <- dplyr::filter(data, type == !!ax3)
  set.seed(seed)
  minSize <- .5
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
    geom_edge <- geom_edge_link(color = sample(color_set()[1:6], 1), alpha = .2, width = width)
  } else {
    geom_edge <- geom_edge_link(aes(edge_color = highlight, edge_width = highlight), alpha = .2)
  }
  if (less.label) {
    geom_label_1 <- ggrepel::geom_label_repel(
      data = dplyr::filter(data, cent > if (sherb) 100 else 1000),
      aes(x = x, y = y, label = name, color = type))
    geom_label_2 <- ggrepel::geom_label_repel(
      data = dplyr::filter(data.tgt, cent >= sort(cent, T)[10]),
      aes(x = x, y = y, label = name, color = type))
  } else {
    geom_label_1 <- ggrepel::geom_label_repel(
      data = data, aes(x = x, y = y, label = name, color = type))
    geom_label_2 <- geom_blank()
  }
  if (is.null(ax2.level)) {
    geom_node_1 <- geom_node_point(data = data, aes(x = x, y = y, color = type, size = cent, shape = type))
    geom_node_2 <- geom_blank()
    scale_fill_gradient <- geom_blank()
  } else {
    data.level <- dplyr::filter(data, .type == !!ax2)
    data.level <- map(data.level, "name", ax2.level, "name", "level", col = "level")
    geom_node_1 <- geom_node_point(data = data.level,
      aes(x = x, y = y, fill = level, size = cent), shape = 21,
      stroke = 0, color = "transparent"
    )
    data.nonlevel <- dplyr::filter(data, .type != !!ax2)
    geom_node_2 <- geom_node_point(data = data.nonlevel,
      aes(x = x, y = y, color = type, size = cent, shape = type))
    scale_fill_gradient <- scale_fill_gradientn(colors = rev(color_set2()))
  }
  p <- ggraph(x) + geom_edge +
    geom_node_1 + geom_node_2 +
    geom_label_1 + geom_label_2 +
    scale_size(range = c(minSize, 15)) +
    guides(size = "none", shape = "none") +
    scale_color_manual(values = rstyle("pal", seed)) +
    scale_edge_color_manual(values = c("Highlight" = "red", "Non-highlight" = "lightblue")) +
    scale_edge_width_manual(values = c("Highlight" = 1, "Non-highlight" = .1)) +
    scale_fill_gradient +
    theme_minimal() +
    labs(color = "Type", edge_color = "Highlight", edge_width = "Highlight", fill = lab.fill) +
    theme(axis.text = element_blank(),
      axis.title = element_blank()) +
    geom_blank()
  if (!is.null(HLs)) {
    data <- dplyr::filter(data, name %in% HLs)
    p <- p + geom_point(data = data, aes(x = x, y = y), shape = 21, color = "red", size = 10) +
      ggrepel::geom_label_repel(data = data, aes(x = x, y = y, label = name),
        size = 7, color = "black")
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

moveToDir_herbs <- function(ids,
  file.pattern = "\\.xlsx$", index.pfun = file_seq.by_time,
  from = "download", to = "herbs_ingredient", suffix = ".xlsx", .id = "herb_id", readFrom = NULL)
{
  if (is.null(readFrom)) {
    if (!file.exists(to)) {
      args <- as.list(environment())
      args$readFrom <- NULL
      files <- do.call(moveToDir, args)
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

moveToDir <- function(ids,
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
    prof <- RSelenium::makeFirefoxProfile(
      list('permissions.default.image' = 2L,
        'browser.download.folderList' = 2L,
        'browser.download.manager.showWhenStarting' = F,
        'browser.download.lastDir' = download.dir,
        'browser.download.dir' = download.dir
      )
    )
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

get_coords.spiral <- function(num = 100, nCir = 3, minRad = 1, maxRad = 3) {
  ang <- seq(0, by = 2 * pi * nCir / num, length.out = num)
  rad <- seq(minRad, maxRad, length.out = num)
  tibble::tibble(x = sin(ang) * rad, y = cos(ang) * rad)
}

