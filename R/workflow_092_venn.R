# ==========================================================================
# workflow of venn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_venn <- setClass("job_venn", 
  contains = c("job"),
  prototype = prototype(
    pg = "venn",
    info = c("..."),
    cite = "",
    method = "",
    tag = "venn",
    analysis = "Venn 交集"
    ))

.job_vennDEGs <- setClass("job_vennDEGs", contains = c("job_venn"))

setGeneric("asjob_venn",
  function(x, ...) standardGeneric("asjob_venn"))

job_vennDEGs <- function(pattern, exclude = NULL,
  mode = c("raw", "pro"), name = "guess", ...,
  pattern_dataset = ".*(GSE[0-9]+).*", gp = NULL)
{
  mode <- match.arg(mode)
  if (mode == "raw") {
    fun_extract <- function(x) x@tables$step2$tops[[1]][, 1:3]
    snapAdd_onExit("x", "数据集的差异表达基因的交集。")
  } else {
    stop(glue::glue("..."))
  }
  degs <- collate_dataset_DEGs(
    pattern, name = name, exclude = exclude, fun_extract = fun_extract, ...
  )
  degs_versus <- collate(
    pattern, function(x) names(x@tables$step2$tops)[1], exclude, ...
  )
  projects <- collate(
    pattern, function(x) x$project, exclude, ...
  )
  names(degs_versus) <- s(names(degs_versus), pattern_dataset, "\\1")
  degs <- split(degs$symbol, degs$Dataset)
  names(degs) <- s(names(degs), pattern_dataset, "\\1")
  x <- .job_vennDEGs(job_venn(lst = degs))
  x$degs_versus <- degs_versus
  if (!is.null(gp)) {
    metadata <- as_df.lst(x$degs_versus, "project", "versus")
    groups <- group_strings(unlist(x$degs_versus), gp, "versus")
    x$metadata <- map(
      metadata, "versus", groups, "versus", "group", col = "group"
    )
    x$metadata <- tibble::as_tibble(x$metadata)
    x$metadata <- set_lab_legend(
      x$metadata, "metadata of mutiple datasets", "数据集的元数据信息 (分组信息) "
    )
    x$groups <- nl(x$metadata$project, x$metadata$group)
  }
  x$projects <- projects
  return(x)
}

job_venn <- function(..., mode = c("key", "candidates"),
  analysis = NULL, lst = NULL, fun_map = function(x) x)
{
  if (is.null(lst)) {
    object <- list(...)
  } else {
    object <- lst
  }
  nature <- "基因集"
  if (all(vapply(object, is, logical(1), "feature"))) {
    snaps <- paste0(
      "- ", vapply(
        object, snap, character(1), enumerate = FALSE, unlist = TRUE
      )
    )
    methodAdd_onExit("x", "数据集为：\n\n{bind(snaps, co = '\n')}\n\n\n\n")
    nature <- object[[1]]@nature
    message(glue::glue("Use nature as: {nature}"))
    object <- lapply(object, function(x) fun_map(unlist(x@.Data)))
    if (is.null(analysis)) {
      mode <- match.arg(mode)
      mode <- switch(mode, key = "关键", candidates = "候选")
      analysis <- glue::glue("{mode}{nature}")
    }
  } else {
    message(glue::glue("Some were not 'feature', be carefull! If you need automatic snap ..."))
  }
  x <- .job_venn(object = lapply(object, unlist))
  x <- methodAdd(x, "以 R 包 `ggVennDiagram` ⟦pkgInfo('ggVennDiagram')⟧ 对{nature}取交集。")
  x$nature <- nature
  x$analysis <- analysis
  return(x)
}

setMethod("step0", signature = c(x = "job_venn"),
  function(x){
    step_message("Prepare your data with function `job_venn`.")
  })

setMethod("step1", signature = c(x = "job_venn"),
  function(x, ...){
    step_message("Intersection.")
    p.venn <- new_venn(lst = object(x), force_upset = FALSE, ...)
    p.venn <- set_lab_legend(
      p.venn,
      glue::glue("{x@sig} intersection of {bind(names(object(x)), co = ' with ')}"),
      glue::glue("{bind(names(object(x)))} 交集维恩图|||不同颜色圆圈代表不同数据集，中间重叠部分表示同时存在多个集合中。图中 {length(p.venn$ins)} 交集为：{less(p.venn$ins, 20)}。")
    )
    x$.append_heading <- FALSE
    if (identical(parent.frame(1), .GlobalEnv)) {
      job_append_heading(
        x, heading = glue::glue(
          "汇总: ", bind(names(object(x)), co = " + ")
        )
      )
    }
    if (length(p.venn$ins) < 10) {
      iter <- glue::glue(" ({bind(p.venn$ins)}) ")
    } else {
      iter <- ""
    }
    x <- snapAdd(x, "对{bind(names(object(x)))} 取交集，得到{length(p.venn$ins)}个交集{iter}{aref(p.venn)}。")
    x <- plotsAdd(x, p.venn)
    if (!is.null(x$analysis)) {
      feature(x) <- as_feature(p.venn$ins, x$analysis, nature = x$nature, ...)
    } else {
      x$.feature <- as_feature(p.venn$ins, x, nature = x$nature, ...)
    }
    return(x)
  })

new_venn <- function(..., lst = NULL, wrap = TRUE,
  fun_pre = rm.no, force_upset = NULL, n = NULL)
{
  if (!is.null(lst) && length(list(...))) {
    lst <- c(lst, list(...))
  }
  if (is.null(lst)) {
    lst <- list(...)
  }
  lst <- lst_clear0(lst)
  lst <- lapply(lst, function(x) as.character(fun_pre(x)))
  if (is.null(force_upset)) {
    if (length(lst) > 3) {
      force_upset <- TRUE
    } else {
      force_upset <- FALSE
    }
  }
  if (force_upset) {
    p <- ggVennDiagram::ggVennDiagram(
      lst, force_upset = TRUE, nintersects = n
    )
    p$plotlist[[2]]$layers[[1]]$geom_params$width <- .4
    p$plotlist[[3]]$layers[[1]]$geom_params$width <- .7
  } else {
    p <- ggVennDiagram::ggVennDiagram(lst, label_percent_digit = 1, edge_size = .3) +
      scale_fill_gradient(low = "grey95", high = sample(color_set(), 1)) +
      theme_void() +
      guides(fill = "none") +
      theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
      theme()
    for (i in rev(seq_along(p$layers))) {
      if (is(p$layers[[i]]$geom, "GeomLabel")) {
        which <- i
      } else if (is(p$layers[[i]]$geom, "GeomText")) {
        whichText <- i
        break
      }
    }
    shouldChange <- !grpl(p$layers[[which]]$data$id, "/")
    p$layers[[which]]$data[shouldChange, ] <- dplyr::mutate(
      p$layers[[which]]$data[shouldChange, ], both = paste0(name, "\n", count, " (", percent, ")")
    )
    p$layers[[which]]$data[!shouldChange, ] <- dplyr::mutate(
      p$layers[[which]]$data[!shouldChange, ], both = paste0(count, " (", percent, ")")
    )
    p$layers[[whichText]] <- NULL
  }
  if (wrap) {
    p <- wrap(p, if (force_upset) 5 else 3, 4)
  }
  attr(p, "ins") <- ins <- ins(lst = lst)
  attr(p, "lich") <- new_lich(list(All_intersection = ins))
  lab(p) <- paste0("Intersection of ", paste0(names(lst), collapse = " with "))
  p <- setLegend(p, "为 {bind(names(lst))} 各自交集。")
  p
}

setMethod("meta", signature = c(x = "job_vennDEGs"),
  function(x, group = NULL, bind = TRUE, arrange = TRUE, get = "project")
  {
    if (is.null(x$metadata)) {
      stop('is.null(x$metadata).')
    }
    if (arrange) {
      x$metadata <- dplyr::arrange(x$metadata, group)
    }
    if (!is.null(group)) {
      res <- dplyr::filter(x$metadata, group %in% !!group)
      if (!nrow(res)) {
        stop('!nrow(res). No any results.')
      }
      res <- res[[ get ]]
    } else {
      res <- x$metadata[[ get ]]
    }
    if (bind) {
      bind(res)
    } else {
      res
    }
  })

setMethod("feature", signature = c(x = "job_vennDEGs"),
  function(x, group, intersect = TRUE){
    if (missing(group)) {
      callNextMethod(x)
    } else {
      projects <- meta(x, group, bind = FALSE)
      sets <- object(x)[ names(object(x)) %in% projects ]
      if (intersect) {
        sets <- ins(lst = sets)
      }
      as_feature(
        sets, x, analysis = " Venn 交集 ({bind(group)}: {bind(projects)})"
      )
    }
  })

alias_intersect_multi <- function(..., mode = c("main", "all", "index"), main = 1L, sep = "///")
{
  require(data.table)

  mode <- match.arg(mode)

  inputs <- list(...)
  k <- length(inputs)

  # Expand all inputs into long alias tables
  dt_list <- lapply(seq_along(inputs), function(i) {
    dt <- .expand_alias_dt(inputs[[i]], paste0("s", i, "_"), sep)
    dt[, set_id := i]
    dt
  })

  # Combine all tables
  dt_all <- rbindlist(dt_list)

  # Count how many distinct sets each alias appears in
  alias_sets <- dt_all[, .(set_count = uniqueN(set_id)), by = alias]

  # Keep only aliases present in all sets
  valid_alias <- alias_sets[set_count == k, alias]

  # If no overlap, return empty structure
  if (length(valid_alias) == 0) {
    if (mode == "main") {
      return(inputs[[main]][ integer(0) ])
    } else {
      return(vector("list", k))
    }
  }

  # Filter hits
  hit <- dt_all[alias %in% valid_alias]

  # Extract unique indices per set
  res_idx <- hit[, .(idx = unique(idx)), by = set_id]

  # Return based on mode
  if (mode == "main") {
    idx <- res_idx[set_id == main, idx]
    inputs[[main]][idx]
  } else if (mode == "index") {
    split(res_idx$idx, res_idx$set_id)
  } else if (mode == "all") {
    lapply(seq_along(inputs), function(i) {
      idx <- res_idx[set_id == i, idx]
      inputs[[i]][idx]
    })
  }
}


# Main function: compute intersection based on alias overlap
alias_intersect <- function(x, y, mode = c("x", "y", "both", "index"), sep = "///")
{
  require(data.table)

  mode <- match.arg(mode)

  # Expand both inputs into long alias tables
  x_dt <- .expand_alias_dt(x, "x_", sep)
  y_dt <- .expand_alias_dt(y, "y_", sep)

  # Set keys for fast join on alias
  setkey(x_dt, alias)
  setkey(y_dt, alias)

  # Perform join to find overlapping aliases
  hit <- x_dt[y_dt, nomatch = 0]

  # Extract matched indices directly (faster than gid parsing)
  x_idx <- unique(hit$idx)
  y_idx <- unique(hit$i.idx)

  # Return results based on mode
  if (mode == "x") {
    x[x_idx]
  } else if (mode == "y") {
    y[y_idx]
  } else if (mode == "both") {
    list(x = x[x_idx], y = y[y_idx])
  } else {
    list(x_index = x_idx, y_index = y_idx)
  }
}


# Expand input into a long-format data.table (robust to vector/list, names optional)
.expand_alias_dt <- function(lst, prefix, sep = "\\s*///\\s*") {

  n <- length(lst)

  # Handle names safely (vector without names → auto-generate)
  nm <- names(lst)
  if (is.null(nm)) {
    nm <- paste0("V", seq_len(n))
  } else {
    empty <- nm == "" | is.na(nm)
    nm[empty] <- paste0("V", which(empty))
  }

  # Build base table
  dt <- data.table(
    name = nm,
    idx  = seq_len(n)
  )

  # Split aliases (vectorized, avoid lapply overhead)
  alias_list <- strsplit(as.character(lst), sep)

  # Flatten once
  alias_vec <- trimws(unlist(alias_list, use.names = FALSE))

  # Map aliases back to original indices
  lens <- lengths(alias_list)
  dt <- dt[rep.int(seq_len(n), lens)]

  dt[, alias := alias_vec]

  # Remove invalid aliases
  dt <- dt[alias != "" & !is.na(alias)]

  # Deduplicate within each group
  dt <- unique(dt, by = c("idx", "alias"))

  # Create group ID (kept for compatibility, though no longer used downstream)
  dt[, gid := paste0(prefix, idx)]

  dt[, .(gid, name, idx, alias)]
}

