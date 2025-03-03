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

job_vennDEGs <- function(pattern, exclude = NULL,
  mode = c("raw", "pro"), name = "guess", ...,
  pattern_dataset = ".*(GSE[0-9]+).")
{
  mode <- match.arg(mode)
  if (mode == "raw") {
    fun_extract <- function(x) x@tables$step2$tops[[1]][, 1:3]
    snapAdd_onExit("x", "取不同数据集的差异表达基因的交集。")
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
  x$projects <- projects
  return(x)
}

job_venn <- function(..., lst = NULL)
{
  if (is.null(lst)) {
    object <- list(...)
  } else {
    object <- lst
  }
  if (all(vapply(object, is, logical(1), "feature"))) {
    snaps <- paste0("- ", vapply(object, snap, character(1)))
    snapAdd_onExit("x", "以下取交集：\n{bind(snaps, co = '\n')}")
    object <- lapply(object, function(x) unlist(x@.Data))
  }
  x <- .job_venn(object = object)
  return(x)
}

setMethod("step0", signature = c(x = "job_venn"),
  function(x){
    step_message("Prepare your data with function `job_venn`.")
  })

setMethod("step1", signature = c(x = "job_venn"),
  function(x, ...){
    step_message("Intersection.")
    p.venn <- new_venn(lst = object(x), ...)
    lab(p.venn) <- paste(sig(x), lab(p.venn))
    x@plots[[ 1L ]] <- namel(p.venn)
    x$.append_heading <- FALSE
    if (identical(parent.frame(1), .GlobalEnv)) {
      job_append_heading(
        x, heading = glue::glue(
          "交集: ", bind(names(object(x)), co = " + ")
        )
      )
    }
    x$.feature <- as_feature(p.venn$ins, x)
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
    p <- ggVennDiagram::ggVennDiagram(lst, label_percent_digit = 1) +
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
    p <- wrap(p, if (force_upset) 5 else 2, 3)
  }
  attr(p, "ins") <- ins <- ins(lst = lst)
  attr(p, "lich") <- new_lich(list(All_intersection = ins))
  lab(p) <- paste0("Intersection of ", paste0(names(lst), collapse = " with "))
  p <- setLegend(p, "将{bind(names(lst))} 取交集。")
  p
}


