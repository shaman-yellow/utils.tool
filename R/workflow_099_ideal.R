# ==========================================================================
# workflow of ideal
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ideal <- setClass("job_ideal", 
  contains = c("job"),
  prototype = prototype(
    pg = "ideal",
    info = c(""),
    cite = "",
    method = "",
    tag = "ideal",
    analysis = "Collate 汇总分析结果"
    ))

job_ideal <- function(mode = c("high", "low"))
{
  mode <- match.arg(mode)
  x <- .job_ideal()
  x$mode <- mode
  return(x)
}

setMethod("step0", signature = c(x = "job_ideal"),
  function(x){
    step_message("Prepare your data with function `job_ideal`.")
  })

is.collate <- function(x) {
  is(x, "df") && !is.null(attr(x, "__COLLATE_NAME__"))
}

validCollates <- function(args) {
  if (!all(vapply(args, is.collate, logical(1)))) {
    stop('!all(vapply(args, is.collate, logical(1))).')
  }
}

setMethod("step1", signature = c(x = "job_ideal"),
  function(x, ..., pattern_dataset = ".*(GSE[0-9]+).")
  {
    step_message("Collate DEGs results.")
    args <- list(...)
    validCollates(args)
    names(args) <- vapply(args, function(x) attr(x, "__COLLATE_NAME__"), character(1))
    t.collated_DEGs <- rbind_list(args)
    t.collated_DEGs <- setLegend(t.collated_DEGs, "为差异分析 DEGs 以及数据集来源汇总。")
    refs <- lapply(
      args, function(x) s(unique(x$Dataset), pattern_dataset, "\\1")
    )
    all_degs <- lapply(args,
      function(object) {
        if (x$mode == "high") {
          dplyr::filter(object, logFC > 0)
        } else if (x$mode == "low") {
          dplyr::filter(object, logFC < 0)
        }
      })
    all_degs_symbols <- lapply(all_degs, function(x) gname(x$symbol))
    x$res_degs <- namel(all_degs, all_degs_symbols, refs)
    x <- tablesAdd(x, t.collated_DEGs)
    return(x)
  })

setMethod("step2", signature = c(x = "job_ideal"),
  function(x, ..., pattern_dataset = ".*(GSE[0-9]+).")
  {
    step_message("Collate survival results.")
    args <- list(...)
    validCollates(args)
    names(args) <- vapply(args, function(x) attr(x, "__COLLATE_NAME__"), character(1))
    refs <- lapply(
      args, function(x) s(unique(x$Dataset), pattern_dataset, "\\1")
    )
    all_genes <- lapply(args,
      function(object) {
        if (x$mode == "high") {
          dplyr::filter(object, group_low_survival == "High")
        } else if (x$mode == "low") {
          dplyr::filter(object, group_low_survival == "Low")
        }
      })
    all_genes_symbols <- lapply(all_genes, function(x) gname(x$name))
    x$res_surv <- namel(all_genes, all_genes_symbols, refs)
    return(x)
  })

setMethod("step3", signature = c(x = "job_ideal"),
  function(x){
    step_message("Gather results.")
    if (TRUE) {
      sign <- switch(x$mode, high = "Up_", low = "Down_")
      degs <- x$res_degs$all_degs_symbols
      names(degs) <- paste0(sign, names(degs))
    }
    if (TRUE) {
      sign <- switch(x$mode, high = "Low_", low = "High_")
      survs <- x$res_surv$all_genes_symbols
      names(survs) <- paste0(sign, names(survs))
    }
    lst <- c(degs, survs)
    p.venn <- new_venn(lst = lst, force_upset = length(lst) >= 3)
    refs <- lapply(c(x$res_degs$refs, x$res_surv$refs), bind)
    names(refs) <- paste0("Refer_dataset_", names(refs))
    p.venn$lich <- .lich(c(p.venn$lich, refs))
    x <- plotsAdd(x, p.intersectionOfConditionalGenes = p.venn)
    sign <- switch(x$mode, high = "高表达", low = "低表达")
    x <- snapAdd(x, "筛选{sign}差异基因 ({bind(names(degs))})，且生存分析为{sign}预后不良。")
    return(x)
  })

