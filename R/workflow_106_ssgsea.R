# ==========================================================================
# workflow of ssgsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ssgsea <- setClass("job_ssgsea", 
  contains = c("job"),
  prototype = prototype(
    pg = "ssgsea",
    info = c("https://www.bioconductor.org/packages/release/bioc/html/GSVA.html"),
    cite = "[@Systematic_RNA_Barbie_2009]",
    method = "",
    tag = "ssgsea",
    analysis = "ssGSEA 单样本GSEA富集分析"
    ))

setGeneric("asjob_ssgsea", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_ssgsea"))

setMethod("asjob_ssgsea", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = .guess_symbol(x), 
    use.format = TRUE, ...)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    cli::cli_alert_info("extract_unique_genes.job_limma")
    object <- extract_unique_genes.job_limma(
      x, use.filter, use, use.format = use.format, ...
    )
    mtx <- object$E
    genes <- object$genes
    if (use.format) {
      rownames(mtx) <- gname(genes[[use]])
    } else {
      rownames(mtx) <- genes[[use]]
    }
    metadata <- object$targets
    if (x@step >= 2L) {
      contrasts <- .get_versus_cell(
        names(x@tables$step2$tops), NULL, unlist = FALSE
      )
    } else {
      contrasts <- NULL
      message("Can not match 'contrasts' in `x`.")
    }
    # mtx <- mtx[!is.na(rownames(mtx)), ]
    x <- job_ssgsea(mtx)
    x$metadata <- metadata
    x$contrasts <- contrasts
    return(x)
  })

job_ssgsea <- function(x)
{
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    BiocManager::install("GSVA")
  }
  if (!is(x, "GsvaExprData")) {
    stop('!is(x, "GsvaExprData").')
  }
  .job_ssgsea(object = x)
}

setMethod("step0", signature = c(x = "job_ssgsea"),
  function(x){
    step_message("Prepare your data with function `job_ssgsea`.")
  })

setMethod("step1", signature = c(x = "job_ssgsea"),
  function(x, mode = c("matrisome"), org = c("human", "mouse"), sets)
  {
    step_message("Calculate ssGSEA enrichment score.")
    if (missing(sets)) {
      mode <- match.arg(mode)
      if (mode == "matrisome") {
        db <- .job_matrisome(sig = x@sig)
        db <- step1(db, org)
        x <- snapAdd(x, db)
        # extract results.
        x$db <- db$db
        sets <- feature(db)
        x <- snapAdd(x, "将{snap(sets)}用于 ssGSEA 富集分析。")
        sets <- list(matrisome = unlist(sets, use.names = FALSE))
      }
    }
    if (is(sets, "list")) {
      sets <- mapply(names(sets), sets, SIMPLIFY = FALSE,
        FUN = function(name, genes) {
          GSEABase::GeneSet(genes, setIdentifier = name) 
        })
      sets <- e(GSEABase::GeneSetCollection(sets))
    }
    param <- e(GSVA::ssgseaParam(object(x), sets))
    res <- e(GSVA::gsva(param))
    data <- dplyr::select(x$metadata, group, sample)
    if (!identical(data$sample, colnames(res))) {
      stop('!identical(data$sample, colnames(res)).')
    }
    data <- cbind(
      data, setNames(data.frame(t(res), check.names = FALSE), "Stromal_score")
    )
    x$data <- tidyr::pivot_longer(data, Stromal_score, names_to = "type", values_to = "score")
    if (!is.null(x$contrasts)) {
      p.scores <- lapply(x$contrasts, 
        function(group) {
          data <- dplyr::filter(x$data, group %in% !!group)
          data <- dplyr::mutate(
            data, group = factor(group, levels = !!group)
          )
          p <- .map_boxplot2(
            x$data, TRUE, y = "score", ylab = "Stromal score", ids = "type"
          )
          wrap(p, 2.5, 3)
        })
      names(p.scores) <- vapply(
        x$contrasts, bind, character(1), co = " "
      )
      p.scores <- .set_lab(
        p.scores, sig(x), names(p.scores), "boxplot of stromal score"
      )
      p.scores <- setLegend(p.scores, glue::glue("为各分组 {names(p.scores)} 基质评分箱形图。"))
      x <- plotsAdd(x, p.scores = p.scores)
    }
    object(x) <- NULL
    x <- methodAdd(x, "以 R 包 `GSVA` ({packageVersion('GSVA')}) 用于 ssGSEA 分析。")
    return(x)
  })
