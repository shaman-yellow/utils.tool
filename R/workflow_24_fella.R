# ==========================================================================
# workflow of fella
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_fella <- setClass("job_fella", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@FellaAnRPacPicart2018]",
    method = "R package `FELLA` used for metabolite enrichment analysis",
    tag = "enrich:fella",
    analysis = "FELLA 代谢物富集分析"
    ))

job_fella <- function(kegg) {
  x <- .job_fella()
  x$ids.lst <- list(ids = kegg)
  x
}

setGeneric("asjob_fella", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_fella"))

setMethod("asjob_fella", signature = c(x = "job_metabo"),
  function(x){
    mapped <- x@tables$step1$mapped
    x <- .job_fella()
    x$mapped <- mapped
    x$ids.lst <- list(ids = x$mapped$KEGG)
    x
  })

setMethod("step0", signature = c(x = "job_fella"),
  function(x){
    step_message("Prepare your data with function `job_fella`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_fella"),
  function(x, org = "hsa", db_fella = .prefix(name = "db")){
    step_message("Use KEGG ID to perform FELLA.
      "
    )
    db_dir <- init_fella(db_fella, org)
    if (is.null(x$db)) {
      db <- load_fella(db_dir)
    } else {
      db <- x$db
    }
    enrich.lst <- enrich_fella(x$ids.lst, db)
    names(enrich.lst) <- names(x$ids.lst)
    x$graph.lst <- graph_fella(enrich.lst, db)
    p.enrich <- lapply(x$graph.lst, function(x) wrap(plotGraph_fella(x), 10, 7))
    p.enrich <- .set_lab(p.enrich, sig(x), names(p.enrich), "enrichment with algorithm PageRank")
    t.enrich <- lapply(x$graph.lst, tibble::as_tibble)
    t.enrich <- .set_lab(t.enrich, sig(x), names(t.enrich), "data of enrichment with algorithm PageRank")
    # results tables
    x$inputs <- lapply(enrich.lst, FELLA::getInput)
    x$resTable <- lapply(enrich.lst,
      function(x) {
        res <- sapply(c("pagerank", "diffusion", "hypergeom"), simplify = F,
          function(method) {
            res <- FELLA::generateResultsTable(
              object = x, method = method, threshold = .1, data = db
            )
            if (method == "hypergeom") {
              input <- FELLA::getInput(x)
              res <- dplyr::mutate(res, Compound_Ratio = CompoundHits / length(input))
            }
            res
          })
        res
      })
    t.hypergeom <- lapply(x$resTable,
      function(x) {
        x <- as_tibble(x$hypergeom)
        dplyr::rename(x, Description = KEGG.name, Count = CompoundHits, pvalue = p.value)
      })
    t.hypergeom <- .set_lab(t.hypergeom, sig(x), "data of enrichment with algorithm Hypergeom")
    sig <- sig(x)
    p.hypergeom <- lapply(t.hypergeom,
      function(x) {
        x <- plot_kegg(x, use = "pvalue", ratio = "Compound_Ratio")
        lab(x) <- paste0(sig, " Compounds hypergeom ", lab(x))
        x
      })
    x@tables[[ 1 ]] <- namel(t.enrich, t.hypergeom)
    x@plots[[ 1 ]] <- namel(p.enrich, p.hypergeom)
    x$org <- org
    x$db_dir <- db_dir
    x$enrich.lst <- enrich.lst
    return(x)
  })

setMethod("clear", signature = c(x = "job_fella"),
  function(x, ...){
    x@params$db <- NULL
    callNextMethod(x, ..., name = "fl")
  })

setMethod("map", signature = c(x = "job_enrich", ref = "job_fella"),
  function(x, ref, use.x = "kegg", use.ref = "hypergeom") {
    if (use.ref == "hypergeom") {
      ref <- ref@tables$step1$t.hypergeom[[1]][[ "KEGG.id" ]]
    } else {
      stop("...")
    }
    data <- x@tables$step1$res.kegg[[1]]
    data <- dplyr::filter(data, ID %in% !!ref)
    data <- .set_lab(data, sig(x), "Co-enriched KEGG pathway data")
    p.kegg <- plot_kegg(data, maxShow = 30, use = "pvalue")
    p.kegg <- .set_lab(p.kegg, sig(x), "Co-enriched KEGG pathway")
    .append_heading("代谢物与基因共同富集的 KEGG 通路")
    return(namel(data, p.kegg))
  })
