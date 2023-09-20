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
    info = c("...")
    ))

setGeneric("asjob_fella", 
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
  function(x, org = "hsa", db_fella = "../"){
    step_message("Use KEGG ID to perform FELLA.
      "
    )
    db_dir <- init_fella(db_fella, org)
    if (is.null(x$db)) {
      db <- load_fella(db_dir)
    } else {
      db <- x$db
    }
    obj.lst <- enrich_fella(x$ids.lst, db)
    names(obj.lst) <- names(x$ids.lst)
    x$graph.lst <- graph_fella(obj.lst, db)
    p.enrich <- lapply(x$graph.lst, plotGraph_fella)
    x@plots[[ 1 ]] <- namel(p.enrich)
    x$org <- org
    x$db_dir <- db_dir
    x$db <- db
    return(x)
  })


