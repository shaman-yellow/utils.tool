# ==========================================================================
# workflow of garnett
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_garnett <- setClass("job_garnett", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://cole-trapnell-lab.github.io/garnett/docs/#1b-train-your-own-classifier")
    ))

setGeneric("asjob_garnett", 
  function(x, ...) standardGeneric("asjob_garnett"))

setMethod("asjob_garnett", signature = c(x = "job_seurat"),
  function(x){
    step_message("
      Transform as `cell_data_set` for train as classifier.
      "
    )
    x <- suppressMessages(asjob_monocle(x))
    .job_garnett(object = object(x))
  })

setMethod("step0", signature = c(x = "job_garnett"),
  function(x){
    step_message("Prepare your data with function `job_garnett`.")
  })

setMethod("step1", signature = c(x = "job_garnett"),
  function(x, marker_file, db = org.Hs.eg.db::org.Hs.eg.db){
    step_message("
      Check the marker file.
      "
    )
    object(x) <- e(monocle3::estimate_size_factors(object(x)))
    markers <- e(garnett::check_markers(object(x), marker_file,
        db = db, cds_gene_id_type = "SYMBOL",
        marker_file_gene_id_type = "SYMBOL"))
    markers <- tibble::as_tibble(markers)
    plot_markers <- e(garnett::plot_markers(markers))
    x@tables[[ 1 ]] <- namel(markers)
    x@plots[[ 1 ]] <- namel(plot_markers)
    x@params <- namel(marker_file, db)
    return(x)
  })

setMethod("step2", signature = c(x = "job_garnett"),
  function(x, seed = 1){
    step_message("Tran markers as  ")
    classifier <- e(garnett::train_cell_classifier(object(x),
        marker_file = x@params$marker_file,
        db = x@params$db,
        cds_gene_id_type = "SYMBOL",
        marker_file_gene_id_type = "SYMBOL"))
    x@params$classifier <- classifier
    return(x)
  })

as_marker_list <- function(x, filter, unique_level = 1:3) {
  .check_columns(x, c("cluster", "gene"), "x")
  stat <- table(x$gene)
  maker_list <- lapply(unique_level,
    function(level) {
      stat <- stat[ stat <= level ]
      x <- filter(x, gene %in% names(!!stat), gene %in% !!filter)
      x <- split(x$gene, x$cluster)
      x
    })
  .marker_list(level = unique_level, marker = maker_list)
}

as_marker_lines <- function(lst) {
  lst <- lapply(names(lst),
    function(name) {
      main <- paste0("expressed: ", paste0(lst[[ name ]], collapse = ", "))
      c(paste0(">", name), main, "")
    })
  unlist(lst)
}
