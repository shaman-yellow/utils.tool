# ==========================================================================
# workflow of markers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_markers <- setClass("job_markers", 
  contains = c("job"),
  prototype = prototype(
    pg = "markers",
    info = c(""),
    cite = "[@CellMarker_a_m_Zhang_2019]",
    method = "",
    tag = "markers",
    analysis = "CellMarker 细胞注释的Marker验证"
    ))

job_markers <- function(tissue = NULL, org = c("Human", "Mouse"), 
  db_dir = .prefix("cellMarker", "db"), distinct = TRUE)
{
  org <- match.arg(org)
  dir.create(db_dir, FALSE)
  url <- "http://xteam.xbio.top/CellMarker/download/all_cell_markers.txt"
  file <- file.path(db_dir, basename(url))
  if (!file.exists(file)) {
    download.file(url, file)
  }
  db <- ftibble(file)
  if (!is.null(tissue)) {
    db <- dplyr::filter(db, grpl(tissueType, !!tissue, TRUE))
  }
  db <- dplyr::filter(db, speciesType == org)
  db <- dplyr::relocate(db, cellType, cellName, geneSymbol)
  if (distinct) {
    db <- dplyr::arrange(
      db, cellName, dplyr::desc(nchar(geneSymbol))
    )
    db <- dplyr::distinct(db, cellName, .keep_all = TRUE)
  }
  if (!is.null(tissue)) {
    message(glue::glue("Get all types:\n{showStrings(db$tissueType, TRUE, TRUE)}"))
  }
  x <- .job_markers()
  x$org <- org
  x$tissue <- tissue
  x$db <- db
  return(x)
}

setMethod("step0", signature = c(x = "job_markers"),
  function(x){
    step_message("Prepare your data with function `job_markers`.")
  })

setMethod("step1", signature = c(x = "job_markers"),
  function(x, cells = NULL, basic = c(
      "^T cell$", "^B cell$", "^Dendritic cell$", "^Monocyte$", "^Natural killer cell$"
    ))
  {
    step_message("Try get Cell Markers by cell names (pattern).")
    if (!is.null(basic)) {
      message(glue::glue('!is.null(basic). Merge default cells: {bind(basic)}'))
      cells <- c(cells, basic)
    }
    data <- dplyr::filter(
      x$db, grpl(cellName, paste0(cells, collapse = "|"))
    )
    if (nrow(data) < length(cells)) {
      message(glue::glue('nrow(data) < length(cells). Not got enough data.'))
    }
    print(data)
    x$data <- data
    return(x)
  })
