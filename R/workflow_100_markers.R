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
    analysis = "CellMarker 细胞注释的Marker"
    ))

job_markers <- function(tissue = NULL, type = c("Normal cell", "Cancer cell"),
  use_numbering = FALSE, org = c("Human", "Mouse"),
  db_dir = .prefix("cellMarker", "db"), distinct = FALSE)
{
  type <- match.arg(type, several.ok = TRUE)
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
  db <- dplyr::relocate(
    db, cellType, cellName, geneSymbol, cellMarker
  )
  db <- dplyr::filter(db, cellType %in% !!type)
  if (!use_numbering) {
    db <- dplyr::filter(db, !grpl(cellName, "[0-9]"))
  }
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
    ), gather = TRUE, duplicated = FALSE, least = 2L)
  {
    step_message("Try get Cell Markers by cell names (pattern).")
    if (!is.null(basic)) {
      message(glue::glue('!is.null(basic). Merge default cells: {bind(basic)}'))
      cells <- sort(c(cells, basic))
    }
    data <- dplyr::filter(
      x$db, grpl(cellName, paste0(cells, collapse = "|"))
    )
    notGot <- lapply(cells,
      function(pat) {
        if (!any(grpl(data$cellName, pat, TRUE))) {
          pat
        } else NULL
      })
    notGot <- unlist(notGot)
    if (nrow(data) < length(cells)) {
      message(glue::glue('nrow(data) < length(cells). Not got enough data.'))
    }
    print(data)
    x$data <- data
    features <- strsplit(gs(x$data$geneSymbol, "\\[|\\]| ", ""), ",")
    features <- nl(x$data$cellName, features)
    if (gather) {
      features <- data.frame(
        cell = rep(names(features), lengths(features)),
        marker = unlist(features, use.names = FALSE)
      )
      features <- split(features$marker, features$cell)
    }
    if (duplicated) {
      features <- sort_by(features, lengths(features))
      features <- tibble::tibble(
        cell = rep(names(features), lengths(features)),
        marker = unlist(features, use.names = FALSE),
        n = unlist(lapply(lengths(features), seq_len))
      )
      features <- features[ !duplicated(features$marker) | features$n <= least, ]
      features <- split(features$marker, features$cell)
    }
    x$.feature <- features
    message(glue::glue("Not got markers of cells: {crayon::red(bind(notGot))}"))
    return(x)
  })

setMethod("ref", signature = c(x = "job_markers"),
  function(x, from = c("db", "data")){
    from <- match.arg(from)
    data <- x[[ from ]]
    if (is.null(data)) {
      stop('is.null(data)')
    }
    data <- dplyr::select(data, cell = cellName, markers = geneSymbol)
    data <- dplyr::mutate(
      data, markers = strsplit(gs(markers, "\\[|\\]| ", ""), ",")
    )
    data <- reframe_col(data, "markers", unlist)
    dplyr::distinct(data)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_markers"),
  function(x, ref, extra = NULL, extra.after = NULL, 
    max = 2, soft = TRUE, group.by = x$group.by,
    markers = feature(ref))
  {
    if (missing(markers) && ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    if (x@step < 3L) {
      stop('x@step < 3L.')
    }
    p.cellMarker <- .plot_marker_heatmap(
      x, markers, group.by, extra, extra.after, max = max, soft = soft
    )
    p.cellMarker <- set_lab_legend(
      wrap(p.cellMarker, 7, 7),
      paste(sig(x), "CellMarker Validation"),
      "为使用 CellMarker 数据库中的特异性 Marker 对单细胞注释结果的验证热图。"
    )
    x$p.cellMarker <- p.cellMarker
    return(x)
  })
