# ==========================================================================
# output compounds identification table 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases format_table
#'
#' @title Format table via dplyr::*
#'
#' @description Format the data.frame via: \code{dplyr::filter}, \code{dplyr::arrange},
#' \code{dplyr::distinct}, \code{dplyr::mutate}, \code{dplyr::select},
#' \code{dplyr::rename}.
#' @param data data.frame. From \code{features_annotation(mcn)}.
#'
#' @name format_table
NULL
#> NULL

#' @export rename_table
#' @aliases rename_table
#' @description \code{rename_table}: ...
#' @rdname format_table
rename_table <- 
  function(data, export_name = .export_name) {
    format_table(data, NULL, NULL, NULL, NULL, NULL, export_name)
  }

#' @export format_table
#' @aliases format_table
#' @description \code{format_table}: ...
#' @rdname format_table
format_table <- 
  function(data, filter = .filter_format, arrange = .arrange_format,
    distinct = .distinct_format, mutate = .mutate_format,
    select = .select_format, export_name = .export_name) {
    if (!is.null(filter))
      data <- dplyr::filter(data, !!!filter)
    if (!is.null(arrange)) {
      if (is.null(data.frame(data)$arrange.rank))
        data <- dplyr::mutate(data, arrange.rank = NA)
      data <- dplyr::arrange(data, !!!arrange)
    }
    if (!is.null(distinct))
      data <- dplyr::distinct(data, !!!distinct, .keep_all = T)
    if (!is.null(mutate))
      data <- dplyr::mutate(data, !!!mutate)
    if (!is.null(select)) {
      select <- select[select %in% colnames(data)]
      if (!is.null(select))
        data <- dplyr::select(data, dplyr::all_of(select))
    }
    if (!is.null(export_name)) {
      export_name <- export_name[names(export_name) %in% colnames(data)]
      export_name <- as.list(turn_vector(export_name))
      data <- dplyr::rename(data, !!!export_name)
    }
    tibble::as_tibble(data)
  }

#' @export .filter_format
#' @aliases .filter_format
#' @description \code{.filter_format}: ...
#' @rdname format_table
.filter_format <- 
  list(quote(tani.score >= .5))

#' @export .arrange_format
#' @aliases .arrange_format
#' @description \code{.arrange_format}: ...
#' @rdname format_table
.arrange_format <- 
  list(
    quote(arrange.rank),
    quote(inchikey2d),
    quote(desc(tani.score))
  )

#' @export .distinct_format
#' @aliases .distinct_format
#' @description \code{.distinct_format}: ...
#' @rdname format_table
.distinct_format <- 
  list(quote(inchikey2d))

#' @export .mutate_format
#' @aliases .mutate_format
#' @description \code{.mutate_format}: ...
#' @rdname format_table
.mutate_format <- 
  list(mz = quote(round(mz, 4)),
    error.mass = quote(floor(error.mass * 10) / 10),
    tani.score = quote(floor(tani.score * 100) / 100),
    rt.min = quote(round(rt.secound / 60, 1))
  )

#' @export .select_format
#' @aliases .select_format
#' @description \code{.select_format}: ...
#' @rdname format_table
.select_format <- c("No.", "synonym", ".features_id", "mz", "error.mass",
  "rt.min", "mol.formula", "adduct", "tani.score", "inchikey2d",
  "class", "logFC", "P.Value", "adj.P.Val"
)

#' @export .export_name
#' @aliases .export_name
#' @description \code{.export_name}: ...
#' @rdname format_table
.export_name <- c(mz = "Precursor m/z",
  rt.min = "RT (min)",
  similarity = "Spectral similarity",
  tani.score = "Tanimoto similarity",
  rel.index = "Relative index",
  rel.int. = "Relative intensity",
  group = "Group",
  .features_id = "ID",
  mol.formula = "Formula",
  inchikey2d = "InChIKey planar",
  error.mass = "Mass error (ppm)",
  synonym = "Synonym",
  adduct = "Adduct",
  class = "Class",
  logFC = "log2(FC)",
  P.Value = "P-value",
  adj.P.Val = "Q-value"
)



