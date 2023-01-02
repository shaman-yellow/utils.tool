# ==========================================================================
# output compounds identification table 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rename_table <- 
  function(data, export_name = .export_name) {
    format_table(data, NULL, NULL, NULL, NULL, NULL, export_name)
  }

format_table <- 
  function(data, filter = .filter_format, arrange = .arrange_format,
           distinct = .distinct_format, mutate = .mutate_format,
           select = .select_format, export_name = .export_name) {
    if (!is.null(filter))
      data <- dplyr::filter(data, !!!filter)
    if (!is.null(arrange))
      data <- dplyr::arrange(data, !!!arrange)
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

.filter_format <- 
  list(quote(tani.score >= .5))

.arrange_format <- 
  list(quote(inchikey2d),
       quote(desc(tani.score))
  )

.distinct_format <- 
  list(quote(inchikey2d))

.mutate_format <- 
  list(mz = quote(round(mz, 4)),
       error.mass = quote(floor(error.mass * 10) / 10),
       tani.score = quote(floor(tani.score * 100) / 100),
       rt.min = quote(round(rt.secound / 60, 1))
  )

.select_format <- c("No.", "synonym", ".features_id", "mz", "error.mass",
                  "rt.min", "mol.formula", "adduct", "tani.score", "inchikey2d",
                  "class"
)

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
                  class = "Class"
)



