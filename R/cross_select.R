# ==========================================================================
# filter data according to colunm within another data.frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

select_features <- function(mcn, classes, q.value = .05, logfc = .3,
                            coef = NULL) {
  if (!requireNamespace("MCnebula2", quietly = T))
    stop("package 'MCnebula2' must be available.")
  .check_data(statistic_set(mcn), list(top_table = "binary_comparison"))
  .check_data(mcn, list(nebula_index = "create_nebula_index"))
  stat <- top_table(statistic_set(mcn))
  if (!is.null(coef)) {
    stat <- stat[coef]
  }
  stat <- data.frame(data.table::rbindlist(stat))
  data.lst <- list(nebula_index(mcn), stat)
  filter.lst <- list(list(rlang::quo(class.name %in% dplyr::all_of(classes))),
                     list(rlang::quo(adj.P.Val < q.value),
                          rlang::quo(abs(logFC) > logfc)))
  cross_select(data.lst, filter.lst, ".features_id", "class.name")
}

cross_select <- function(data.lst, filter.lst, target, split = NULL) {
  if (!is.list(data.lst) | !is.list(filter.lst))
    stop("`data.lst` and `filter.lst` must be 'list'.")
  if (length(data.lst) != length(filter.lst))
    stop("`data.lst` and `filter.lst` must be 'list' with the same length.")
  lst <- lapply(1:length(data.lst),
                function(n) {
                  if (!is.null(filter.lst[[n]]))
                    dplyr::filter(data.lst[[n]], !!!(filter.lst[[n]]))
                  else
                    data.lst[[n]]
                })
  fun <- function(res, lst) {
    for (i in 2:length(lst)) {
      res <- res[res %in% lst[[ i ]]]
    }
    return(res)
  }
  if (is.null(split)) {
    lst <- lapply(lst, function(data) data[[ target ]])
    res <- fun(lst[[1]], lst)
  } else {
    res <- lapply(split(lst[[1]], lst[[1]][[ split ]]),
                  function(data) data[[ target ]])
    lst <- lapply(lst, function(data) data[[ target ]])
    res <- lapply(res, fun, lst = lst)
  }
  return(res)
}
