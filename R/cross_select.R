# ==========================================================================
# filter data according to colunm within another data.frame
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases select_features
#'
#' @title Select 'features' for MCnebula2
#'
#' @description Select significant 'features' from MCnebula2 with
#' statistic results for downstream analysis of metabolomics.
#'
#' @name select_features
NULL
#> NULL

#' @export select_features
#' @aliases select_features
#' @description \code{select_features}: ...
#' @rdname select_features
select_features <- function(
  mcn, classes = unique(nebula_index(mcn)$class.name),
  q.value = .05, logfc = .3, coef = NULL, tani.score_cutoff = NULL,
  order_by_coef = NULL, togather = F) 
{
  if (!requireNamespace("MCnebula2", quietly = T))
    stop("package 'MCnebula2' must be available.")
  .check_data(statistic_set(mcn), list(top_table = "binary_comparison"))
  .check_data(mcn, list(nebula_index = "create_nebula_index",
      features_annotation = "create_features_annotation"))
  stat <- top_table(statistic_set(mcn))
  if (!is.null(coef)) {
    stat <- stat[coef]
  }
  stat <- data.frame(data.table::rbindlist(stat))
  data.lst <- list(nebula_index(mcn), stat)
  filter.lst <- list(
    rlang::quos(class.name %in% dplyr::all_of(classes)),
    rlang::quos(adj.P.Val < q.value, abs(logFC) > logfc)
  )
  if (!is.null(tani.score_cutoff)) {
    data.lst[[3]] <- features_annotation(mcn)
    filter.lst[[3]] <- rlang::quos(tani.score >= tani.score_cutoff)
  }
  res <- cross_select(data.lst, filter.lst, ".features_id", "class.name")
  if (!is.null(order_by_coef)) {
    ranks <- top_table(statistic_set(mcn))[[ order_by_coef ]]$.features_id
    res <- lapply(res,
      function(ids) {
        ranks[ ranks %in% ids ]
      })
  }
  if (togather) {
    res <- unlist(res, use.names = F)
    res <- ranks[ ranks %in% res ]
  }
  return(res)
}

#' @export cross_select
#' @aliases cross_select
#' @description \code{cross_select}: ...
#' @rdname select_features
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

#' @importFrom rlang as_label
.check_data <- 
  function(object, lst, tip = "(...)"){
    target <- rlang::as_label(substitute(object))
    mapply(lst, names(lst), FUN = function(value, name){
             obj <- match.fun(name)(object)
             if (is.null(obj)) {
               stop(paste0("is.null(", name, "(", target, ")) == T. ",
                           "use `", value, tip, "` previously."))
             }
             if (is.list(obj)) {
               if (length(obj) == 0) {
                 stop(paste0("length(", name, "(", target, ")) == 0. ",
                             "use `", value, tip, "` previously."))
               }
             }
           })
  }

