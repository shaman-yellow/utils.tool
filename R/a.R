
gs <- function(x, pattern, replace, ...) {
  gsub(pattern, replace, x, ...)
}

n.l <- function(name, object) {
  object <- list(object)
  names(object) <- name
  return(object)
}

nl <- function(names, values, as.list = TRUE, ...) {
  .as_dic(values, names, as.list = as.list, ...)
}

#' @export obj.size
#' @aliases obj.size
#' @description \code{obj.size}: ...
#' @rdname utilites
obj.size <- function(x, ...) {
  format(object.size(x), units = "MB", ...)
}

#' @export .as_dic
#' @aliases .as_dic
#' @description \code{.as_dic}: ...
#' @rdname utilites
.as_dic <- function(vec, names, default,
  fill = TRUE, as.list = TRUE, na.rm = FALSE)
{
  if (is.null(names(vec))) {
    names(vec) <- names[seq_along(vec)]
  }
  if (fill) {
    if (any(!names %in% names(vec))) {
      ex.names <- names[!names %in% names(vec)]
      ex <- rep(default, length(ex.names))
      names(ex) <- ex.names
      vec <- c(vec, ex)
    }
  }
  if (as.list) {
    if (!is.list(vec))
      vec <- as.list(vec)
  }
  if (na.rm) {
    vec <- vec[!is.na(names(vec))]
  }
  vec
}

#' @export grouping_vec2list
#' @aliases grouping_vec2list
#' @description \code{grouping_vec2list}: ...
#' @rdname query_synonyms
grouping_vec2list <- function(vector, group_number, byrow = FALSE){
  if(length(vector) < group_number){
    attr(vector, "name") <- "G1"
    return(list(vector))
  }
  rest <- length(vector) %% group_number
  group <- matrix(vector[1:(length(vector) - rest)],
    ncol = group_number,
    byrow = byrow)
  group <- apply(group, 1, c, simplify = FALSE)
  group <- c(group, list(tail(vector, n = rest)))
  group <- lapply(seq_along(group),
    function(n) {
      vec <- group[[n]]
      attr(vec, "name") <- paste0("G", n)
      vec
    })
  if(rest == 0)
    group <- group[1:(length(group) - 1)]
  return(group)
}

