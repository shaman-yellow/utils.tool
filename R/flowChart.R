# ==========================================================================
# help with drawing flow chart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sortSugi <- function(values, dec = T) {
  num <- values
  i <- 1L
  order <- order(values, decreasing = dec)
  num[order] <- vapply(num[order],
    FUN.VALUE = integer(1),
    function(value) {
      x <- i
      i <<- i + 1L
      return(x)
    })
  num
}

as_network <- function(lst, layout = 'sugiyama', seed = 100)
{
  lst <- vapply(lst, function(ch) gsub(" ", "", ch), character(1))
  if (any(grepl("::", lst))) {
    lst <- unlist(lapply(lst,
        function(ch) {
          if (grepl("::", ch)) {
            ch <- strsplit(ch, "::")[[ 1 ]]
            ch <- c(paste0(ch, collapse = ":"), paste0(rev(ch), collapse = ":"))
          } else ch
        }))
  }
  lst <- strsplit(lst, ":")
  lst <- lapply(lst,
    function(ch) {
      ch <- strsplit(ch, ",")
      names(ch) <- c("from", "to", n(attr, length(ch) - 2))
      do.call(data.frame, ch)
    })
  data <- data.table::rbindlist(lst)
  data <- dplyr::mutate(data, from = Hmisc::capitalize(from), to = Hmisc::capitalize(to))
  nodes <- data.frame(name = unique(unlist(apply(data, 1, c, simplify = F), use.names = F)))
  if (!is.null(seed)) {
    data <- lapply(seed,
      function(s) {
        set.seed(s)
        fast_layout(data, layout, nodes)
      })
    if (length(seed) == 1)
      data <- data[[1]]
  }
  data
}

#' @import ggplot2
flowChart <- function(graph, scale.x = 1.2, scale.y = 1.2, node.size = 4,
  sc = 8, ec = 8, arr.len = 2, edge.color = 'lightblue', edge.width = 1)
{
  if (is(graph, "layout_tbl_graph")) {
    graphs <- list(graph)
    num <- 1L
  } else {
    graphs <- graph
    num <- 2L
  }
  p.lst <- lapply(graphs,
    function(graph) {
      data <- ggraph::get_edges()(graph)
      devX <- data$node1.x - data$node2.x
      devY <- data$node1.y - data$node2.y
      devCos <- devX / sqrt(devX ^ 2 + devY ^ 2)
      p <- ggraph(graph) +
        geom_node_label(aes(label = paste0(sortSugi(y), ". ", name)),
          size = node.size, label.padding = u(.5, lines)) +
        geom_edge_fan(
          aes(x = x, y = y,
            start_cap = circle(nchar(node1.name) * devCos + 5, 'mm'),
            end_cap = circle(nchar(node2.name) * devCos + 5, 'mm')),
          arrow = arrow(length = unit(arr.len, 'mm'), type = "closed"),
          color = edge.color, width = edge.width) +
        scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
        scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
        theme_void()
      p
    })
  if (num == 1L) {
    return(p.lst[[1]])
  } else {
    preview.gl(p.lst)
  }
  p.lst
}

preview.gl <- function(p.lst) {
  lst <- lapply(p.lst, as_grob)
  names(lst) <- n(p, length(lst))
  panel <- frame_col(fill_list(names(lst), 1), lst)
  legend <- sapply(names(lst), simplify = F, gtext, gp_arg = list(cex = 2))
  legend <- frame_col(fill_list(names(legend), 1), legend)
  frame <- frame_row(c(legend = 1, panel = 5), namel(panel, legend))
  setdev(width = 20)
  draw(frame)
}

zoRange <- function (x, factor) 
{
  x <- x[ !is.na(x) ]
  x <- range(x)
  ex <- abs(x[2] - x[1]) * (factor - 1)
  x[1] <- x[1] - ex
  x[2] <- x[2] + ex
  return(x)
}

ll <- function(str, sep = ":", link = ":") {
  strs <- strsplit(str, split = sep)[[ 1 ]]
  lst <- list()
  length(lst) <- length(strs) - 1
  for (i in 1:(length(strs) - 1)) {
    lst[[ i ]] <- paste0(strs[i], ":", strs[i + 1])
  }
  return(lst)
}
