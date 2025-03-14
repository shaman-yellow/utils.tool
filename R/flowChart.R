# ==========================================================================
# help with drawing flow chart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sortSugi <- function(values, dec = TRUE) {
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

as_network <- function(lst, layout = 'sugiyama', seed = 100, env = parent.frame(1))
{
  lst <- lapply(lst, function(x) glue::glue(x, .envir = env))
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
  if (any(isMulti <- lengths(lst) > 2)) {
    lst[isMulti] <- lapply(lst[isMulti],
      function(x) {
        x <- lapply(seq_along(x)[-1], 
          function(n) {
            c(x[n-1], x[n])
          })
        unlist(x)
      })
    lst <- unlist(lst)
    lst <- lapply(seq(1, length(lst), 2), 
      function(n) {
        c(lst[[n]], lst[[n + 1]])
      })
  }
  lst <- lapply(lst,
    function(ch) {
      ch <- strsplit(ch, ",")
      names(ch) <- c("from", "to", n(attr, length(ch) - 2))
      do.call(data.frame, ch)
    })
  data <- data.table::rbindlist(lst)
  data <- dplyr::mutate(data, from = Hmisc::capitalize(from), to = Hmisc::capitalize(to))
  nodes <- data.frame(name = unique(unlist(apply(data, 1, c, simplify = FALSE), use.names = FALSE)))
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
  sc = 6, ec = 6, consider.string.length = TRUE,
  arr.len = 2, edge.color = 'lightblue', edge.width = .5)
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
      if (consider.string.length) {
        devX <- data$node1.x - data$node2.x
        devY <- data$node1.y - data$node2.y
        devCos <- abs(devX) / sqrt(devX ^ 2 + devY ^ 2)
        devCos <- devCos ^ 2
        start_cap_sleng <- circle(nchar(data$node1.name) * 2.5 + sc, 'mm')
        start_cap <- circle(nchar(data$node1.name) * devCos + sc, 'mm')
        end_cap_sleng <- circle(nchar(data$node2.name) * 2.5 + ec, 'mm')
        end_cap <- circle(nchar(data$node2.name) * devCos + ec, 'mm')
      } else {
        start_cap_sleng <- start_cap <- circle(sc, 'mm')
        end_cap_sleng <- end_cap <- circle(ec, 'mm')
      }
      data <- .set_whether_pack(data)
      data <- dplyr::mutate(
        data, start_cap = start_cap, end_cap = end_cap,
        sc = circle(sc, 'mm'), ec = circle(ec, 'mm'),
        start_cap_sleng = start_cap_sleng, end_cap_sleng = end_cap_sleng,
        strength = ifelse(isMulti, (1 - abs(levelOutside)) * 3, abs(levelOutside) * 2)
      )
      data_edge_diagonal <- dplyr::filter(
        data, !isMulti | (abs(levelOutside) > 0) & (abs(levelOutside) < 1)
      )
      data_edge_bend <- dplyr::anti_join(data, data_edge_diagonal)
      data_edge_bend <- split(data_edge_bend, ~ isPack)
      if (!is.null(data_edge_bend[[ "TRUE" ]])) {
        geom_edge_bend_1 <- geom_edge_bend(data = data_edge_bend[[ "TRUE" ]],
          aes(x = x, y = y, start_cap = sc, end_cap = end_cap_sleng),
          arrow = arrow(length = unit(arr.len, 'mm'), type = "closed"),
          color = edge.color, width = edge.width, flipped = TRUE,
          strength = 1)
      } else {
        geom_edge_bend_1 <- geom_blank()
      }
      if (!is.null(data_edge_bend[[ "FALSE" ]])) {
        geom_edge_bend_2 <- geom_edge_bend(data = data_edge_bend[[ "FALSE" ]],
          aes(x = x, y = y, start_cap = start_cap_sleng, end_cap = ec),
          arrow = arrow(length = unit(arr.len, 'mm'), type = "closed"),
          color = edge.color, width = edge.width,
          strength = 1)
      } else {
        geom_edge_bend_2 <- geom_blank()
      }
      p <- ggraph(graph) +
        geom_node_label(aes(label = paste0(sortSugi(y), ". ", name)),
          size = node.size, label.padding = u(.5, lines)) +
        geom_edge_bend_1 +
        geom_edge_bend_2 +
        geom_edge_diagonal(
          data = data_edge_diagonal,
          aes(x = x, y = y, start_cap = sc, end_cap = ec),
          arrow = arrow(length = unit(arr.len, 'mm'), type = "closed"),
          color = edge.color,
          width = edge.width, flipped = FALSE,
          strength = data_edge_diagonal$strength) +
        guides(edge_alpha = "none") +
        scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
        scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
        theme_void()
      p
    })
  if (num == 1L) {
    return(wrap(p.lst[[1]], 10, min(10, max(p.lst[[1]]$data$y) * .6), showtext = TRUE))
  } else {
    preview.gl(p.lst)
  }
  p.lst
}

.set_whether_pack <- function(data) {
  .check_columns(
    data, c("node1.x", "node2.x", "node1.y", "node2.y"), "data"
  )
  nodes <- dplyr::bind_rows(
    dplyr::select(data, name = from, x = node1.x, y = node1.y),
    dplyr::select(data, name = to, x = node2.x, y = node2.y)
  )
  nodes <- dplyr::distinct(nodes)
  res <- lapply(seq_len(nrow(data)),
    function(n) {
      setFrom <- dplyr::filter(data, node1.y == data$node1.y[n], to == data$to[n])
      setTo <- dplyr::filter(data, node2.y == data$node2.y[n], from == data$from[n])
      isPack <- nrow(setFrom) > nrow(setTo)
      fun_dist <- function(x) (x[2] - x[1]) / 2
      if (isPack) {
        range <- range(setFrom$node1.x)
        distCenter <- data$node1.x[n] - mean(range)
        levelOutside <- distCenter / fun_dist(range)
        isMulti <- nrow(setFrom) > 3
      } else {
        range <- range(setTo$node2.x)
        distCenter <- data$node2.x[n] - mean(range)
        levelOutside <- distCenter / fun_dist(range)
        isMulti <- nrow(setTo) > 3
      }
      levelOutside <- ifelse(is.infinite(levelOutside) | is.nan(levelOutside), 0, levelOutside)
      namel(isPack, levelOutside, isMulti, distCenter)
    })
  dplyr::bind_cols(data, do.call(dplyr::bind_rows, res))
}

preview.gl <- function(p.lst) {
  lst <- lapply(p.lst, as_grob)
  names(lst) <- n(p, length(lst))
  panel <- frame_col(fill_list(names(lst), 1), lst)
  legend <- sapply(names(lst), simplify = FALSE, gtext, gp_arg = list(cex = 2))
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
