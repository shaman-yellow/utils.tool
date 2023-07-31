# ==========================================================================
# heat map with ggplot2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases plot_heatmap
#'
#' @title Plot heat map with ggplot2
#'
#' @description According to list of 'ID' to draw mutiple heatmap...
#'
#' @name plot_heatmap
NULL
#> NULL

#' @export plot_heatmap
#' @aliases plot_heatmap
#' @description \code{plot_heatmap}: ...
#' @rdname plot_heatmap
plot_heatmap <- function(id.lst, data, metadata,
  pal_class = ggsci::pal_futurama()(12), pal_group,
  clust_row = T, clust_col = T, method = 'complete')
{
  if (is.null(names(id.lst))) {
    stop("is.null(names(id.lst)) == T. The names of `id.lst` should be chemical classes.")
  }
  if (is.null(names(pal_class))) {
    pal_class <- pal_class[1:length(id.lst)]
    names(pal_class) <- names(id.lst)
  }
  .check_columns(metadata, c("sample", "group"), "metadata")
  .check_columns(data, c(".features_id", "sample", "value"), "data")
  lst <- sapply(names(id.lst), simplify = F,
    function(class.name) {
      ## basic heatmap
      ids <- id.lst[[ class.name ]]
      data <- dplyr::filter(data, .data$.features_id %in% dplyr::all_of(ids))
      p <- tile_heatmap(data)
      ## chemical classes
      data.class <- data.frame(class = class.name, .features_id = ids)
      pal_class <- pal_class[names(pal_class) == class.name]
      p <- add_ygroup.tile.heatmap(data.class, p, pal_class)
      ## cluster tree
      if (clust_row | clust_col) {
        data.w <- tidyr::spread(data, .data$sample, .data$value)
        data.w <- data.frame(data.w)
        rownames(data.w) <- data.w$.features_id
        data.w <- dplyr::select(data.w, dplyr::all_of(metadata[[ "sample" ]]))
        p <- add_tree.heatmap(
          data.w, p, method = method,
          clust_row = clust_row, clust_col = clust_col
        )
      }
      ## sample metadata
      p <- add_xgroup.tile.heatmap(metadata, p, pal_group)
      return(p)
    })
  return(lst)
}

#' @export handling_na
#' @aliases handling_na
#' @description \code{handling_na}:
#' For each subset of data, the missing values will be filled with the average
#' value; if the set is all missing values, they will be filled with zero.
#' @rdname plot_heatmap
handling_na <- function(data, id.cols = c(".features_id"),
  metadata, sample.col = "sample", group.col = "group")
{
  metadata <- metadata[, c(sample.col, group.col)]
  metadata <- split(metadata, metadata[[ group.col ]])
  id.cols <- data[, id.cols]
  data <- lapply(names(metadata),
    function(group) {
      meta <- metadata[[ group ]]
      df <- data[, meta[[ sample.col ]]]
      lst <- apply(df, 1, simplify = F,
        function(vec) {
          if (all(is.na(vec))) {
            vec[] <- 0
          } else if (any(!is.na(vec))) {
            vec[is.na(vec)] <- mean(vec, na.rm = T)
          }
          dplyr::bind_rows(vec)
        })
      data.table::rbindlist(lst)
    })
  data <- do.call(dplyr::bind_cols, data)
  dplyr::bind_cols(id.cols, data)
}

#' @export log_trans
#' @aliases log_trans
#' @description \code{log_trans}:
#' Convert wide data to long data; log transform the values; if there is a
#' value 0, replace it with 1/10 of the minimum value of the value column.
#' @rdname plot_heatmap
log_trans <- function(data, id.cols = c(".features_id"),
  key = "sample", value = "value",
  set_min = T, factor = 10, fun = log2, center = T)
{
  data <- tidyr::gather(data, !!key, !!value, -dplyr::all_of(id.cols))
  if (set_min) {
    min <- min(dplyr::filter(data, .data[[ value ]] != 0)[[ value ]])
    data[[ value ]] <- ifelse(data[[ value ]] == 0, min / factor, data[[ value ]])
  }
  data[[ value ]] <- fun(data[[ value ]])
  if (center) {
    data[[ value ]] <- scale(data[[ value ]], scale = F)[, 1]
  }
  return(data)
}

# dot_heatmap <- function(df){
#   p <- ggplot(df, aes(x = sample, y = .features_id)) +
#     geom_point(aes(size = abs(value), color = value), shape = 16) +
#     theme_minimal() +
#     guides(size = "none") +
#     scale_color_gradient2(low = "#3182BDFF", high = "#A73030FF") +
#     theme(text = element_text(family = .font),
#       axis.text.x = element_text(angle = 90))
#     return(p)
# }

dot_heatmap <- function(data, x = "sample", y = ".features_id",
  color = "value", size = "value", shape = NULL,
  lab_x = "Sample", lab_y = "Feature ID", lab_color = "log2 (Feature level)",
  lab_size = "", lab_shape = NULL, ...)
{
  scale_color <- if ( 0L > min(data[[ color ]]) & 0L < max(data[[ color ]])) {
    scale_color_gradient2(
      low = "#3182BDFF", high = "#A73030FF",
      limits = range(data[[ color ]])) 
  } else {
    scale_color_gradientn(
      colors = c("#3182BDFF", "white", "#A73030FF"),
      limits = range(data[[ color ]]))
  }
  if (is.null(shape)) {
    geom_point <- geom_point(aes(color = !!rlang::sym(color), size = !!rlang::sym(size)), shape = 16)
    guides <- geom_blank()
  } else {
    geom_point <- geom_point(aes(color = !!rlang::sym(color), size = !!rlang::sym(size),
        shape = !!rlang::sym(shape)))
    guides <- guides(shape = guide_legend(override.aes = list(size = 5)))
  }
  p <- ggplot(data, aes(x = !!rlang::sym(x), y = !!rlang::sym(y))) +
    geom_point +
    theme_minimal() +
    scale_color +
    scale_size(limits = range(data[[ size ]])) +
    labs(x = lab_x, y = lab_y, color = lab_color, size = lab_size, shape = lab_shape) +
    theme(text = element_text(family = .font, face = "bold"),
      axis.text = element_text(face = "plain"),
      axis.text.x = element_blank()) +
    guides
  return(p)
}

## long data
tile_heatmap <- function(data, x = "sample", y = ".features_id", fill = "value",
  lab_x = "Sample", lab_y = "Feature ID", lab_fill = "log2 (Feature level)", ...)
{
  scale_fill <- if ( 0L > min(data[[ fill ]]) & 0L < max(data[[ fill ]])) {
    scale_fill_gradient2(
      low = "#3182BDFF", high = "#A73030FF",
      limits = c(min(data[[ fill ]]), max(data[[ fill ]]))) 
  } else {
    scale_fill_gradientn(
      colors = c("#3182BDFF", "white", "#A73030FF"),
      limits = c(min(data[[ fill ]]), max(data[[ fill ]])))
  }
  p <- ggplot(data, aes(x = !!rlang::sym(x), y = !!rlang::sym(y))) +
    geom_tile(aes(fill = !!rlang::sym(fill)),
      color = "white", height = 1, width = 1, size = 0.2) +
    theme_minimal() +
    scale_fill +
    labs(x = lab_x, y = lab_y, fill = lab_fill) +
    theme(text = element_text(family = .font, face = "bold"),
      axis.text = element_text(face = "plain"),
      axis.text.x = element_blank()) +
    geom_blank()
  return(p)
}

plot_ytree <- function(data, method = "complete") {
  p <- hclust(dist(data), method)
  p <- ggtree::ggtree(p, layout = "rectangular", branch.length = "none", hang = 0) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  p
}

plot_xtree <- function(data, method = "complete") {
  p <- plot_ytree(t(data), method) +
    ggtree::layout_dendrogram()
  p
}

#' @export add_tree.heatmap
#' @aliases add_tree.heatmap
#' @description \code{add_tree.heatmap}: ...
#' @rdname plot_heatmap
add_tree.heatmap <- function(data, p, clust_row = T, clust_col = T, method = 'complete', ...){
    if (clust_row) {
      phr <- plot_ytree(data)
      p <- aplot::insert_left(p, phr, width = 0.3)
    }
    if (clust_col) {
      phc <- plot_xtree(data)
      p <- aplot::insert_top(p, phc, height = 0.3)
    }
    return(p)
  }

#' @export add_xgroup.heatmap
#' @aliases add_xgroup.heatmap
#' @description \code{add_xgroup.heatmap}: ...
#' @rdname plot_heatmap
add_xgroup.heatmap <- function(df, p){
    p.xgroup <- ggplot(df, aes(y = "Group", x = sample)) +
      geom_point(aes(color = group), size = 6) +
      ggsci::scale_color_simpsons() +
      labs(x = "", y = "", fill = "Group") +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
      com <- aplot::insert_bottom(p, p.xgroup, height = 0.05)
      return(com)
  }

add_xgroup.dot.heatmap <- function(data, p = NULL, pal = NA, x = "sample", y = "Group", color = "group",
    lab_x = "", lab_y = "", lab_color = "Group")
  {
    scale_color <- if (is.null(pal)) {
      ggsci::scale_color_simpsons()
    } else {
      scale_color_manual(values = pal)
    }
    p.xgroup <- ggplot(data, aes(y = !!y, x = !!rlang::sym(x))) +
      geom_point(aes(color = !!rlang::sym(color)), shape = 16) +
      scale_color +
      labs(x = lab_x, y = lab_y, color = lab_color) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_blank()
    if (is.null(p))
      return(p.xgroup)
    com <- aplot::insert_bottom(p, p.xgroup, height = 0.05)
    return(com)
  }

#' @export add_xgroup.tile.heatmap
#' @aliases add_xgroup.tile.heatmap
#' @description \code{add_xgroup.tile.heatmap}: ...
#' @rdname plot_heatmap
add_xgroup.tile.heatmap <- function(data, p = NULL, pal = NA, x = "sample", y = "Group", fill = "group",
    lab_x = "", lab_y = "", lab_fill = "Group")
  {
    scale_fill <- if (is.null(pal)) {
      ggsci::scale_fill_simpsons()
    } else {
      scale_fill_manual(values = pal)
    }
    p.xgroup <- ggplot(data, aes(y = !!y, x = !!rlang::sym(x))) +
      geom_tile(aes(fill = !!rlang::sym(fill)), 
        color = "white", height = 1, width = 1, size = 0.2) +
      scale_fill +
      labs(x = lab_x, y = lab_y, fill = lab_fill) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_blank()
    if (is.null(p))
      return(p.xgroup)
    com <- aplot::insert_bottom(p, p.xgroup, height = 0.05)
    return(com)
  }

add_ygroup.dot.heatmap <- function(data, p = NULL, pal = NULL,
  x = "Class", y = ".features_id", color = "class",
    lab_x = "", lab_y = "", lab_color = "From")
  {
    scale_color <- if (is.null(pal)) {
      ggsci::scale_color_simpsons()
    } else {
      scale_color_manual(values = pal)
    }
    p.ygroup <- ggplot(data, aes(x = !!x, y = !!rlang::sym(y))) +
      geom_point(aes(color = !!rlang::sym(color)), shape = 16, size = 5) +
      labs(x = lab_x, y = lab_y, color = lab_color) +
      scale_color +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_blank()
    if (is.null(p))
      return(p.ygroup)
    com <- aplot::insert_left(p, p.ygroup, width = 0.02) 
    return(com)
  }

#' @export add_ygroup.tile.heatmap
#' @aliases add_ygroup.tile.heatmap
#' @description \code{add_ygroup.tile.heatmap}: ...
#' @rdname plot_heatmap
add_ygroup.tile.heatmap <- function(data, p = NULL, pal = NULL, x = "Class",
  y = ".features_id", fill = "class",
    lab_x = "", lab_y = "", lab_fill = "From")
  {
    scale_fill <- if (is.null(pal)) {
      ggsci::scale_fill_simpsons()
    } else {
      scale_fill_manual(values = pal)
    }
    p.ygroup <- ggplot(data, aes(x = !!x, y = !!rlang::sym(y))) +
      geom_tile(aes(fill = !!rlang::sym(fill)),
        color = "white", height = 1, width = 1, size = 0.2) +
      labs(x = lab_x, y = lab_y, fill = lab_fill) +
      scale_fill +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      geom_blank()
    if (is.null(p))
      return(p.ygroup)
    com <- aplot::insert_left(p, p.ygroup, width = 0.02) 
    return(com)
  }
