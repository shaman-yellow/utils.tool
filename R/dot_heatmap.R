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

dot_heatmap <- function(df){
  p <- ggplot(df, aes(x = sample, y = .features_id)) +
    geom_point(aes(size = abs(value), color = value), shape = 16) +
    theme_minimal() +
    guides(size = "none") +
    scale_color_gradient2(low = "#3182BDFF", high = "#A73030FF") +
    theme(text = element_text(family = .font),
      axis.text.x = element_text(angle = 90))
    return(p)
}

## long data
tile_heatmap <- 
  function(df){
    p <- ggplot(df, aes(x = sample, y = .features_id)) +
      geom_tile(aes(fill = value),
        color = "white", height = 1, width = 1, size = 0.2) +
      theme_minimal() +
      scale_fill_gradient2(low = "#3182BDFF", high = "#A73030FF",
        limits = c(min(df$value), max(df$value))) +
      labs(x = "Sample", y = "Feature ID", fill = "log2 (Feature level)") +
      theme(text = element_text(family = .font, face = "bold"),
        axis.text = element_text(face = "plain"),
        axis.text.x = element_blank()
      )
      return(p)
  }

#' @export add_tree.heatmap
#' @aliases add_tree.heatmap
#' @description \code{add_tree.heatmap}: ...
#' @rdname plot_heatmap
add_tree.heatmap <- 
  function(df, p, clust_row = T, clust_col = T, method = 'complete'){
    if (clust_row) {
      phr <- hclust(dist(df), method)
      phr <- ggtree::ggtree(phr, layout = "rectangular", branch.length = "branch.length") +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      p <- aplot::insert_left(p, phr, width = 0.3)
    }
    if (clust_col) {
      phc <- hclust(dist(t(df)), method)
      phc <- ggtree::ggtree(phc, layout = "rectangular", branch.length = "branch.length") +
        ggtree::layout_dendrogram() +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
      p <- aplot::insert_top(p, phc, height = 0.3)
    }
    return(p)
  }

#' @export add_xgroup.heatmap
#' @aliases add_xgroup.heatmap
#' @description \code{add_xgroup.heatmap}: ...
#' @rdname plot_heatmap
add_xgroup.heatmap <- 
  function(df, p){
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

#' @export add_xgroup.tile.heatmap
#' @aliases add_xgroup.tile.heatmap
#' @description \code{add_xgroup.tile.heatmap}: ...
#' @rdname plot_heatmap
add_xgroup.tile.heatmap <- 
  function(df, p, pal = NA){
    expr.pal <- ifelse(is.na(pal),
      'ggsci::scale_fill_simpsons()',
      'scale_fill_manual(values = pal)')
    p.xgroup <- ggplot(df, aes(y = "Group", x = sample)) +
      geom_tile(aes(fill = group), 
        color = "white", height = 1, width = 1, size = 0.2) +
      eval(parse(text = expr.pal)) +
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

#' @export add_ygroup.tile.heatmap
#' @aliases add_ygroup.tile.heatmap
#' @description \code{add_ygroup.tile.heatmap}: ...
#' @rdname plot_heatmap
add_ygroup.tile.heatmap <-
  function(df, p, pal = NA){
    expr.pal <- ifelse(is.na(pal),
      'ggsci::scale_fill_npg()',
      'scale_fill_manual(values = pal)')
    p.ygroup <- ggplot(df, aes(x = "Class", y = .features_id)) +
      geom_tile(aes(fill = class),
        color = "white", height = 1, width = 1, size = 0.2) +
      labs(x = "", y = "", fill = "From") +
      eval(parse(text = expr.pal)) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        text = element_text(family = .font, face = "bold"),
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      )
      com <- aplot::insert_left(p, p.ygroup, width = 0.02) 
      return(com)
  }
