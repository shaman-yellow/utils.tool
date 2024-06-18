# ==========================================================================
# workflow of stringdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_stringdb <- setClass("job_stringdb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html"),
    cite = "[@TheStringDataSzklar2021; @CytohubbaIdenChin2014]",
    method = "R package `STEINGdb` used for PPI network construction",
    tag = "ppi:stringdb",
    analysis = "STRINGdb PPI 分析"
    ))

job_stringdb <- function(data)
{
  if (is.character(data)) {
    data <- data.frame(Symbol = rm.no(data))
  }
  data <- dplyr::distinct(data, Symbol)
  .job_stringdb(object = data)
}

setGeneric("asjob_stringdb", 
  function(x, ...) standardGeneric("asjob_stringdb"))

setMethod("asjob_stringdb", signature = c(x = "job_herb"),
  function(x){
    job_stringdb(data = x@params$ppi_used)
  })

setMethod("asjob_stringdb", signature = c(x = "character"),
  function(x){
    message("`x` should be gene Symbols.")
    job_stringdb(data.frame(Symbol = x))
  })

setMethod("step0", signature = c(x = "job_stringdb"),
  function(x){
    step_message("Prepare your data with function `job_stringdb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_stringdb"),
  function(x, tops = 30, layout = "kk", species = 9606,
    network_type = "phy", input_directory = .prefix("stringdb_physical_v12.0", name = "db"),
    version = "12.0", label = F, HLs = NULL, use.anno = T,
    file_anno = .prefix("stringdb_physical_v12.0/9606.protein.physical.links.full.v12.0.txt.gz", name = "db"),
    filter.exp = 0, filter.text = 0)
  {
    step_message("Create PPI network.")
    if (is.null(x@params$sdb)) {
      message("Use STRINGdb network type of '", network_type, "'")
      sdb <- new_stringdb(species = species, network_type = network_type,
        input_directory = input_directory, version = version
      )
      x@params$sdb <- sdb
    } else {
      sdb <- x@params$sdb 
    }
    if (is.null(x@params$res.str)) {
      res.str <- create_interGraph(sdb, data.frame(object(x)), col = "Symbol")
      x@params$res.str <- res.str
    } else {
      res.str <- x@params$res.str
    }
    if (is.null(x@params$graph)) {
      graph <- fast_layout.str(res.str, sdb, layout = layout, seed = x$seed)
      x@params$graph <- graph
    } else {
      graph <- x@params$graph 
    }
    edges <- as_tibble(igraph::as_data_frame(res.str$graph))
    edges <- dplyr::distinct(edges, from, to, .keep_all = T)
    if (species == 9606) {
      message("`use.anno` not available for non hsa.")
      use.anno <- F
    }
    if (use.anno) {
      message("Get PPI annotation from:\n\t", file_anno)
      anno <- ftibble(file_anno)
      edges <- tbmerge(edges[, 1:2], anno, by.x = c("from", "to"), by.y = paste0("protein", 1:2),
        all.x = T, sort = F)
      if (filter.exp || filter.text) {
        edges <- dplyr::filter(edges, experiments >= !!filter.exp, textmining >= !!filter.text)
      }
    }
    edges <- map(edges, "from", res.str$mapped, "STRING_id", "Symbol", rename = F)
    edges <- map(edges, "to", res.str$mapped, "STRING_id", "Symbol", rename = F)
    des_edges <- list("STRINGdb network type:" = match.arg(network_type, c("physical", "full")))
    if (use.anno) {
      des_edges <- c(des_edges,
        list("Filter experiments score:" = paste0("At least score ", filter.exp),
          "Filter textmining score:" = paste0("At least score ", filter.text)
        ))
    }
    edges <- .set_lab(edges, sig(x), "PPI annotation")
    attr(edges, "lich") <- new_lich(des_edges)
    x$edges <- edges
    p.ppi <- plot_network.str(graph, label = label)
    p.ppi <- .set_lab(wrap(p.ppi, 4.5, 3), sig(x), "raw PPI network")
    ## hub genes
    hub_genes <- cal_mcc.str(res.str, "Symbol", F)
    graph_mcc <- get_subgraph.mcc(res.str$graph, hub_genes, top = tops)
    x$graph_mcc <- graph_mcc <- fast_layout(graph_mcc, layout = "linear", circular = T)
    p.mcc <- plot_networkFill.str(graph_mcc, label = "Symbol", HLs = HLs)
    p.mcc <- .set_lab(wrap(p.mcc), sig(x), paste0("Top", tops, " MCC score"))
    x@plots[[1]] <- namel(p.ppi, p.mcc)
    x@tables[[1]] <- c(list(mapped = relocate(res.str$mapped, Symbol, STRING_id)),
      namel(hub_genes)
    )
    x@params$tops <- tops
    return(x)
  })

setMethod("filter", signature = c(x = "job_stringdb"),
  function(x, ref.x, ref.y, lab.x = "Source", lab.y = "Target",
    use = "preferred_name", data = x@params$graph, level.x = NULL,
    lab.fill = "log2FC",
    ## this top is used for 'from' or 'to'
    top = 10, use.top = c("from", "to"),
    top_in = NULL, keep.ref = T,
    arrow = T, ...)
  {
    message("Search and filter: ref.x in from, ref.y in to; or, reverse.")
    use.top <- match.arg(use.top)
    data <- tibble::as_tibble(get_edges()(data), .name_repair = "minimal")
    data <- dplyr::select(data, dplyr::ends_with(use))
    data <- dplyr::rename(data, from = 1, to = 2)
    if (keep.ref) {
      data <- dplyr::filter(data, from %in% c(ref.x, ref.y), to %in% c(ref.x, ref.y))
    } else {
      data <- dplyr::filter(data,
        (from %in% ref.x & to %in% ref.y) | (from %in% ref.y & to %in% ref.x)
      )  
    }
    data <- dplyr::mutate(data, needRev = ifelse(from %in% ref.x, F, T))
    data <- apply(data, 1,
      function(x) {
        if (x[[ "needRev" ]]) {
          x[2:1]
        } else {
          x[1:2]
        }
      })
    data <- tibble::as_tibble(t(data))
    edges <- data <- dplyr::rename(data, from = 1, to = 2)
    nodes <- tibble::tibble(name = unique(c(data$from, data$to)))
    nodes <- dplyr::mutate(nodes,
      type = ifelse(name %in% ref.x, "from", "to")
    )
    data <- dplyr::mutate(data, id = "pseudo")
    data <- dplyr::relocate(data, id)
    p.ppi <- plot_network.pharm(data, edge_width = 1, ax2 = lab.x, ax3 = lab.y,
      ax2.level = level.x, lab.fill = lab.fill, ...)
    p.ppi <- .set_lab(p.ppi, sig(x), "filtered and formated PPI network")
    mcc <- cal_mcc(edges)
    nodes <- map(nodes, "name", mcc, "name", "MCC_score", col = "MCC_score")
    fun_tops <- function() {
      ## tops from ref.x or ref.y
      tops <- dplyr::filter(nodes, type == !!use.top)
      if (!is.null(top_in)) {
        fun <- function() {
          name <- switch(use.top, from = lab.x, to = lab.y)
          if (!is(top_in, "list")) {
            message("`top_in` is not 'list' with names, converted as 'list'.")
            top_in <- nl("Set", list(top_in))
          }
          ## venn plot
          new_venn(lst = c(nl(name, list(tops$name)), top_in))
        }
        p.top_in <- fun()
        p.top_in <- .set_lab(p.top_in, sig(x), "intersection with pre-filter data")
        if (length(p.top_in$ins)) {
          tops <- dplyr::filter(tops, name %in% !!unlist(top_in))
        } else {
          stop("length(p.top_in$ins) == 0, no features in the `top_in`.")
        }
      } else {
        p.top_in <- NULL
      }
      if (!is.null(top) && !keep.ref) {
        tops <- dplyr::slice_max(tops, MCC_score, n = top)
        edges <- dplyr::filter(edges, !!rlang::sym(use.top) %in% tops$name)
      }
      nodes <- dplyr::slice(nodes, c(which(name %in% ref.x), which(name %in% ref.y)))
      nodes <- dplyr::distinct(nodes, name, .keep_all = T)
      graph <- igraph::graph_from_data_frame(edges, vertices = nodes)
      graph <- fast_layout(graph, layout = "linear", circular = T)
      p.mcc <- plot_networkFill.str(graph, label = "name",
        arrow = arrow, shape = T,
        levels = level.x,
        lab.fill = if (is.null(level.x)) "MCC score" else "Log2(FC)", ...
      )
      p.mcc <- .set_lab(p.mcc, sig(x), "Top MCC score")
      colnames(edges) <- c(lab.x, lab.y)
      namel(p.mcc, nodes, edges, p.top_in)
    }
    p.mcc <- fun_tops()
    namel(p.ppi, nodes, edges = edges, p.top_in = p.mcc$p.top_in,
      p.mcc = p.mcc$p.mcc, nodes_mcc = p.mcc$nodes, edges_mcc = p.mcc$edges
    )
  })

setMethod("asjob_enrich", signature = c(x = "job_stringdb"),
  function(x, tops = x@params$tops){
    ids <- head(x@tables$step1$hub_genes$Symbol, tops)
    job_enrich(list(hub_genes = ids), x@tables$step1$hub_genes)
  })

setMethod("vis", signature = c(x = "job_stringdb"),
  function(x, HLs){
    p <- ggraph(x@params$graph_mcc) +
      geom_edge_arc(
        aes(color = ifelse(node1.Symbol == !!HLs | node2.Symbol == !!HLs,
            paste0(node1.Symbol, " <-> ", node2.Symbol), "...Others"))) +
      geom_node_point(aes(x = x, y = y, fill = ifelse(is.na(MCC_score), 0, MCC_score)),
        size = 12, shape = 21, stroke = .3) +
      geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = Symbol,
          angle = -((-node_angle(x,  y) + 90) %% 180) + 90), size = 4) +
      scale_fill_gradient(low = "lightyellow", high = "red") +
      scale_edge_color_manual(values = c("grey92", color_set())) +
      labs(edge_color = "Link", fill = "MCC score") +
      theme_void()
    p <- wrap(p, 12, 9)
    p <- .set_lab(p, sig(x), "MCC score of PPI top feature")
    data <- get_edges()(x@params$graph_mcc)
    data <- dplyr::filter(data, node1.Symbol %in% !!HLs | node2.Symbol %in% !!HLs)
    data <- dplyr::select(data, node1.Symbol, node2.Symbol)
    namel(p.mcc = p, data)
  })

new_stringdb <- function(
  score_threshold = 200,
  species = 9606,
  network_type = c("physical", "full"),
  input_directory = .prefix("stringdb_physical_v12.0", name = "db"),
  version = "12.0")
{
  e(STRINGdb::STRINGdb$new(score_threshold = score_threshold,
    species = species, network_type = match.arg(network_type),
    input_directory = input_directory, version = version
  ))
}

create_interGraph <- function(sdb, data, col = "name", rm.na = T) {
  cli::cli_alert_info("sdb$map")
  mapped <- sdb$map(data.frame(data), col, removeUnmappedRows = rm.na)
  message()
  cli::cli_alert_info("sdb$get_subnetwork")
  graph <- sdb$get_subnetwork(mapped$STRING_id)
  list(mapped = mapped, graph = graph)
}

fast_layout.str <- function(res.str, sdb, layout = "fr", seed = runif(1, max = 1000000000)) {
  ## get annotation
  target <- paste0(sdb$input_directory, "/", sdb$species, ".protein.info.v",
    sdb$version, ".txt")
  dir.create(dir <- paste0(sdb$input_directory, "/temp"), F)
  if (!file.exists(file <- paste0(dir, "/", basename(target)))) {
    file.copy(gfile <- paste0(target, ".gz"), dir)
    R.utils::gunzip(paste0(file, ".gz"))
  }
  anno <- data.table::fread(file)
  ## merge with annotation
  igraph <- add_attr.igraph(res.str$graph, anno, by.y = "#string_protein_id")
  set.seed(seed)
  graph <- fast_layout(igraph, layout)
  attr(graph, "igraph") <- igraph 
  graph
}

cal_pagerank <- function(igraph) {
  res <- igraph::page_rank(igraph)$vector
  data <- data.frame(name = names(res), weight = unname(res))
  igraph <- add_attr.igraph(igraph, data, by.y = "name")
  igraph
}

get_nodes <- function(igraph, from = "vertices") {
  tibble::as_tibble(igraph::as_data_frame(igraph, from))
} 

dedup.edges <- function(igraph){
  add_attr.igraph(igraph)
}

add_attr.igraph <- function(igraph, data, by.x = "name", by.y, dedup.edges = T)
{
  comps <- igraph::as_data_frame(igraph, "both")
  if (!missing(data)) {
    nodes <- merge(comps$vertices, data, by.x = by.x, by.y = by.y, all.x = T)
  } else {
    nodes <- comps$vertices
  }
  if (dedup.edges) {
    edges <- dplyr::distinct(comps$edges)
  } else {
    edges <- comps$edges
  }
  igraph <- igraph::graph_from_data_frame(edges, vertices = nodes)
  igraph
}

output_graph <- function(igraph, file, format = "graphml", toCyDir = T) {
  igraph::write_graph(igraph, file, format = format)
}

plot_network.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, sc = 5, ec = 5, 
  arr.len = 2, edge.color = 'grey70', edge.width = .4, label = F)
{
  if (label) {
    layer.nodes <- geom_node_label(aes(label = preferred_name), size = label.size)
  } else {
    layer.nodes <- geom_node_point(aes(x = x, y = y, color = centrality_degree))
  }
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      color = edge.color, width = edge.width) +
    layer.nodes +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  p
} 

plot_networkFill.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, node.size = 12, sc = 5, ec = 5, 
  arr.len = 2, edge.color = NULL, edge.width = 1,
  lab.fill = if (is.null(levels)) "MCC score" else "Levels",
  label = "genes", HLs = NULL, arrow = F, shape = F, levels = NULL,
  label.shape = c(from = "from", to = "to"), netType = c("physical", "functional"), ...)
{
  dataNodes <- ggraph::get_nodes()(graph)
  if (is.null(levels)) {
    fill <- "MCC_score"
    dataNodes <- dplyr::mutate(dataNodes, Levels = ifelse(is.na(MCC_score), 0, MCC_score))
    pal <- rev(color_set2())
    scale_fill_gradient <- scale_fill_gradient(low = pal[1], high = pal[2])
  } else {
    fill <- "Levels"
    if (!is.data.frame(levels)) {
      stop("!is.data.frame(levels) == F")
    }
    message("The first and second columns of `levels` were used as name and levels.")
    dataNodes <- map(tibble::tibble(dataNodes), "name",
      levels, colnames(levels)[1], colnames(levels)[2], col = "Levels")
    dataNodes <- dplyr::mutate(dataNodes, Levels = ifelse(is.na(Levels), 0, Levels))
    pal <- color_set2()
    scale_fill_gradient <- ggplot2::scale_fill_gradient2(low = pal[2], high = pal[1])
  }
  if (shape) {
    geom_node_point <- geom_node_point(
      data = dataNodes,
      aes(x = x, y = y, fill = !!rlang::sym(fill), shape = type),
      size = node.size, stroke = .3, alpha = .7)
  } else {
    geom_node_point <- geom_node_point(
      data = dataNodes,
      aes(x = x, y = y, fill = !!rlang::sym(fill)),
      size = node.size, shape = 21, stroke = .3, alpha = .7)
  }
  if (is.null(edge.color)) {
    edge.color <- sample(color_set()[1:10], 1)
  }
  netType <- match.arg(netType)
  p <- ggraph(graph) +
    ggraph::geom_edge_arc(
      aes(x = x, y = y, edge_linetype = !!netType),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      arrow = if (arrow) arrow(length = unit(arr.len, 'mm')) else NULL,
      color = edge.color, width = edge.width, alpha = .5) +
    geom_node_point +
    geom_node_text(aes(label = !!rlang::sym(label)), size = label.size) +
    scale_fill_gradient +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    labs(fill = lab.fill, shape = "Type", edge_linetype = "Interaction") +
    theme_void() +
    theme(plot.margin = margin(r = .05, unit = "npc")) +
    geom_blank()
  if (shape) {
    p <- p + scale_shape_manual(values = c(24, 21, 22, 23), labels = label.shape)
  }
  if (!is.null(HLs)) {
    data <- dplyr::filter(dataNodes, name %in% !!HLs)
    p <- p + geom_point(data = data, aes(x = x, y = y), shape = 21, color = "red", size = 20)
  }
  p
} 

cal_mcc.str <- function(res.str, name = "name", rename = T){
  hubs_score <- cal_mcc(res.str$graph)
  hubs_score <- tbmerge(res.str$mapped, hubs_score, by.x = "STRING_id", by.y = "name", all.x = T)
  hubs_score <- dplyr::relocate(hubs_score, !!rlang::sym(name), MCC_score)
  hubs_score <- dplyr::arrange(hubs_score, dplyr::desc(MCC_score))
  if (rename)
    hubs_score <- dplyr::rename(hubs_score, genes = name)
  hubs_score
}

cal_mcc <- function(edges) 
{
  if (is(edges, "igraph")) {
    igraph <- edges
    edges <- igraph::as_data_frame(edges, "edges")
  } else if (is(edges, "data.frame")) {
    igraph <- igraph::graph_from_data_frame(edges, F)
  }
  nodes <- unique(unlist(c(edges[, 1], edges[ , 2])))
  maxCliques <- igraph::max_cliques(igraph)
  scores <- vapply(nodes, FUN.VALUE = double(1), USE.NAMES = F,
    function(node) {
      if.contains <- vapply(maxCliques, FUN.VALUE = logical(1), USE.NAMES = F,
        function(clique) {
          members <- attributes(clique)$names
          if (any(members == node)) T else F
        })
      in.cliques <- maxCliques[ if.contains ]
      scores <- vapply(in.cliques, FUN.VALUE = double(1),
        function(clique) {
          num <- length(attributes(clique)$names)
          factorial(num - 1)
        })
      sum(scores)
    })
  res <- data.frame(name = nodes, MCC_score = scores)
  res <- tibble::as_tibble(dplyr::arrange(res, dplyr::desc(MCC_score)))
  if (exists(".add_internal_job")) {
    .add_internal_job(
      .job(method = "The MCC score was calculated referring to algorithm of `CytoHubba`",
        cite = "[@CytohubbaIdenChin2014]"))
  }
  res
}

get_subgraph.mcc <- function(igraph, resMcc, top = 10) 
{
  tops <- dplyr::arrange(resMcc, dplyr::desc(MCC_score))
  tops <- head(tops$STRING_id, n = top)
  data <- igraph::as_data_frame(igraph, "both")
  nodes <- dplyr::filter(data$vertices, name %in% !!tops)
  nodes <- merge(nodes, resMcc, by.x = "name", by.y = "STRING_id", all.x = T)
  edges <- dplyr::filter(data$edges, (from %in% !!tops) & (to %in% !!tops))
  nodes <- dplyr::distinct(nodes, name, .keep_all = T)
  igraph <- igraph::graph_from_data_frame(edges, F, nodes)
  igraph <- dedup.edges(igraph)
  igraph
}

sortDup_edges <- function(edges) {
  edges.sort <- apply(dplyr::select(edges, 1:2), 1,
    function(vec) {
      sort(vec)
    })
  edges.sort <- tibble::as_tibble(data.frame(t(edges.sort)))
  edges.sort <- dplyr::distinct(edges.sort)
  edges.sort
}

getBelong_edges <- function(edges) {
  nodes <- unique(unlist(c(edges[, 1], edges[, 2])))
  edges.rev <- edges
  edges.rev[, 1:2] <- edges.rev[, 2:1]
  links.db <- rbind(edges, edges.rev)
  lst.belong <- split(links.db, unlist(links.db[, 1]))
  lst.belong <- lapply(lst.belong,
    function(data) unlist(data[, 2], use.names = F))
  lst.belong
}
