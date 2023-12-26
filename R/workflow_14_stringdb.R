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
    method = "R package `STEINGdb` used for PPI network construction"
    ))

job_stringdb <- function(data)
{
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
    job_stringdb(data.frame(hgnc_symbol = x))
  })

setMethod("step0", signature = c(x = "job_stringdb"),
  function(x){
    step_message("Prepare your data with function `job_stringdb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_stringdb"),
  function(x, tops = 30, mcc_layout = "circle", layout = "kk", species = 9606,
    network_type = "phy", input_directory = "../", version = "11.5", label = F)
  {
    step_message("Create PPI network.
      "
    )
    if (is.null(x@params$sdb)) {
      sdb <- new_stringdb(species = species, network_type = network_type,
        input_directory = input_directory, version = version
      )
      x@params$sdb <- sdb
    } else {
      sdb <- x@params$sdb 
    }
    if (is.null(x@params$res.str)) {
      res.str <- create_interGraph(sdb, data.frame(object(x)), col = "hgnc_symbol")
      x@params$res.str <- res.str
    } else {
      res.str <- x@params$res.str
    }
    if (is.null(x@params$graph)) {
      graph <- fast_layout.str(res.str, sdb, layout = layout)
      x@params$graph <- graph
    } else {
      graph <- x@params$graph 
    }
    edges <- as_tibble(igraph::as_data_frame(res.str$graph))
    edges <- map(edges, "from", res.str$mapped, "STRING_id", "hgnc_symbol", rename = F)
    edges <- map(edges, "to", res.str$mapped, "STRING_id", "hgnc_symbol", rename = F)
    x$edges <- edges
    p.ppi <- plot_network.str(graph, label = label)
    ## hub genes
    hub_genes <- cal_mcc.str(res.str, "hgnc_symbol", F)
    graph_mcc <- get_subgraph.mcc(res.str$graph, hub_genes, top = tops)
    graph_mcc <- fast_layout(graph_mcc, layout = mcc_layout)
    p.mcc <- plot_networkFill.str(graph_mcc, label = "hgnc_symbol")
    x@plots[[1]] <- namel(p.ppi, p.mcc)
    x@tables[[1]] <- c(list(mapped = relocate(res.str$mapped, hgnc_symbol, STRING_id)),
      namel(hub_genes)
    )
    x@params$tops <- tops
    return(x)
  })

setMethod("asjob_enrich", signature = c(x = "job_stringdb"),
  function(x, tops = x@params$tops){
    ids <- head(x@tables$step1$hub_genes$hgnc_symbol, tops)
    job_enrich(list(hub_genes = ids), x@tables$step1$hub_genes)
  })

new_stringdb <- function(
  score_threshold = 200,
  species = 9606,
  network_type = c("physical", "full"),
  input_directory = "../",
  version = "11.5")
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

fast_layout.str <- function(res.str, sdb, layout = "fr", seed = 10) {
  ## get annotation
  target <- paste0(sdb$input_directory, "/", sdb$species, ".protein.info.v",
    sdb$version, ".txt")
  dir.create(dir <- paste0(sdb$input_directory, "/temp"), F)
  if (!file.exists(file <- paste0(dir, "/", get_filename(target)))) {
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
  arr.len = 2, edge.color = 'lightblue', edge.width = 1, lab.fill = "MCC score",
  label = "genes")
{
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      # arrow = arrow(length = unit(arr.len, 'mm')),
      color = edge.color, width = edge.width) +
    geom_node_point(aes(x = x, y = y, fill = ifelse(is.na(MCC_score), 0, MCC_score)),
      size = node.size, shape = 21, stroke = .3) +
    geom_node_text(aes(label = !!rlang::sym(label)), size = label.size) +
    scale_fill_gradient(low = "lightyellow", high = "red") +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    labs(fill = "MCC score") +
    theme_void() +
    theme(plot.margin = margin(r = .05, unit = "npc")) +
    geom_blank()
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
