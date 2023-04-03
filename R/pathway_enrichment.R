# ==========================================================================
# Combining multiple tools for pathway enrichment analysis.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases pathway_enrichment
#'
#' @title Perform pathway enrichment via package of 'FELLA'
#'
#' @description Pathway enrichment analysis was performed using KEGG ID
#' via package of 'FELLA'.
#' (Convert CID to KEGG ID using the 'MetaboAnalystR' package.
#' See <https://github.com/xia-lab/MetaboAnalystR> for installation.)
#'
#' @name pathway_enrichment
NULL
#> NULL

#' @export init_fella
#' @aliases init_fella
#' @description \code{init_fella}: ...
#' @seealso [FELLA::buildDataFromGraph()], [FELLA::buildGraphFromKEGGREST()]
#' @rdname pathway_enrichment
init_fella <- 
  function(dir, org = c("hsa", "mmu", "rno"), seed = 1, rebuild = F) {
    if (!file.exists(dir))
      stop("file.exists(dir) == F")
    dir <- paste0(dir, "/fella_pathway")
    dir.create(dir, F)
    org <- match.arg(org)
    db.dir <- paste0(dir, "/", org, ".db.dir")
    if (file.exists(db.dir) & !rebuild) {
      return(db.dir)
    } else {
      graph.file <- paste0(dir, "/", org, ".graph.Rdata")
      unlink(db.dir, T)
      set.seed(seed)
      graph <- FELLA::buildGraphFromKEGGREST(organism = org)
      save(graph, file = graph.file)
      FELLA::buildDataFromGraph(
        keggdata.graph = graph,
        databaseDir = db.dir, internalDir = FALSE,
        matrices = c("hypergeom", "diffusion", "pagerank"),
        normality = c("diffusion", "pagerank"),
        dampingFactor = 0.85, niter = 100)
    }
    return(db.dir)
  }

#' @export load_fella
#' @aliases load_fella
#' @description \code{load_fella}: ...
#' @rdname pathway_enrichment
load_fella <- function(dir) {
  if(!file.exists(dir)){
    stop("file.exists(dir) == F")
  }
  FELLA::loadKEGGdata(
    databaseDir = dir, internalDir = FALSE, 
    loadMatrix = c("hypergeom", "diffusion", "pagerank")
  )
}

#' @export enrich_fella
#' @aliases enrich_fella
#' @description \code{enrich_fella}: ...
#' @rdname pathway_enrichment
enrich_fella <- function(id.lst, data) {
  if (!is.list(id.lst)) {
    id.lst <- list(id.lst)
  }
  lapply(1:length(id.lst),
    function(n) {
      message("\n=========", "Enrichment:", n, "=========")
      id <- id.lst[[n]]
      res <- try(
        FELLA::enrich(
          id, data = data,
          method = FELLA::listMethods(),
          approx = "normality"
        )
      )
      if (inherits(res, "try-error"))
        return(NULL)
      else
        res
    })
}

#' @export graph_fella
#' @aliases graph_fella
#' @description \code{graph_fella}: ...
#' @rdname pathway_enrichment
graph_fella <- function( obj.lst, data, method = c("pagerank", "diffusion", "hypergeom"),
  threshold = .1)
{
  method <- match.arg(method)
  graph.lst <-
    lapply(obj.lst,
      function(obj) {
        if (is.null(obj))
          return()
        inmap <- FELLA::getInput(obj)
        graph <- FELLA::generateResultsGraph(
          object = obj,
          method = method,
          threshold = threshold,
          data = data
        )
        graph <- tidygraph::as_tbl_graph(graph)
        graph <- dplyr::select(graph, -entrez)
        graph <- dplyr::mutate(
          graph, NAME = vapply(NAME, function(c) c[1], ""),
          abbrev.name = stringr::str_trunc(NAME, 15),
          input = ifelse(input, "Input", "Others"),
          type = vapply(
            name, FUN.VALUE = "", USE.NAMES = F,
            function(str){
              str <- stringr::str_extract(str, "^[^[0-9]]{1,3}|\\.")
              str <- ifelse(nchar(str) > 1, "pathway", str)
              switch(
                str, pathway = "Pathway",
                M = "Module",
                "." = "Enzyme",
                C = "Compound",
                R = "Reaction")
            }))
      })
}

#' @import ggraph
#' @export plotGraph_fella
#' @aliases plotGraph_fella
#' @description \code{plotGraph_fella}: Draw the graph via
#' package of 'ggplot2'.
#' @rdname pathway_enrichment
plotGraph_fella <- function(
  graph, layout = "graphopt", seed = 1,
  shape = c(Input = 15, Others = 16),
  color = c(
    Pathway = "#E64B35FF",
    Module = "#E377C2",
    Enzyme = "#EFC000",
    Reaction = "#4DBBD5FF",
    Compound = "#00A087FF"),
  size = c(
    Pathway = 7,
    Module = 5,
    Enzyme = 6,
    Reaction = 5,
    Compound = 10))
{
  set.seed(seed)
  layout <- ggraph::create_layout(graph, layout = layout)
  ggraph(layout) + 
    geom_edge_fan(
      aes(edge_width = weight),
      color = "black",
      show.legend = F,
      end_cap = ggraph::circle(3, 'mm'),
      arrow = arrow(length = unit(2, 'mm'))) + 
    geom_node_point(
      aes(color = type,
        shape = input,
        size = type),
      stroke = 0.1) + 
    ggraph::geom_node_text(
      aes(label = stringr::str_wrap(
          abbrev.name, 15)),
      size = 3,
      family = .font,
      color = "black") +
    scale_shape_manual(values = shape) +
    scale_color_manual(values = color) +
    scale_size_manual(values = size) +
    scale_edge_width(range = c(0.1, 0.3)) + 
    guides(
      size = "none",
      shape = guide_legend(override.aes = list(size = 4)),
      color = guide_legend(override.aes = list(size = 4))) +
    labs(color = "Category", shape = "Type") +
    theme_void() +
    theme(text = element_text(family = .font))
}

#' @export cid.to.kegg
#' @aliases cid.to.kegg
#' @description \code{cid.to.kegg}: ...
#' @rdname pathway_enrichment
cid.to.kegg <- 
  function(cids){
    if (!requireNamespace("MetaboAnalystR", quietly = T)) {
      stop("package 'MetaboAnalystR' not available.",
        "See <https://github.com/xia-lab/MetaboAnalystR> for installation.")
    }
    obj <- MetaboAnalystR::InitDataObjects("conc", "msetora", F)
    obj <- MetaboAnalystR::Setup.MapData(obj, cids)
    obj <- MetaboAnalystR::CrossReferencing(obj, "pubchem")
    obj <- MetaboAnalystR::CreateMappingResultTable(obj)
    obj <- dplyr::as_tibble(obj$dataSet$map.table)
    dplyr::filter(obj, KEGG != "NA")
  }


