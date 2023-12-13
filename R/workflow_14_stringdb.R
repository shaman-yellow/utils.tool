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


