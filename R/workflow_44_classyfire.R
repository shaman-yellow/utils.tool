# ==========================================================================
# workflow of classyfire
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_classyfire <- setClass("job_classyfire", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@PubchemSubstanKimS2015; @ClassyfireAutDjoumb2016]",
    method = paste0("Database `PubChem` used for querying information (e.g., InChIKey, CID) of chemical compounds; ",
      "Tools of `Classyfire` used for get systematic classification of chemical compounds")
    ))

setGeneric("asjob_classyfire", 
  function(x, ...) standardGeneric("asjob_classyfire"))

setMethod("asjob_classyfire", signature = c(x = "job_herb"),
  function(x){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    cpdsID <- unique(x@tables$step2$compounds_targets$Ingredient_id)
    allcpds <- dplyr::filter(x@object$component,
      !is.na(PubChem_id), Ingredient_id %in% !!cpdsID
    )
    x <- .job_classyfire(object = allcpds$PubChem_id)
    x$herb_cpds <- allcpds
    return(x)
  })

job_classyfire <- function(query, type = NULL) {
  if (is.null(type)) {
    if (is.numeric(query)) {
      type <- "cid"
    } else if (is.character(query)) {
      type <- "inchikey"
    }
  }
  x <- .job_classyfire(object = query)
  x$type <- type
  return(x)
}

setMethod("step0", signature = c(x = "job_classyfire"),
  function(x){
    step_message("Prepare your data with function `asjob_classyfire`.")
  })

setMethod("step1", signature = c(x = "job_classyfire"),
  function(x, cl = 5, savedir = "inchikey")
  {
    step_message("Use inchikey2d or Pubchem CID query inchikey.")
    x$rd.inchi <- query_inchikey(object(x), savedir, type = x$type, curl_cl = cl)
    db_inchi <- frbind(extract_rdata_list(x$rd.inchi))
    x$db_inchi <- dplyr::mutate(db_inchi,
      inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]+"))
    return(x)
  })

setMethod("step2", signature = c(x = "job_classyfire"),
  function(x, cl = 5, savedir = "classify"){
    step_message("Use inchikey query classification.")
    x$rd.class <- query_classification(unique(x$db_inchi$inchikey2d), savedir,
      inchikey.rdata = x$rd.inchi)
    db_class <- frbind(extract_rdata_list(x$rd.class), idcol = T)
    db_class <- tbmerge(x$db_inchi, db_class,
      by.x = "inchikey2d", by.y = ".id", all.x = T, allow.cartesian = T)
    x$db_class <- db_class
    if (!is.null(x$herb_cpds)) {
      cpds <- dplyr::select(x$herb_cpds, PubChem_id, Ingredient_id, Ingredient_name)
      x$db_class <- tbmerge(cpds, x$db_class, by.x = "PubChem_id", by.y = "CID", allow.cartesian = T)
    }
    plot_fun <- function(data) {
      nodes <- dplyr::summarize(dplyr::group_by(data, Classification),
        freq = length(Classification))
      nodes <- dplyr::rename(nodes, name = Classification)
      edges <- dplyr::distinct(data, Classification, inchikey2d)
      edges <- split_lapply_rbind(edges, ~ inchikey2d,
        function(data) {
          data <- dplyr::select(data, Classification)
          data <- cbind(data[-nrow(data), ], data[-1, ])
          colnames(data) <- c("from", "to")
          data
        })
      graph <- tidygraph::tbl_graph(nodes, edges)
      p <- ggraph(graph, 'partition', circular = T, weight = freq) +
        geom_node_arc_bar(aes(fill = freq)) +
        scale_fill_gradientn(colors = color_set()[1:3]) +
        geom_text(aes(x = x, y = y,
            angle = ifelse(freq == max(freq), 0,
              node_angle(x, y) + sign(node_angle(x, y) - 180) * 90),
            label = ifelse(freq < fivenum(freq)[4],
              "", stringr::str_trunc(name, 20)))) +
        theme_minimal() +
        labs(fill = "Number") +
        ggtitle("Classification hierarchy") +
        theme(axis.text = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(size = 15, face = "bold", hjust = .5)) +
        geom_blank()
      wrap(p, 14, 12)
    }
    p.classes_freq <- plot_fun(x$db_class)
    p.classes_freq <- .set_lab(p.classes_freq, sig(x), "classification hierarchy")
    x@plots[[ 2 ]] <- namel(p.classes_freq)
    return(x)
  })
