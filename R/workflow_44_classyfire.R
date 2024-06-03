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
      "Tools of `Classyfire` used for get systematic classification of chemical compounds"),
    tag = "cpd:class"
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
    x <- job_classyfire(allcpds$PubChem_id, type = "cid")
    x$herb_cpds <- allcpds
    return(x)
  })

job_classyfire <- function(query, type = NULL) {
  query <- unique(query)
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
    db_inchi <- frbind(extract_rdata_list(x$rd.inchi, object(x)))
    db_inchi <- dplyr::filter(db_inchi, !is.na(as.integer(CID)))
    x$db_inchi <- dplyr::mutate(db_inchi,
      inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]+"))
    return(x)
  })

setMethod("step2", signature = c(x = "job_classyfire"),
  function(x, pattern_match = NULL, cl = 5, savedir = "classify"){
    step_message("Use inchikey query classification.")
    if (x$type == "inchikey") {
      query <- unique(x$db_inchi$inchikey2d)
    } else if (x$type == "cid") {
      query <- unique(x$db_inchi$CID)
    }
    x$rd.class <- query_classification(query, savedir,
      inchikey.rdata = x$rd.inchi)
    if (x$type == "cid") {
      query <- unique(x$db_inchi$inchikey2d)
    }
    db_class <- frbind(extract_rdata_list(x$rd.class, query), idcol = "inchikey2d")
    if (x$type == "cid") {
      db_class <- map(db_class, "inchikey2d", x$db_inchi, "inchikey2d", "CID", col = "cid")
    }
    x$db_class <- db_class
    plot_fun <- function(data, pattern_match = NULL) {
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
      graph <- ggraph::create_layout(graph, 'partition', circular = T, weight = freq)
      p <- ggraph(graph) +
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
      if (!is.null(pattern_match)) {
        data <- dplyr::as_tibble(graph)
        data <- dplyr::filter(data, grpl(name, pattern_match, T))
        if (nrow(data)) {
          p <- p + ggrepel::geom_label_repel(data = data,
            aes(x = x, y = y, label = name), color = "red")
        }
      }
      wrap(p, 14, 12)
    }
    p.classes_freq <- plot_fun(x$db_class, pattern_match)
    p.classes_freq <- .set_lab(p.classes_freq, sig(x), "classification hierarchy")
    if (x$type == "inchikey") {
      inchikey2d <- object(x)
    } else if (x$type == "cid") {
      inchikey2d <- x$db_inchi$inchikey2d[match(as.integer(object(x)), as.integer(x$db_inchi$CID))]
      attr(object(x), "inchikey2d") <- inchikey2d
    }
    isClassified <- inchikey2d %in% unique(x$db_class$inchikey2d)
    p.isClassified <- new_pie(isClassified, title = "Compounds classify",
      fun_text = ggrepel::geom_label_repel)
    p.isClassified <- .set_lab(p.isClassified, sig(x), "Compounds classify")
    if (!is.null(pattern_match)) {
      t.match <- dplyr::filter(db_class, grpl(Classification, pattern_match, T))
      match <- unique(t.match$inchikey2d)
      if (x$type == "cid") {
        x$match <- list(cids = object(x)[inchikey2d %in% match],
          inchikey2d = inchikey2d[ inchikey2d %in% match ])
      } else {
        x$match <- list(inchikey2d = match)
      }
    } else {
      t.match <- NULL
    }
    x@plots[[ 2 ]] <- namel(p.classes_freq, p.isClassified)
    x@tables[[ 2 ]] <- namel(t.class = db_class, t.match)
    return(x)
  })
