# ==========================================================================
# workflow of seurat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.job_cellchat <- setClass("job_cellchat", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html"),
    cite = "[@InferenceAndAJinS2021]",
    method = "CellChat used for cell communication analysis"
    ))

setGeneric("asjob_cellchat", 
  function(x, ...) standardGeneric("asjob_cellchat"))

setMethod("asjob_cellchat", signature = c(x = "job_seurat"),
  function(x, group.by = x@params$group.by, assay = "SCT", ...){
    step_message("Prarameter red{{`group.by`}}, red{{`assay`}}, 
      and red{{`...`}} would passed to `CellChat::createCellChat`.",
      show_end = NULL
    )
    if (x@step < 2) {
      stop("x@step < 2. At least, preprocessed assay data should be ready.")
    }
    if (is.null(group.by))
      stop("is.null(group.by) == T")
    object(x)@meta.data[[ group.by ]] %<>% droplevels()
    object <- e(CellChat::createCellChat(
        object = object(x), group.by = group.by, assay = assay,
        ...
        ))
    .job_cellchat(object = object, params = list(group.by = group.by))
  })

setMethod("step0", signature = c(x = "job_cellchat"),
  function(x){
    step_message("Prepare your data with methods `asjob_cellchat`. ",
      "A processed 'Seurat' object is needed."
    )
  })

setMethod("step1", signature = c(x = "job_cellchat"),
  function(x, db = CellChat::CellChatDB.human,
    ppi = CellChat::PPI.human, cl = 4, python = "/usr/bin/python3")
  {
    step_message("One step forward computation of most.
      yellow{{In a nutshell, this do:
      Infer cell communication;
      Infer signaling role of cell groups;
      Clustering signaling.}}
      By default, red{{`CellChat::CellChatDB.human`}}
      grey{{(Consider use `CellChat::subsetDB` to filter `db`)}}
      was used as reference database.
      The `future::plan('multisession', worker = cl)` were called.
      By default, red{{`CellChat::PPI.human`}} were used by `CellChat::projectData`
      to smooth gene expression level (set NULL to cancel).
      "
    )
    if (!is.null(python)) {
      e(base::Sys.setenv(RETICULATE_PYTHON = python))
      e(reticulate::py_config())
    }
    p.showdb <- e(CellChat::showDatabaseCategory(db))
    p.showdb <- wrap(p.showdb, 8, 4)
    object(x)@DB <- db
    future::plan("multisession", workers = cl)
    object(x) <- e(CellChat::subsetData(object(x)))
    ## cell communication
    object(x) <- e(CellChat::identifyOverExpressedGenes(object(x)))
    object(x) <- e(CellChat::identifyOverExpressedInteractions(object(x)))
    if (!is.null(ppi)) {
      object(x) <- e(CellChat::projectData(object(x), ppi))
    }
    object(x) <- e(CellChat::computeCommunProb(object(x)))
    object(x) <- e(CellChat::filterCommunication(object(x), min.cells = 10))
    lp_net <- as_tibble(CellChat::subsetCommunication(object(x)))
    object(x) <- e(CellChat::computeCommunProbPathway(object(x)))
    pathway_net <- as_tibble(CellChat::subsetCommunication(object(x), slot.name = "netP"))
    object(x) <- e(CellChat::aggregateNet(object(x)))
    p.comms <- plot_communication.cellchat(object(x))
    p.comms <- .set_lab(p.comms, sig(x), paste("overall communication", c("count", "weight", "individuals")))
    ## Signaling role of cell groups
    object(x) <- e(CellChat::netAnalysis_computeCentrality(object(x), slot.name = "netP"))
    ## Clustering
    object <- object(x)
    res <- try(
      e({
        for (i in c("functional", "structural")) {
          object(x) <- CellChat::computeNetSimilarity(object(x), type = i)
          object(x) <- CellChat::netEmbedding(object(x), type = i)
          object(x) <- CellChat::netClustering(object(x), type = i, do.parallel = F)
        }
      })
    )
    if (inherits(res, "try-error")) {
      message("Due to error, escape from clustering; But the object was returned.")
      object(x) <- object
    }
    x@plots[[ 1 ]] <- c(namel(p.showdb), p.comms)
    x@tables[[ 1 ]] <- namel(lp_net, pathway_net)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cellchat"),
  function(x, pathways = NULL){
    step_message("This step visualize all results computed in previous step.
      These are:
      1. Cells yellow{{communication}} in heatmap grey{{(and L-R contribution)}};
      2. L-R in cells yellow{{communication}} in bubble;
      3. Gene yellow{{expression}} of signaling in violin;
      4. Overview of yellow{{Signaling roles}} of composition in heatmap;
      5. Weight of yellow{{signaling roles}} grey{{('outgoing' or 'incomming')}} in scatter;
      6. L-R in yellow{{signaling roles}} of cells in heatmap.
      "
    )
    if (is.null(pathways)) {
      pathways <- object(x)@netP$pathways
    }
    sapply <- pbapply::pbsapply
    cell_comm_heatmap <- e(sapply(c("ALL", pathways), simplify = F,
      function(name) {
        if (name == "ALL")
          name <- NULL
        main <- try(CellChat::netVisual_heatmap(
          object(x), color.heatmap = "Reds", signaling = name), T)
        if (inherits(main, "try-error")) {
          return(list(main = NULL, contri = NULL))
        }
        if (!is.null(name)) {
          contri <- CellChat::extractEnrichedLR(object(x),
            signaling = name, geneLR.return = F)
          return(namel(main, contri))
        } else {
          return(namel(wrap(main)))
        }
      }))
    lr_comm_bubble <- e(CellChat::netVisual_bubble(object(x), remove.isolate = FALSE))
    gene_expr_violin <- e(sapply(pathways, simplify = F,
      function(name) {
        CellChat::plotGeneExpression(object(x),
            signaling = name, group.by = NULL) +
          theme(legend.position = "none")
      }))
    role_comps_heatmap <- e(sapply(pathways, simplify = F,
        function(name) {
          res <- try({CellChat::netAnalysis_signalingRole_network(object(x), signaling = name,
            width = 8, height = 2.5, font.size = 8, cluster.rows = T)
            recordPlot()}, T)
          if (inherits(res, "try-error"))
            return(NULL)
          else res
        }))
    role_weight_scatter <- e(sapply(pathways, simplify = F,
        function(name) {
          CellChat::netAnalysis_signalingRole_scatter(object(x),
            signaling = name)
        }))
    res <- try(lr_role_heatmap <- e(sapply(c("outgoing", "incoming", "all"), simplify = F,
          function(name) {
            p <- CellChat::netAnalysis_signalingRole_heatmap(object(x), pattern = name,
              height = 1 + length(object(x)@netP$pathways) * .35
            )
            wrap(grid::grid.grabExpr(print(p)))
          })))
    if (inherits(res, "try-error")) {
      lr_role_heatmap <- NULL
      message("Due to error, escape from `CellChat::netAnalysis_signalingRole_heatmap`; ",
        "But the object was returned.")
    }
    lr_role_heatmap <- .set_lab(lr_role_heatmap, sig(x), names(lr_role_heatmap), "ligand-receptor role")
    cell_comm_heatmap$ALL$main <- .set_lab(cell_comm_heatmap$ALL$main, sig(x), "Cell communication heatmap")
    x@plots[[ 2 ]] <- namel(cell_comm_heatmap, lr_comm_bubble, gene_expr_violin,
      role_comps_heatmap, role_weight_scatter, lr_role_heatmap)
    return(x)
  })

setMethod("step3", signature = c(x = "job_cellchat"),
  function(x){
    step_message("Select pattern number for identify communication patterns.")
    require(NMF)
    lst <- e(sapply(c("outgoing", "incoming"), 
      function(pattern) {
        p <- CellChat::selectK(object(x), pattern = pattern)
        wrap(p + theme(legend.position = "none"), 8, 4)
      }))
    x@plots[[ 3 ]] <- lst
    return(x)
  })

setMethod("step4", signature = c(x = "job_cellchat"),
  function(x, k.out, k.in){
    step_message("Identification of major signals for specific
      cell groups and general communication
      patterns.
      "
    )
    if (missing(k.out) | missing(k.in))
      stop("missing(k.out) | missing(k.in)")
    object(x) <- e(CellChat::identifyCommunicationPatterns(
        object(x), pattern = "outgoing", k = k.out, heatmap.show = F
        ))
    object(x) <- CellChat::identifyCommunicationPatterns(
      object(x), pattern = "incoming", k = k.in, heatmap.show = F
    )
    require(ggalluvial)
    lst.p <- e(sapply(c("outgoing", "incoming"), simplify = F,
        function(pattern) {
          p.alluvial <- CellChat::netAnalysis_river(object(x), pattern = pattern)
          p.alluvial <- wrap(p.alluvial, 12, 7)
          p.dot <- CellChat::netAnalysis_dot(object(x), pattern = pattern)
          p.dot <- wrap(p.dot, 7, 5)
          namel(p.alluvial, p.dot)
        }))
    x@plots[[ 4 ]] <- lst.p
    return(x)
  })

plot_communication.cellchat <- function(x) {
  groupSize <- as.integer(table(x@idents))
  p.aggre_count <- CellChat::netVisual_circle(x@net$count, vertex.weight = groupSize,
    weight.scale = T, label.edge = F, title.name = "Number of interactions")
  p.aggre_weight <- CellChat::netVisual_circle(x@net$weight, vertex.weight = groupSize,
    weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")
  mat <- x@net$weight
  scale <- cal_panelScale(length(levels(x@idents)))
  par(mfrow = scale, xpd = TRUE, omi = rep(0, 4), mar = rep(1.5, 4))
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    CellChat::netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T,
      edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
  p.commSep <- recordPlot()
  p.commSep <- wrap(p.commSep, scale[1] * 4, scale[2] * 4)
  namel(p.aggre_count, p.aggre_weight, p.commSep)
}

cal_panelScale <- function(num) {
  ncol <- round(sqrt(num))
  if (ncol ^ 2 < num) {
    nrow <- ncol + 1
  } else {
    nrow <- ncol
  }
  c(nrow, ncol)
}

setMethod("map", signature = c(x = "job_cellchat", ref = "character"),
  function(x, ref, ref2, width = 20, height = 5, cap = 5, layout = "sugiyama"){
    data <- select_pathway(x, ref, "lps")
    data <- dplyr::filter(data, grepl(ref2, source) | grepl(ref2, target))
    p <- plot_lps_interaction(data, cap, layout)
    p <- wrap(p, width, height)
    p <- .set_lab(p, sig(x), "ligand-receptor of", paste(ref, "communicate with", ref2))
    namel(p, data)
  })

plot_lps_interaction <- function(edges, cap = 5, layout = "sugiyama") {
  nodes1 <- dplyr::select(edges, name = ligand, belong = source)
  nodes1 <- dplyr::mutate(nodes1, role = "ligand")
  nodes2 <- dplyr::select(edges, name = receptor, belong = target)
  nodes2 <- dplyr::mutate(nodes2, role = "receptor")
  nodes <- dplyr::bind_rows(nodes1, nodes2)
  nodes <- dplyr::distinct(nodes, name, role)
  nodes <- split_lapply_rbind(nodes, ~ name,
    function(data) {
      if (nrow(data) == 2) {
        dplyr::mutate(data[1, ], role = "ligand_or_receptor")
      } else {
        data
      }
    })
  edges <- dplyr::relocate(edges, ligand, receptor)
  graph <- fast_layout(edges, layout, nodes)
  plot_sc_interaction <- function(graph, sc = cap, ec = cap, arr.len = 1)
  {
    p <- ggraph(graph) +
      geom_edge_arc(aes(x = x, y = y, width = -log2(pval + .001), color = source),
        start_cap = circle(sc, 'mm'),
        end_cap = circle(ec, 'mm'),
        arrow = arrow(length = unit(arr.len, 'mm')), alpha = .5) +
      geom_node_point(
        aes(x = x, y = y, shape = role,
          size = centrality_degree, fill = centrality_degree),
        stroke = .3) +
      ggrepel::geom_text_repel(aes(x = x, y = y, label = name), angle = 90, size = 2) +
      guides(size = "none", shape = guide_legend(override.aes = list(size = 4))) +
      scale_edge_width(range = c(.5, 1)) +
      scale_size(range = c(3, 6)) +
      scale_shape_manual(values = 21:25) +
      scale_fill_gradient2(low = "blue", high = "red") +
      scale_color_manual(values = color_set()) +
      theme_void() +
      theme(legend.position = "right")
    p
  }
  p <- plot_sc_interaction(graph)
  p
}

setGeneric("select_pathway", 
  function(x, ...) standardGeneric("select_pathway"))

setMethod("select_pathway", signature = c(x = "job_cellchat"),
  function(x, pattern, get = c("pathways", "lps", "intersect")){
    get <- match.arg(get)
    if (x@step < 1) {
      stop("x@step != 1")
    }
    pathways <- filter(x@tables$step1$pathway_net,
      grepl(pattern, source) | grepl(pattern, target),
      pval < .05
    )
    if (get == "pathways") {
      return(pathways)
    }
    lps <- filter(x@tables$step1$lp_net,
      grepl(pattern, source) | grepl(pattern, target),
      pval < .05
    )
    if (get == "lps") {
      return(lps)
    }
    ints <- intersect(unique(pathways$pathway_name), unique(lps$pathway_name))
    if (get == "intersect") {
      return(ints)
    }
  })

unique.lps <- function(x) {
  unique(unlist(stringr::str_extract_all(x, "[^_]+")))
}
