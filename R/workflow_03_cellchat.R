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
    info = c("Tutorial: https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html")
    ))

setGeneric("asjob_cellchat", 
  function(x, ...) standardGeneric("asjob_cellchat"))

setMethod("asjob_cellchat", signature = c(x = "job_seurat"),
  function(x, group.by = "SingleR_cell", assay = "SCT", ...){
    step_message("Prarameter red{{`group.by`}}, red{{`assay`}}, 
      and red{{`...`}} would passed to `CellChat::createCellChat`.",
      show_end = NULL
    )
    if (x@step < 2) {
      stop("x@step < 2. At least, preprocessed assay data should be ready.")
    }
    if (group.by == "SingleR_cell") {
      if (x@step < 4L)
        stop("`step4(x)` has not been performed.")
    }
    object <- e(CellChat::createCellChat(
        object = object(x), group.by = group.by, assay = assay,
        ...
        ))
    .job_cellchat(object = object)
  })

setMethod("step0", signature = c(x = "job_cellchat"),
  function(x){
    step_message("Prepare your data with methods `asjob_cellchat`. ",
      "A processed 'Seurat' object is needed."
    )
  })

setMethod("step1", signature = c(x = "job_cellchat"),
  function(x, db = CellChat::CellChatDB.human,
    ppi = CellChat::PPI.human, cl = 4)
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
    ## Signaling role of cell groups
    object(x) <- e(CellChat::netAnalysis_computeCentrality(object(x), slot.name = "netP"))
    ## Clustering
    e({
      for (i in c("functional", "structural")) {
        object(x) <- CellChat::computeNetSimilarity(object(x), type = i)
        object(x) <- CellChat::netEmbedding(object(x), type = i)
        object(x) <- CellChat::netClustering(object(x), type = i, do.parallel = F)
      }
    })
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
    cell_comm_heatmap <- e(sapply(c("ALL", pathways), simplify = F,
      function(name) {
        if (name == "ALL")
          name <- NULL
        main <- CellChat::netVisual_heatmap(
          object(x), color.heatmap = "Reds", signaling = name)
        if (!is.null(name)) {
          contri <- CellChat::extractEnrichedLR(object(x),
            signaling = name, geneLR.return = F)
          return(namel(main, contri))
        } else {
          return(namel(main))
        }
      }))
    lr_comm_bubble <- e(CellChat::netVisual_bubble(object(x), remove.isolate = FALSE))
    gene_expr_violin <- e(sapply(pathways, simplify = F,
      function(name) {
        CellChat::plotGeneExpression(object(x),
            signaling = name, group.by = NULL)
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
    lr_role_heatmap <- e(sapply(c("outgoing", "incoming", "all"), simplify = F,
      function(name) {
        p <- CellChat::netAnalysis_signalingRole_heatmap(object(x), pattern = name)
        grid::grid.grabExpr(print(p))
      }))
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


