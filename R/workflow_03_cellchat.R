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
      and red{{`...`}} would passed to `CellChat::createCellChat`."
    )
    if (group.by == "SingleR_cell") {
      if (x@step < 4L)
        stop("`step4(x)` has not been performed.")
    }
    object <- CellChat::createCellChat(
      object = object(x), group.by = group.by, assay = assay,
      ...
    )
    .job_cellchat(object = object)
  })

setMethod("step0", signature = c(x = "job_cellchat"),
  function(x){
    callNextMethod()
    step_message("Prepare your data with methods `asjob_cellchat`. ",
      "A processed 'Seurat' object is needed."
    )
  })

setMethod("step1", signature = c(x = "job_cellchat"),
  function(x, db = CellChat::CellChatDB.human,
    ppi = CellChat::PPI.human, cl = 4)
  {
    x <- callNextMethod(x)
    step_message("One step forward computation of all. By default,
      red{{`CellChat::CellChatDB.human`}}
      grey{{(Consider use `CellChat::subsetDB` to filter `db`)}}
      was used as reference database.
      The `future::plan('multisession', worker = cl)` were called.
      Then run: `CellChat::subsetData`;
      `CellChat::identifyOverExpressedGenes`;
      `CellChat::identifyOverExpressedInteractions`.
      By default, red{{`CellChat::PPI.human`}} were used by `CellChat::projectData`
      to smooth gene expression level (set NULL to cancel).
      Then run: `CellChat::computeCommunProb`;
      `CellChat::filterCommunication`; 
      `CellChat::computeCommunProbPathway`;
      `CellChat::aggregateNet`.
      "
    )
    p.showdb <- CellChat::showDatabaseCategory(db)
    p.showdb <- wrap(p.showdb, 8, 4)
    object(x)@DB <- db
    future::plan("multisession", workers = cl)
    object(x) <- CellChat::subsetData(object(x))
    object(x) <- CellChat::identifyOverExpressedGenes(object(x))
    object(x) <- CellChat::identifyOverExpressedInteractions(object(x))
    if (!is.null(ppi)) {
      object(x) <- CellChat::projectData(object(x), ppi)
    }
    object(x) <- CellChat::computeCommunProb(object(x))
    object(x) <- CellChat::filterCommunication(object(x), min.cells = 10)
    lp_net <- as_tibble(CellChat::subsetCommunication(object(x)))
    object(x) <- CellChat::computeCommunProbPathway(object(x))
    pathway_net <- as_tibble(CellChat::subsetCommunication(object(x), slot.name = "netP"))
    object(x) <- CellChat::aggregateNet(object(x))
    p.comms <- plot_communication.cellchat(object(x))
    x@plots[[ 1 ]] <- c(namel(p.showdb), p.comms)
    x@tables[[ 1 ]] <- namel(lp_net, pathway_net)
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
