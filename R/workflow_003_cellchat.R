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
    info = c("https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html"),
    cite = "[@InferenceAndAJinS2021]",
    method = "R package `CellChat` used for cell communication analysis",
    tag = "scrna:cellchat",
    analysis = "CellChat 细胞通讯分析"
    ))

setGeneric("asjob_cellchat", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_cellchat"))

setMethod("asjob_cellchat", signature = c(x = "job_seurat"),
  function(x, group.by = x$group.by, assay = SeuratObject::DefaultAssay(object(x)), 
    sample = .5, ...)
  {
    step_message("Prarameter red{{`group.by`}}, red{{`assay`}}, 
      and red{{`...`}} would passed to `CellChat::createCellChat`.",
      show_end = NULL
    )
    if (x@step < 2) {
      stop("x@step < 2. At least, preprocessed assay data should be ready.")
    }
    if (!missing(sample) || ncol(object(x)) > 5e4) {
      message(glue::glue("Too many cells, sampling for analysis."))
      x <- asjob_seurat_sub(x, sample = sample)
    }
    if (is.null(group.by))
      stop("is.null(group.by) == TRUE")
    if (is.character(object(x)@meta.data[[ group.by ]])) {
      object(x)@meta.data[[ group.by ]] %<>% as.factor()
    } else {
      object(x)@meta.data[[ group.by ]] %<>% droplevels()
    }
    object <- e(CellChat::createCellChat(
        object = object(x), group.by = group.by, assay = assay,
        ...
        ))
    x <- snapAdd(
      .job_cellchat(object = object, params = list(group.by = group.by)), x
    )
    x <- methodAdd(x, "以 `CellChat` R 包 ({packageVersion('CellChat')}) {cite_show('InferenceAndAJinS2021')} 对单细胞数据进行细胞通讯分析。以 `CellChat::createCellChat` 将 `Seurat` 对象的 {assay} Assay 转化为 CellChat 对象。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_cellchat"),
  function(x){
    step_message("Prepare your data with methods `asjob_cellchat`. ",
      "A processed 'Seurat' object is needed."
    )
  })

setMethod("step1", signature = c(x = "job_cellchat"),
  function(x, workers = 4, python = NULL, py_config = FALSE, debug = FALSE, 
    org = c("human", "mouse"), ...)
  {
    step_message("One step forward computation of most. ")
    org <- match.arg(org)
    if (is.remote(x)) {
      ## need python umap-learn, the best way is to use reticulate python.
      ## reticulate::py_install(packages = 'umap-learn')
      x <- run_job_remote(x, wait = 3, ..., debug = debug,
        {
          options(future.globals.maxSize = 5e10)
          x <- step1(
            x, workers = "{workers}", org = "{org}",
            debug = "{debug}", python = "{python}", py_config = "{py_config}"
          )
        }
      )
      x@plots$step1$p.commHpAll@data <- e(
        CellChat::netVisual_heatmap(object(x), color.heatmap = "Reds", signaling = NULL)
      )
      return(x)
    }
    if (org == "human") {
      db <- CellChat::CellChatDB.human
      ppi <- CellChat::PPI.human
    } else if (org == "mouse") {
      db <- CellChat::CellChatDB.mouse
      ppi <- CellChat::PPI.mouse
    }
    if (!is.null(python)) {
      e(base::Sys.setenv(RETICULATE_PYTHON = python))
      e(reticulate::py_config())
    } else if (py_config) {
      e(reticulate::py_config())
    }
    if (!debug) {
      p.showdb <- e(CellChat::showDatabaseCategory(db))
      p.showdb <- wrap(p.showdb, 8, 4)
      object(x)@DB <- db
      future::plan("multisession", workers = workers)
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
      types <- c("count", "weight", "individuals")
      p.comms <- .set_lab(p.comms, sig(x), paste("overall communication", types))
      p.comms <- setLegend(p.comms, glue::glue("为细胞通讯 {types} 统计。"))
      ## Signaling role of cell groups
      object(x) <- e(CellChat::netAnalysis_computeCentrality(object(x), slot.name = "netP"))
    }
    ## Clustering
    object <- object(x)
    res <- try(
      e({
        for (i in c("functional", "structural")) {
          object(x) <- CellChat::computeNetSimilarity(object(x), type = i)
          object(x) <- CellChat::netEmbedding(object(x), type = i)
          object(x) <- CellChat::netClustering(object(x), type = i, do.parallel = FALSE)
        }
      })
    )
    if (inherits(res, "try-error")) {
      message("Due to error, escape from clustering; But the object was returned.")
      object(x) <- object
    }
    p.commHpAll <- CellChat::netVisual_heatmap(object(x), color.heatmap = "Reds", signaling = NULL)
    p.commHpAll <- .set_lab(wrap(p.commHpAll), sig(x), "All Cell communication heatmap")
    p.commHpAll <- setLegend(p.commHpAll, "为所有细胞的通讯热图。")
    #  Interaction net count plot and interaction weight
    #  plot. The thicker the line represented, the more
    #  the number of interactions, and the stronger the interaction
    #  weights/strength between the two cell types
    if (!debug) {
      x@plots[[ 1 ]] <- c(namel(p.showdb, p.commHpAll), p.comms)
      x@tables[[ 1 ]] <- namel(lp_net, pathway_net)
    }
    x <- methodAdd(x, "参照 {x@info} 分析 scRNA-seq 数据。")
    x <- snapAdd(x, "对 {showStrings(object(x)@meta[[x$group.by]], trunc = FALSE)} 细胞通讯分析。")
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
    if (FALSE) {
      cell_comm_heatmap <- e(sapply(pathways, simplify = FALSE,
        function(name) {
          main <- try(CellChat::netVisual_heatmap(
            object(x), color.heatmap = "Reds", signaling = name), TRUE)
          if (inherits(main, "try-error")) {
            return(list(main = NULL, contri = NULL))
          }
          if (!is.null(name)) {
            contri <- CellChat::extractEnrichedLR(object(x),
              signaling = name, geneLR.return = FALSE)
            return(namel(main, contri))
          } else {
            return(namel(wrap(main)))
          }
        }))
      cell_comm_heatmap$ALL$main <- .set_lab(cell_comm_heatmap$ALL$main, sig(x), "Cell communication heatmap")
    } else {
      cell_comm_heatmap <- NULL
    }
    lr_comm_bubble <- e(CellChat::netVisual_bubble(object(x), remove.isolate = FALSE))
    lr_comm_bubble <- set_lab_legend(
      lr_comm_bubble,
      glue::glue("{x@sig} communication probability and significant"),
      glue::glue("为细胞通讯概率以及显著性。")
    )
    gene_expr_violin <- e(sapply(pathways, simplify = FALSE,
      function(name) {
        CellChat::plotGeneExpression(object(x),
            signaling = name, group.by = NULL) +
          theme(legend.position = "none")
      }))
    if (FALSE) {
      # this plot saved as RDS too large.
      role_comps_heatmap <- e(sapply(pathways, simplify = FALSE,
          function(name) {
            res <- try({CellChat::netAnalysis_signalingRole_network(object(x), signaling = name,
              width = 8, height = 2.5, font.size = 8, cluster.rows = TRUE)
              recordPlot()}, TRUE)
            if (inherits(res, "try-error"))
              return(NULL)
            else res
          }))
    } else {
      role_comps_heatmap <- NULL
    }
    role_weight_scatter <- e(sapply(pathways, simplify = FALSE,
        function(name) {
          CellChat::netAnalysis_signalingRole_scatter(object(x),
            signaling = name)
        }))
    res <- try(lr_role_heatmap <- e(sapply(c("outgoing", "incoming", "all"), simplify = FALSE,
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
    lr_role_heatmap <- setLegend(
      lr_role_heatmap, glue::glue("为细胞间的 {names(lr_role_heatmap)} 类型通路信号强度 ")
    )
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
        object(x), pattern = "outgoing", k = k.out, heatmap.show = FALSE
        ))
    object(x) <- CellChat::identifyCommunicationPatterns(
      object(x), pattern = "incoming", k = k.in, heatmap.show = FALSE
    )
    require(ggalluvial)
    lst.p <- e(sapply(c("outgoing", "incoming"), simplify = FALSE,
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
    weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  p.aggre_weight <- CellChat::netVisual_circle(x@net$weight, vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  mat <- x@net$weight
  scale <- cal_panelScale(length(levels(x@idents)))
  par(mfrow = scale, xpd = TRUE, omi = rep(0, 4), mar = rep(1.5, 4))
  for (i in seq_len(nrow(mat))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    CellChat::netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE,
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
  function(x, pattern.x, pattern.y = pattern.x, get = c("lps", "pathways", "intersect")){
    message("`Pattern` used for match cells.")
    get <- match.arg(get)
    if (x@step < 1) {
      stop("x@step != 1")
    }
    if (pattern.x == pattern.y) {
      `%with%` <- `|`
    } else {
      `%with%` <- `&`
    }
    pathways <- filter(x@tables$step1$pathway_net,
      grepl(pattern.x, source) %with% grepl(pattern.y, target),
      pval < .05
    )
    pathways <- .set_lab(pathways, sig(x), "interaction of pathways")
    if (get == "pathways") {
      return(pathways)
    }
    lps <- filter(x@tables$step1$lp_net,
      grepl(pattern.x, source) %with% grepl(pattern.y, target),
      pval < .05
    )
    lps <- .set_lab(lps, sig(x), "interaction of ligand and receptor")
    if (get == "lps") {
      return(lps)
    }
    ints <- intersect(unique(pathways$pathway_name), unique(lps$pathway_name))
    if (get == "intersect") {
      return(ints)
    }
  })

unique.lps <- function(x, use = c("ligand", "receptor")) {
  if (is.data.frame(x)) {
    x <- unlist(lapply(use, function(n) x[[ n ]]))
  }
  unique(unlist(stringr::str_extract_all(x, "[^_]+")))
}

setMethod("set_remote", signature = c(x = "job_cellchat"),
  function(x, wd = glue::glue("~/cellchat_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

