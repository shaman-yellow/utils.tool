# ==========================================================================
# workflow of cellchatn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_cellchatn <- setClass("job_cellchatn", 
  contains = c("job"),
  prototype = prototype(
    info = c("https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html"),
    cite = "[@InferenceAndAJinS2021]",
    method = "R package `CellChat` used for cell communication analysis",
    tag = "scrna:cellchat",
    analysis = "CellChat 细胞通讯分析"
    ))

setGeneric("asjob_cellchatn",
  function(x, ...) standardGeneric("asjob_cellchatn"))

setMethod("asjob_cellchatn", signature = c(x = "list"),
  function(x, nms = names(x))
  {
    if (any(!vapply(x, is, logical(1), "job_cellchat"))) {
      message("All elements of list should be 'job_cellchat'.")
    }
    if (any(vapply(x, function(x) x@step < 2, logical(1)))) {
      stop('any(vapply(x, function(x) x@step < 2, logical(1))).')
    }
    meth <- meth(x[[1]])$step1
    p.lr_comm_bubbles <- lapply(x, 
      function(x) {
        x@plots$step2$lr_comm_bubble
      })
    t.lr_comm_bubble <- lapply(x, 
      function(x) {
        x@tables$step2$t.lr_comm_bubble
      })
    group.by <- vapply(x, function(x) x$group.by, character(1))
    if (!length(group.by <- unique(group.by))) {
      stop('!length(unique(group.by)).')
    }
    args_inters <- sapply(c("count", "weight"), simplify = FALSE,
      function(type) {
        weight.max <- e(CellChat::getMaxWeight(
            lapply(x, function(x) object(x)), attribute = c("idents", type)
            ))
        nets <- lapply(x, function(x) object(x)@net[[type]])
        list(nets = nets, wmax = weight.max, type = type)
      })
    object <- e(CellChat::mergeCellChat(lapply(x, function(x) object(x)), add.names = nms))
    x <- .job_cellchatn(object = object)
    x$.pre_meth <- meth
    x$each_lr_comm <- namel(p.lr_comm_bubbles, t.lr_comm_bubble)
    x$args_inters <- args_inters
    x$group.by <- group.by
    return(x)
  })


setMethod("step0", signature = c(x = "job_cellchatn"),
  function(x){
    step_message("Prepare your data with function `job_cellchatn`.")
  })

setMethod("step1", signature = c(x = "job_cellchatn"),
  function(x){
    step_message("Plot group comparison.")
    p.inters_counts <- funPlot(
      .plot_interactions_across_datasets, x$args_inters$count
    )
    p.inters_counts <- set_lab_legend(
      wrap(p.inters_counts, 14, 7),
      glue::glue("{x@sig} Comparison Number communication network"),
      glue::glue("数量通讯网络|||不同细胞的连线表示潜在的通讯关系，通讯强度数量越多，连线越粗。")
    )
    p.inters_weights <- funPlot(
      .plot_interactions_across_datasets, x$args_inters$weight
    )
    p.inters_weights <- set_lab_legend(
      wrap(p.inters_weights, 14, 7),
      glue::glue("{x@sig} Comparison Strength communication network"),
      glue::glue("强度通讯网络|||不同细胞的连线表示潜在的通讯关系，通讯强度数量越多，连线越粗。")
    )
    x <- methodAdd(x, x$.pre_meth)
    x <- methodAdd(x, "参照 <https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html> 对多数据集单细胞数据进行组间比较分析。")
    x <- plotsAdd(x, p.inters_counts, p.inters_weights)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cellchatn"),
  function(x, cell){
    step_message("Cell as source or target interactions with other cells.")
    allCells <- levels(object(x)@meta[[ x$group.by ]])
    if (length(cell) > 1) {
      stop('length(cell) > 1.')
    }
    if (!any(isThat <- grpl(allCells, cell))) {
      stop('!any(isThat <- grpl(allCells, cell)), can not found.')
    }
    keyCell <- allCells[ isThat ]
    otherCells <- allCells[ !isThat ]
    p.keyCellAsSource <- e(
      CellChat::netVisual_bubble(
        object(x), comparison = c(1, 2),
        sources.use = keyCell, targets.use = otherCells, 
        angle.x = 45
      )
    )
    data.keyCellAsSource <- .get_ggplot_content(p.keyCellAsSource)
    p.keyCellAsSource <- set_lab_legend(
      p.keyCellAsSource,
      glue::glue("{x@sig} key cell as source comparison"),
      glue::glue("关键细胞作为通讯发送方与其他细胞之间的通讯在组间的比较|||横坐标的不同颜色代表不同组别。")
    )
    p.keyCellAsTarget <- e(
      CellChat::netVisual_bubble(
        object(x), comparison = c(1, 2),
        sources.use = otherCells, targets.use = keyCell, 
        angle.x = 45
      )
    )
    x <- methodAdd(x, "利用 netVisual_bubble 函数绘制气泡图对各个细胞类型中配体受体介导的相互作用。")
    data.keyCellAsTarget <- .get_ggplot_content(p.keyCellAsTarget)
    p.keyCellAsTarget <- set_lab_legend(
      p.keyCellAsTarget,
      glue::glue("{x@sig} key cell as target comparison"),
      glue::glue("关键细胞作为通讯接收方与其他细胞之间的通讯在组间的比较|||横坐标的不同颜色代表不同组别。")
    )
    p.allCells_LP_comm_each_group <- x$each_lr_comm$p.lr_comm_bubbles
    p.allCells_LP_comm_each_group <- set_lab_legend(
      p.allCells_LP_comm_each_group,
      glue::glue("{x@sig} {names(p.allCells_LP_comm_each_group)} ligand receptor interactions bubble plot"),
      glue::glue("Group {names(p.allCells_LP_comm_each_group)}: 不同细胞类型之间的配体-受体对相互作用气泡图。|||纵坐标为配体-受体对，横坐标为细胞-细胞相互作用方向，颜色代表交互的可能性，颜色越红代表通讯可能越高，气泡的点大小代表显著性。")
    )
    x <- methodAdd(x, "以 `CellChat::netVisual_bubble` 比较关键细胞相关的组间差异配体与受体相互作用 (P value &lt; 0.05)。")
    lr_comm <- x$each_lr_comm$t.lr_comm_bubble
    s.com <- glue::glue(
      "{names(lr_comm)} 组包含 {vapply(lr_comm, nrow, integer(1))} 对唯一互作"
    )
    x$keyCell_data <- namel(data.keyCellAsSource, data.keyCellAsTarget)
    x <- snapAdd(x, "各组细胞间的配体、受体通讯统计{aref(p.allCells_LP_comm_each_group)}，{bind(s.com)} (P value &lt; 0.05)。")
    # x <- snapAdd(x, "其中，{keyCell} 作为互作发送方。")
    x <- plotsAdd(x, p.keyCellAsSource, p.keyCellAsTarget, p.allCells_LP_comm_each_group)
    return(x)
  })

.plot_interactions_across_datasets <- function(nets, 
  wmax, type, nms = names(nets))
{
  if (is.null(nms)) {
    stop('is.null(nms).')
  }
  par(mfrow = c(1, length(nets)), xpd = TRUE)
  type <- switch(type, count = "Number", weight = "Strength")
  lapply(seq_along(nets), 
    function(i) {
      CellChat::netVisual_circle(
        nets[[i]], weight.scale = TRUE,
        label.edge= FALSE, edge.weight.max = wmax[2],
        edge.width.max = 12, title.name = paste0(
          type, " of interactions - ", nms[i]
        )
      )
    })
}


