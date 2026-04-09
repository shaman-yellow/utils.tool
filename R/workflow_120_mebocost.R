# ==========================================================================
# workflow of mebocost
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mebocost <- setClass("job_mebocost", 
  contains = c("job"),
  prototype = prototype(
    pg = "mebocost",
    info = c("https://github.com/kaifuchenlab/MEBOCOST"),
    cite = "",
    method = "",
    tag = "mebocost",
    analysis = "MEBOCOST 细胞代谢通讯分析"
    ))

setGeneric("asjob_mebocost",
  function(x, ...) standardGeneric("asjob_mebocost"))

setMethod("asjob_mebocost", signature = c(x = "job_seurat"),
  function(x, dir_cache = create_job_cache_dir(x, "mebocost"), 
    conda_env = pg("mebocostEnv"))
  {
    message("Convert data.")
    cli::cli_alert_info("scCustomize::Convert_Assay")
    object <- scCustomize::Convert_Assay(
      object(x), "RNA", convert_to = "V3"
    )
    hash <- digest::digest(
      list(cells = colnames(object(x)), genes = rownames(object(x))), 
      "xxhash64", serializeVersion = 3
    )
    file_anndata <- file.path(dir_cache, glue::glue("anndata_{hash}.h5ad"))
    if (!file.exists(file_anndata)) {
      activate_env(conda_env, pg("conda"))
      cli::cli_alert_info("sceasy::convertFormat")
      sceasy::convertFormat(
        object, from = "seurat", to = "anndata",
        outFile = file_anndata
      )
    } else {
      message(glue::glue('file.exists(file_anndata): {file_anndata}'))
    }
    metadata <- as_tibble(object(x)@meta.data, idcol = "cell")
    levels <- .guess_levels_from_job_seurat(x)
    group.by <- x$group.by
    x <- .job_mebocost()
    x$dir_cache <- dir_cache
    x$file_anndata <- file_anndata
    x$metadata <- metadata
    x$levels <- levels
    x$group.by <- group.by
    x$conda_env <- conda_env
    x$hash <- hash
    return(x)
  })

setMethod("step0", signature = c(x = "job_mebocost"),
  function(x){
    step_message("Prepare your data with function `job_mebocost`.")
  })

setMethod("activate", signature = c(x = "job_mebocost"),
  function(x){
    if (is.null(x$conda_env)) {
      stop('is.null(x$conda_env).')
    }
    activate_env(x$conda_env, pg("conda"))
    x$scanpy <- e(reticulate::import("scanpy"))
    x$mebocost <- e(reticulate::import("mebocost")$mebocost)
    return(x)
  })

setMethod("step1", signature = c(x = "job_mebocost"),
  function(x, workers = 10, species = "human", path_mebocost = pg("path_mebocost"),
    try_compass = FALSE)
  {
    step_message("Create mebocost object.")
    # scrub <- import("scrublet")
    x <- activate(x)
    anndata <- x$scanpy$read_h5ad(x$file_anndata)
    x$file_config <- file.path(path_mebocost, "mebocost.conf")
    if (!is.null(x$metadata$group)) {
      x$compare.by <- "group"
    } else {
      x$compare.by <- NULL
    }
    cli::cli_alert_info("Run: mebocost.create_obj")
    object(x) <- x$mebocost$create_obj(
      adata = anndata,
      group_col = x$group.by,
      condition_col = x$compare.by,
      met_est = "mebocost",
      config_path = x$file_config,
      exp_mat = NULL, cell_ann = NULL, species = species,
      met_pred = NULL, met_enzyme = NULL, met_sensor = NULL,
      met_ann = NULL, scFEA_ann = NULL, compass_met_ann = NULL,
      compass_rxn_ann = NULL,
      cutoff_exp = 'auto',
      cutoff_met = 'auto',
      cutoff_prop = 0.15,
      sensor_type = 'All',
      thread = workers
    )
    if (try_compass) {
      x$file_avgExp <- file.path(x$dir_cache, glue::glue("avgExp_{x$hash}.tsv"))
      if (file.exists(x$file_avgExp)) {
        message(glue::glue("file.exists: {x$file_avgExp}"))
      } else {
        pd <- reticulate::import("pandas")
        np <- reticulate::import("numpy")
        avg_exp <- x$scanpy$get$aggregate(
          anndata, by = list("group", x$group.by), func = "mean"
        )
        cols <- avg_exp$var_names$to_list()
        rows <- paste0(avg_exp$obs[[ "group" ]], " ~ ", avg_exp$obs[[ x$group.by ]])
        avg_exp <- as.data.frame(avg_exp$layers[[ "mean" ]])
        rownames(avg_exp) <- rows
        colnames(avg_exp) <- cols
        avg_exp <- t(avg_exp)
        avg_exp <- expm1(avg_exp)
        data.table::fwrite(
          avg_exp, x$file_avgExp, sep = "\t", row.names = TRUE
        )
      }
    }
    x <- methodAdd(x, "**MEBOCOST** 是一种基于单细胞转录组数据推断代谢物介导细胞间通讯的计算方法，其主要目的是系统性解析不同细胞类型之间通过代谢物–受体轴所形成的潜在相互作用网络。该方法结合代谢酶表达信息与受体表达谱，推断细胞产生特定代谢物的能力及其被其他细胞感知的可能性，从而构建代谢通讯关系，并在此基础上识别具有显著性的通讯通路及关键调控分子。\n\n")
    x <- methodAdd(x, "将单细胞数据集以 Python 工具 MEBOCOST (1.2.0) (<https://github.com/kaifuchenlab/MEBOCOST>) 分析细胞间代谢通讯。")
    x <- methodAdd(x, "借助 MEBOCOST 内置的代谢物-聚合酶、代谢物-传感器先验知识库计算各细胞群体中聚合酶、传感器基因的平均表达量，作为该群体对该代谢物发送、接收潜力的评估指标。\n\n")
    x <- methodAdd(x, "对于每个代谢物，将⟦mark$blue('聚合酶表达得分高于默认阈值')⟧（通常为所有细胞 25%）的细胞类型定义为潜在 Sender 亚群；将对应⟦mark$blue('传感器基因平均表达高于默认阈值')⟧的细胞类型定义为潜在 Receiver 亚群。同时要求⟦mark$blue('代谢物相关酶与传感器基因在相应细胞群体中的表达细胞比例不低于15%')⟧，以确保群体代表性。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_mebocost"),
  function(x, species = "homo_sapiens", workers = 1L)
  {
    step_message("Run compass")
    if (x$.args$step1$try_compass && !is.null(x$file_avgExp)) {
      x$dir_compass <- file.path(x$dir_cache, "compass_output")
      x$dir_compass_tmp <- file.path(x$dir_cache, "compass_tmp")
      input <- glue::glue("--data {x$file_avgExp} ")
      output <- glue::glue("--output-dir {x$dir_compass} --temp-dir {x$dir_compass_tmp} ")
      setting <- glue::glue("--num-thread {workers} --species {species} --calc-metabolites --lambda 0")
      system(glue::glue("{pg('compass')} {input} {setting} {output}"))
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_mebocost"),
  function(x, min_cell_number = 10L, cut.p = .05, use.p = "permutation_test_fdr", rerun = FALSE)
  {
    step_message("Infer communication.")
    x <- activate(x)
    x$use.p <- use.p
    fun_infer <- function(min_cell_number, cut.p, metadata) {
      object(x)$infer_commu(
        n_shuffle = 1000L,
        seed = as.integer(x$seed),
        Return = FALSE,
        thread = as.integer(x$.args$step1$workers),
        save_permuation = TRUE,
        min_cell_number = as.integer(min_cell_number),
        pval_method = x$use.p,
        pval_cutoff = cut.p
      )
      object(x)
    }
    object(x) <- expect_local_data(
      x$dir_cache, "mebocost", fun_infer, list(min_cell_number, cut.p, x$metadata),
      fun_read = x$mebocost$load_obj, fun_save = x$mebocost$save_obj, 
      ext = "pk", rerun = rerun
    )
    t.commu_res <- as_tibble(object(x)$commu_res)
    t.commu_res <- dplyr::filter(t.commu_res, !!rlang::sym(x$use.p) < !!cut.p)
    t.commu_res <- .mutate_get_chain_in_mebocost_table(t.commu_res)
    t.commu_res <- set_lab_legend(
      t.commu_res,
      glue::glue("{x@sig} cell metabolic communication results"),
      glue::glue("酶和传感器共表达检测的显著代谢物介导的细胞间通讯")
    )
    ps.heatmaps <- sapply(x$levels, simplify = FALSE,
      function(group) {
        p <- vis(x, "commu_dotmap", group = group, cut.p = cut.p)
        p <- set_lab_legend(p,
          glue::glue("{x@sig} group {group} significant communication heatmap"),
          glue::glue("Group: {group} 细胞间代谢通讯气泡图|||展示显著的“代谢物–感受器”对在不同细胞对之间的通讯关系，其中纵轴为代谢物及其对应感受器的组合，横轴为具体的发送细胞与接收细胞配对，每个气泡代表一条显著的代谢通讯事件；气泡大小表示通讯强度，数值越大代表该通讯关系越强，气泡颜色表示统计显著性水平。")
        )
      })
    snap_commu <- .stat_table_by_pvalue(
      t.commu_res, n = 5, split = "Condition", use.p = x$use.p, 
      colName = "Chain", target = "细胞代谢通讯", by = "组中检测到"
    )
    x <- snapAdd(x, "通过 MEBOCOST `infer_commu` 在单细胞数据集中一共检测到{aref(ps.heatmaps)} {snap_commu}")
    p.eventnum_bar <- vis(x, "eventnum_bar")
    p.eventnum_bar <- set_lab_legend(
      p.eventnum_bar,
      glue::glue("{x@sig} bar plot of communication events"),
      glue::glue("通讯事件柱状图|||柱状图展示发送方与接收方的通讯数量，横轴为各细胞类型，纵轴为通讯事件数。")
    )
    x$cut.p <- cut.p
    x <- methodAdd(x, "对每一对 Sender-Receiver 细胞类型及每一组代谢物-传感器对，计算酶-传感器共表达得分（Sender细胞中代谢物聚合酶表达均值 × Receiver 细胞中传感器表达均值）作为原始通讯强度。通过 1000 次细胞标签置换检验构建零分布，计算经验 p 值，并经 Benjamini‑Hochberg 法进行 FDR 校正，⟦mark$blue('以 FDR &lt; {cut.p} 为阈值筛选出显著的代谢物-传感器结合概率')⟧。")
    x <- tablesAdd(x, t.commu_res)
    x <- plotsAdd(x, ps.heatmaps, p.eventnum_bar)
    return(x)
  })

setMethod("step4", signature = c(x = "job_mebocost"),
  function(x, flux_pass = TRUE, sig_mccc_only = TRUE, 
    cut.p = x$cut.p, cut.fc = .5, rerun = FALSE)
  {
    step_message("Differential analysis.")
    x <- activate(x)
    compare <- paste0(x$levels[1], "_vs_", x$levels[2])
    cli::cli_alert_info("CommDiff")
    if (flux_pass && is.null(x$is_compass_run)) {
      flux_pass <- FALSE
    }
    fun_diff <- function(...) {
      object(x)$CommDiff(
        comps = as.list(compare),
        sig_mccc_only = sig_mccc_only,
        flux_pass = flux_pass,
        thread = as.integer(x$.args$step1$workers)
      )
      object(x)
    }
    object(x) <- expect_local_data(
      x$dir_cache, "diff", fun_diff, list(flux_pass, sig_mccc_only),
      fun_read = x$mebocost$load_obj, fun_save = x$mebocost$save_obj, 
      ext = "pk", rerun = rerun
    )
    x <- methodAdd(x, "对 {bind(x$levels)} 的细胞代谢通讯进行组间比较。")
    ts.diff_commu <- sapply(compare, simplify = FALSE,
      function(com) {
        data <- tibble::as_tibble(object(x)$diffcomm_res[[com]])
        data <- dplyr::filter(data, abs(Log2FC) > !!cut.fc, !!rlang::sym(x$use.p) < !!cut.p)
        data <- .mutate_get_chain_in_mebocost_table(data)
        data <- set_lab_legend(
          data,
          glue::glue("{x@sig} {com} differential communication"),
          glue::glue("{com} 组间细胞代谢通讯差异分析表格。")
        )
      })
    x <- tablesAdd(x, ts.diff_commu)
    p.diff_flow <- vis(
      x, "diff_flow", compare = compare, cut.p = cut.p, cut.fc = cut.fc
    )
    p.diff_flow <- set_lab_legend(
      p.diff_flow,
      glue::glue("{x@sig} Cellular differential metabolic communication network"),
      glue::glue("细胞差异代谢通讯网络|||{.mebocost_network_note}")
    )
    maxShow <- 5L
    snap_diff <- vapply(names(ts.diff_commu), FUN.VALUE = character(1),
      function(name) {
        data <- ts.diff_commu[[name]]
        data <- dplyr::arrange(data, dplyr::desc(abs(Log2FC)))
        data <- head(data, n = maxShow)
        cp <- .setup_compare_pvalue_with_table(
          data, "Chain", x$use.p, "Log2FC", levels = x$levels
        )
        snap <- .stat_compare_by_pvalue(cp, x$levels, "", mode = "communication")
        glue::glue("对显著性结果按 Log2FC 降序排序，如图{aref(p.diff_flow)}，排名前{nrow(data)}的细胞代谢通讯中，{snap}")
      })
    snap_diff <- bind(snap_diff, co = "\n\n")
    x <- snapAdd(x, "⟦mark$red('对 {bind(x$levels)} 通讯差异分析，一共检测到 {nrow(ts.diff_commu[[1]])} 个显著差异的细胞代谢通讯')⟧。{snap_diff}")
    x <- plotsAdd(x, p.diff_flow)
    return(x)
  })

setMethod("step5", signature = c(x = "job_mebocost"),
  function(x, use.score = c("scale", "raw"), key = "Metabolite_Name",
    group_by = c("Metabolite_Name", "Receiver"),
    axis = c("Sender", "Metabolite_Name", "Sensor", "Receiver"))
  {
    step_message("PageRank.")
    use.score <- match.arg(use.score)
    if (use.score == "scale") {
      # x$levels[1] is the 'disease' group
      use.score <- glue::glue("Scaled_Commu_Score_{x$levels[1]}")
    } else {
      use.score <- glue::glue("Commu_Score_{x$levels[1]}")
    }
    raw <- x@tables$step4$ts.diff_commu[[1]]
    raw <- dplyr::filter(raw, !!rlang::sym(x$use.p) < .05)
    data <- dplyr::select(raw, dplyr::all_of(axis), Log2FC, !!rlang::sym(use.score))
    fshow <- function(x) strx(x, "[^_]+")
    name_axis <- bind(fshow(group_by), co = "——")
    x <- methodAdd(
      x, "以多维度综合策略锁定关键 {name_axis} 通讯轴。"
    )
    snap_score <- .get_mebocost_pagerank_note(x$levels[1], group_by, fshow(key))
    x <- methodAdd(x, "该综合策略通过计算综合得分实现：\n\n{snap_score}\n\n")
    # calculate score
    data <- dplyr::mutate(
      data, weight = !!rlang::sym(use.score) * abs(Log2FC)
    )
    nodes <- dplyr::select(data, dplyr::all_of(axis))
    nodes <- tidyr::pivot_longer(
      nodes, dplyr::everything(), names_to = "type", values_to = "name"
    )
    nodes <- dplyr::distinct(nodes)
    edges <- lapply(seq_len(length(axis) - 1),
      function(n) {
        dplyr::select(data, from = !!n, to = !!(n + 1), weight)
      })
    edges <- dplyr::bind_rows(edges)
    data_score <- .get_page_rank_score(edges, "weight")
    nodes <- map(
      nodes, "name", data_score, "name", "score", col = "score"
    )
    # nodes <- dplyr::arrange(nodes, dplyr::desc(score))
    data <- dplyr::group_by(data, !!!rlang::syms(group_by))
    data <- dplyr::summarize(data, weight_sum = sum(weight))
    # get pagerank score of `key`
    data <- map(data, key, nodes, "name", "score", col = "score")
    data <- dplyr::mutate(data, overall_score = score * weight_sum)
    data <- dplyr::arrange(data, dplyr::desc(overall_score))
    if (length(group_by) == 2) {
      p.score <- ggplot(data, aes(x = !!rlang::sym(group_by[1]), y = !!rlang::sym(group_by[2]))) +
        geom_point(aes(size = overall_score)) +
        labs(x = fshow(group_by[1]), y = fshow(group_by[2]), size = "Overall score") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      p.score <- set_lab_legend(
        wrap_scale_heatmap(p.score, data[[group_by[1]]], data[[group_by[2]]]),
        glue::glue("{x@sig} Dotplot for overall score"),
        glue::glue("综合得分气泡图|||根据 MEBOCOST 网络计算得到的综合得分。横纵坐标对应于 {bind(fshow(group_by))}。")
      )
      x <- plotsAdd(x, p.score)
    } else {
      p.score <- NULL
    }
    t.overallScore <- set_lab_legend(data,
      glue::glue("{x@sig} Overall score for mebocost network"),
      glue::glue("根据 MEBOCOST 网络计算得到的综合得分。")
    )
    x <- tablesAdd(x, t.overallScore)
    snap <- bind(c(data[[1]][1], data[[2]][1]), co = " -> ")
    for (i in group_by) {
      x[[ glue::glue(".feature_{tolower(fshow(i))}") ]] <- as_feature(
        data[[i]][1], "关键通讯轴", nature = fshow(i)
      )
    }
    x <- snapAdd(x, "根据多维度综合策略，如图{aref(p.score)}，得分最高的 {name_axis} 通讯轴为 {snap}。")
    return(x)
  })

.get_mebocost_pagerank_note <- function(model, group_by, key)
{
  dic <- c(Metabolite_Name = "m", Sender = "s", Sensor = "k", Receiver = "r")
  group_by <- dplyr::recode(group_by, !!!dic)
  others <- dic[ !dic %in% group_by ]
  formula <- glue_ex("\n\n$$\nOverallScore(⟦bind(group_by)⟧) = \\left( \\sum_{⟦bind(others)⟧} CommuScore_{⟦model⟧}(s,m,k,r) \\times |Log2FC(s,m,k,r)| \\right) \\times PR(⟦key⟧)\n$$\n\n")
  glue_ex(
    "⟦formula⟧\n其中，$s$ 表示 Sender（信号发送细胞类型），$m$ 表示 Metabolite（代谢物），$k$ 表示 Sensor（代谢物对应的受体或转运分子），$r$ 表示 Receiver（信号接收细胞类型）；$CommuScore_{⟦model⟧}(s,m,k,r)$ 表示在 ⟦model⟧ 组中针对完整通讯链（Sender–Metabolite–Sensor–Receiver）计算得到的通讯强度评分；$Log2FC(s,m,k,r)$ 表示 ⟦model⟧ 组相对于对照组的通讯强度对数倍数变化，其绝对值用于刻画差异幅度；$\\sum_{⟦bind(others)⟧}$ 表示对对应通讯节点进行求和，得到局部通讯强度；$PR(⟦key⟧)$ 表示在以差异加权权重为边权构建的有向网络（Sender→Metabolite→Sensor→Receiver）中计算得到的$⟦key⟧$节点 PageRank 值 (`igraph::page_rank`)，用于表征其在全局网络中的拓扑重要性。\n\n"
  )
}

glue_ex <- function(string, envir = parent.frame(1)) {
  glue::glue(string, .open = "⟦", .close = "⟧", .envir = envir)
}

.get_page_rank_score <- function(edges, weight, directed = TRUE) {
  graph <- igraph::graph_from_data_frame(edges, directed = directed)
  score <- igraph::page_rank(
    graph, weights = igraph::edge_attr(graph, weight)
    )$vector
  data <- as_df.lst(score, "name", "score")
  dplyr::arrange(data, dplyr::desc(score))
}

setMethod("vis", signature = c(x = "job_mebocost"),
  function(x, mode = c("eventnum_bar", "diff_flow", "commu_dotmap"), ...){
    mode <- match.arg(mode)
    if (mode == "eventnum_bar") {
      plot <- .plot_mebocost_eventnum_bar(x, ...)
    } else if (mode == "diff_flow") {
      plot <- .plot_mebocost_diff_flow(x, ...)
    } else if (mode == "commu_dotmap") {
      plot <- .plot_mebocost_commu_dotmap(x, ...)
    }
    file <- tempfile(mode, fileext = ".pdf")
    plot$savefig(file, bbox_inches = "tight")
    as_data_binary(.file_fig(file))
  })

.mutate_get_chain_in_mebocost_table <- function(data) {
  fun_fix <- function(x) s(x, "^[^~]+~ ", "")
  dplyr::mutate(
    data, Chain = paste0(
      fun_fix(Sender), " -> ", Metabolite_Name, " -> ", Sensor, " -> ", fun_fix(Receiver)
    )
  )
}

.plot_mebocost_eventnum_bar <- function(x, ...) {
  celltypes <- unique(x$metadata[[x$group.by]])
  groups <- unique(x$metadata$group)
  orders <- unlist(lapply(celltypes, function(x) paste0(groups, " ~ ", x)))
  object(x)$eventnum_bar(
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    ## uncomment and set to focus on one condition
    # conditions  =  ['Primary'],
    xorder = as.list(orders),
    and_or = "and",
    pval_method = x$use.p,
    pval_cutoff = 0.05,
    comm_score_col = "Norm_Commu_Score",
    comm_score_cutoff = 0,
    cutoff_prop = 0.25,
    figsize = c(1 + length(celltypes) * .5, 5),
    save = NULL,
    show_plot = FALSE,
    show_num = TRUE,
    include = list("sender-receiver"),
    group_by_cell = TRUE,
    colorcmap = "tab20",
    return_fig = TRUE
  )
}

.plot_mebocost_commu_dotmap <- function(x, group, cut.p = .05, flux_pass = TRUE)
{
  if (flux_pass && is.null(x$is_compass_run)) {
    flux_pass <- FALSE
  }
  object(x)$commu_dotmap(
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    conditions = list(group),
    and_or = 'and',
    flux_pass = flux_pass,
    pval_method = x$use.p,
    pval_cutoff = cut.p, 
    cmap = 'Reds',
    cellpair_order = c(),
    met_sensor_order = c(),
    show_plot = FALSE,
    comm_score_col = 'Commu_Score',
    comm_score_range = NULL,
    comm_score_cutoff = NULL,
    cutoff_prop = NULL,
    return_fig = TRUE
  )
}

.plot_mebocost_diff_flow <- function(x, compare, cut.p = .05, cut.fc = .5) {
  object(x)$DiffFlowPlot(
    comp_cond = compare, 
    pval_method = x$use.p,
    pval_cutoff = cut.p,
    Log2FC_threshold = cut.fc,
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    remove_unrelevant = TRUE,
    and_or = 'and',
    node_label_size = 8,
    node_alpha = .8,
    figsize = 'auto',
    node_cmap = 'Set1',
    line_color_col = 'Log2FC',
    line_cmap = 'bwr',
    line_cmap_vmin = -2,
    line_cmap_vmax = 2,
    line_cmap_center = 0,
    linewidth = 1.5,
    node_size_norm = c(10, 150),
    node_value_range = NULL,
    save = NULL, 
    save_plot = FALSE, 
    show_plot = FALSE,
    text_outline = FALSE,
    return_fig = TRUE
  )
}

.mebocost_network_note <- "基于 MEBOCOST 的细胞间差异代谢通讯网络可视化，按照“发送细胞（Sender）–代谢物（Metabolite）–感受器（Sensor）–接收细胞（Receiver）”四层结构展示完整的通讯路径，其中左侧为分泌代谢物的细胞类型，中间依次为参与通讯的小分子及其对应的受体或转运蛋白，右侧为表达感受器并接收信号的细胞类型；连线表示一条代谢通讯关系，颜色根据组间差异分析的 log2FC 显示变化方向与幅度（红色表示上调，蓝色表示下调，颜色越深代表变化越显著），节点大小表示该节点参与的连接数量（即连接度），用于反映其在网络中的参与程度；该网络仅包含经过差异分析筛选后的显著代谢通讯关系，用于整体呈现不同细胞类型之间通过代谢物介导的通讯模式。"

setMethod("set_remote", signature = c(x = "job_mebocost"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
