# ==========================================================================
# workflow of scissor  (AI) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scissor <- setClass("job_scissor",
  contains = "job",
  prototype = prototype(
    pg = "scissor",
    info = c(
      "https://github.com/sunduanchen/Scissor",
      "https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html"
    ),
    cite = "",
    method = "",
    tag = "scissor",
    analysis = "Scissor 表型关联单细胞亚群鉴定"
  )
)

setGeneric("do_scissor",
  function(x, ref, ...) standardGeneric("do_scissor"))

setMethod("do_scissor", signature = c(x = "job_seurat", ref = "job_deseq2"),
  function(x, ref, group = "group", family = "binomial", label = group,
    assay = SeuratObject::DefaultAssay(object(x)))
  {
    if (x@step < 3L) {
      stop('x@step < 3L.')
    }
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    levels <- NULL
    if (family == "binomial") {
      levels <- rev(.guess_compare_deseq2(ref))
      meta_bulk <- as.integer(factor(object(ref)@colData[[group]], levels = levels)) - 1L
      methodAdd_onExit(
        "x", "选择 Scissor 回归类型为 {family}，即构建二分类变量回归模型，计算每个细胞与表型的回归系数。将负回归系数的细胞描述为 “Scissor-”，在本研究中，与{label}为 {levels[1]} 的组高度相关；将正回归系数的细胞描述为 “Scissor+”，与{label}为 {levels[2]} 的组高度相关。回归系数为零的细胞为背景细胞。"
      )
    } else {
      stop('family != "binomial", not yet ready.')
    }
    expr_bulk <- ref$vst@assays@data[[1L]]
    object <- object(x)
    assayObj <- object[[ assay ]]
    if (is(assayObj, "Assay5") || is(assayObj, "SCTAssay")) {
      message(glue::glue("Detected 'Assay5' or 'SCTAssay', convert to 'Assay'."))
      assayObj <- as(assayObj, "Assay")
      object[[ "RNA" ]] <- assayObj
    }
    if (assay != "RNA") {
      # Scissor only extracts "RNA" assay internally!!!
      nameSnn <- glue::glue("{assay}_snn")
      object@graphs[[ "RNA_snn" ]] <- object@graphs[[ nameSnn ]]
      SeuratObject::DefaultAssay(object) <- "RNA"
      object[[ assay ]] <- NULL
    }
    pr <- params(x)
    project_x <- x$project
    if (is.null(project_x)) {
      project_x <- .guess_geo_project(x)
    }
    project_ref <- .guess_geo_project(ref)
    x <- .job_scissor(object = list(bulk_dataset = expr_bulk, sc_dataset = object, phenotype = meta_bulk))
    x <- methodAdd(x, "Scissor 是一种用于整合单细胞转录组数据与群体水平表型信息（如临床性状或 bulk 转录组数据）的分析方法，其主要目的是在单细胞分辨率下识别与特定表型显著相关的细胞亚群。该方法通过构建单细胞与 bulk 样本之间的表达相关性，并结合回归模型（如惩罚回归）筛选出与目标表型正相关或负相关的细胞群体（分别称为 Scissor+ 和 Scissor− 细胞）。在此基础上，可进一步对相关细胞群进行差异分析与功能注释。\n\n")
    x <- methodAdd(x, "选择 {project_x} 作为单细胞表达矩阵来源，以 {project_ref} 为 bulk 表达矩阵输入，且以 {project_ref} 的表型数据{label}作为目标表型，使用 R 包 `Scissor` ⟦pkgInfo('Scissor')⟧ 实现表型关联单细胞亚群分析。")
    x@params <- append(x@params, pr)
    x$family <- family
    x$tag <- levels
    return(x)
  }
)

setMethod("step1", signature = "job_scissor",
  function(x, mode = c("normal", "small", "middle"), workers = 4L,
    file_save = file.path(create_job_cache_dir(x), "scissor_inputs.qs"),
    rerun = FALSE, ...)
  {
    step_message("Run Scissor...")
    mode <- match.arg(mode)
    if (mode == "small") {
      alphas <- c(5e-3, 5e-4, 1e-4, 5e-5)
    } else if (mode == "middle") {
      alphas <- c(3e-3, 1e-3, 7e-4, 9e-4)
    } else if (mode == "normal") {
      alphas <- c(.005, .05, .1, .2)
    }
    res <- .optimize_scissor_with_speed_up_and_score_alpha(
      x, alphas, x$family, rerun = rerun, workers = workers, ...
    )
    x$stat_alphas_rough <- res$stat_alphas
    x$res_scissor_rough <- res$res_scissor
    x$alphas_rough <- alphas
    x$file_save <- file_save
    # hierarchical_res_scissor
    return(x)
  }
)

setMethod("step2", signature = c(x = "job_scissor"),
  function(x, alpha_range, k = 8L, alphas = NULL, workers = x$.args$step1$workers, 
    rerun = FALSE, ...)
  {
    step_message("Optimize")
    if (is.null(alphas)) {
      alphas <- .refine_alphas_by_index(
        range(alpha_range), c(1L, 2L), k = k
      )
    }
    res <- .optimize_scissor_with_speed_up_and_score_alpha(
      x, alphas, x$family, rerun = rerun, workers = workers, ...
    )
    x$stat_alphas_hierar <- res$stat_alphas
    x$res_scissor_hierar <- res$res_scissor
    x$alphas_hierar <- alphas
    alls <- dplyr::bind_rows(
      x$stat_alphas_rough, x$stat_alphas_hierar
    )
    alls <- dplyr::select(alls, -dplyr::starts_with("score"))
    x$stat_alphas_all <- .score_alpha_from_scissor(alls)
    t.stat_alphas_all <- set_lab_legend(
      x$stat_alphas_all,
      glue::glue("{x@sig} data scissor alpha selection"),
      glue::glue("Scissor alpha 选择统计表。")
    )
    x <- tablesAdd(x, t.stat_alphas_all)
    snap_alpha <- readLines(file.path(.expath, "description", "scissor_alpha.md"))
    p.alpha <- .plot_scissor_alphas_selection(x$stat_alphas_all)
    p.alpha <- set_lab_legend(
      wrap(p.alpha, 8, 5),
      glue::glue("{x@sig} scissor alpha selection"),
      glue::glue("Scissor alpha 参数的选择|||{.note_legend_scissor_alpha}")
    )
    x <- plotsAdd(x, p.alpha)
    x <- methodAdd(x, "{snap_alpha}")
    x <- methodAdd(x, "实际应用中，分析将首先从初步拟定的 α 区间开始，随后对细胞选择数波动较大，且趋向于选择较少细胞数的 α 区间以对数尺度增加梯度，最终得到所有 α 值对应的细胞选择数以及各细胞类型对应的 Scissor+ 或 Scissor- 数量，再计算综合得分。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_scissor"),
  function(x, focus = NULL, alpha = NULL, cell_num = NULL, n = 10L, 
    nfold = 10L, qs_nthreads = 5L)
  {
    step_message("Significant test.")
    if (is.null(alpha) || is.null(cell_num)) {
      which <- which.max(x$stat_alphas_all$score)
      use <- as.list(x$stat_alphas_all[ which,  ])
      alpha <- use$alpha
      score <- use$score
      cell_num <- use$ncell_select
      x <- snapAdd(
        x, "如图{aref(x@plots$step2$p.alpha)}，选择具有最高综合得分 (Score ={round(score, 2)}) 的 α 值 (α = {alpha}) 对表型关联分析最终定性。在该 α 值下，⟦mark$red('Scissor 细胞选择数为 {cell_num}，占所有细胞数比例为 {use$ratio_select_vs_all}，Scissor+ 细胞所占 Scissor 细胞选择数的比例为 {use$ratio_pos_vs_select}')⟧。"
      )
      if (!is.null(focus)) {
        snap_focus <- vapply(focus, FUN.VALUE = character(1),
          function(name) {
            ratio_pos <- use[[ glue::glue("ratio_pos_vs_select_by_{name})") ]]
            if (is.null(ratio_pos)) {
              stop('is.null(ratio_pos), can not extract of celltype: ', name)
            }
            ratio_neg <- use[[ glue::glue("ratio_neg_vs_select_by_{name})") ]]
            glue::glue(
              "{name} 有 {ratio_pos} 的比例被选择为 Scissor+ 细胞，有 {ratio_neg} 的比例被选择为 Scissor- 细胞"
            )
          })
        x <- snapAdd(x, "其中，⟦mark$red('{bind(snap_focus, co = '；')}')⟧。")
      }
    }
    args <- .qload_multi(x$file_save, nthreads = qs_nthreads)
    fun_test <- function(...) {
      .replace_fun_diag_for_scissor(x$family)
      e(Scissor::reliability.test(
          args$X, args$Y, args$network, alpha = alpha,
          family = x$family, cell_num = cell_num, n = n, nfold = nfold
          ))
    }
    x$res_test <- expect_local_data(
      "tmp", "scissor_reliability", fun_test,
      list(
        colnames(object(x)$bulk_dataset),
        colnames(object(x)$sc_dataset),
        object(x)$phenotype,
        alpha, n, nfold
      )
    )
    return(x)
  })

.refine_alphas_by_index <- function(alphas, index_range, k = 4L) {
  from <- alphas[ index_range[1] ]
  to  <- alphas[ index_range[2] ]
  alphas <- exp(seq(log(from), log(to), length.out = k + 2L))
  alphas[-c(1L, k + 2L)]
}

.optimize_scissor_with_speed_up_and_score_alpha <- function(x, alphas,
  family = x$family, file_save = x$file_save, workers = 4L,
  qs_nthreads = 5L, rerun = FALSE)
{
  if (!is(x, "job_scissor")) {
    stop('!is(x, "job_scissor").')
  }
  fun_scissor <- function(...) {
    if (!file.exists(file_save)) {
      args <- .prepare_scissor_X_Y_network_KeepSparse(
        object(x)$bulk_dataset, 
        object(x)$sc_dataset, 
        object(x)$phenotype,
        family = x$family,
        Save_file = file_save,
        qs_nthreads = qs_nthreads
      )
    } else {
      args <- .qload_multi(file_save, nthreads = qs_nthreads)
    }
    options(future.globals.maxSize = Inf)
    old_plan <- future::plan()
    future::plan(future::multicore, workers = workers)
    on.exit(future::plan(old_plan))
    lst_alphas <- grouping_vec2list(alphas, workers, TRUE)
    lst_res <- lapply(lst_alphas,
      function(alphas) {
        message(glue::glue("Run scissor with alpha: {bind(alphas)}"))
        cli::cli_alert_info("future.apply::future_lapply")
        res_lst <- future.apply::future_lapply(
          alphas, future.seed = TRUE,
          function(alpha) {
            .run_scissor_with_X_Y_network_KeepSparse(
              args$X, args$Y, args$network, alpha = alpha, 
              family = x$family
            )
          }
        )
        setNames(res_lst, paste0("alpha_", alphas))
      })
    res_scissor <- unlist(lst_res, recursive = FALSE)
  }
  res_scissor <- expect_local_data(
    "tmp", "scissor", fun_scissor, rerun = rerun,
    list(
      colnames(object(x)$bulk_dataset),
      colnames(object(x)$sc_dataset),
      object(x)$phenotype,
      alphas
    )
  )
  fun_stat <- function(cells, filter = NULL) {
    meta <- x$metadata
    if (!is.null(filter)) {
      meta <- dplyr::filter(meta, cell %in% !!filter)
    }
    stat <- prop.table(
      table(meta[[ x$group.by ]], meta$cell %in% cells), 1
    )
    stat[, which(colnames(stat) == "TRUE")]
  }
  stat_alphas <- lapply(res_scissor, 
    function(res) {
      cell_pos <- res$Scissor_pos
      cell_neg <- res$Scissor_neg
      cell_select <- c(cell_pos, cell_neg)
      ncell_select <- length(cell_select)
      ncell_all <- length(res$Coefs)
      stats_id <- lapply(c("pos", "neg"),
        function(i) {
          lapply(c("select", "all"),
            function(j) {
              cell_nume <- get(glue::glue("cell_{i}"))
              global <- setNames(
                length(cell_nume) / get(glue::glue("ncell_{j}")),
                glue::glue("ratio_{i}_vs_{j}")
              )
              local <- fun_stat(
                cell_nume,
                switch(j, select = cell_select, all = NULL)
              )
              names(local) <- glue::glue("ratio_{i}_vs_{j}_by_{names(local)}")
              c(global, local)
            })
        })
      stats_select <- fun_stat(cell_select)
      names(stats_select) <- glue::glue("RATIO_SELECT_VS_ALL_BY_{names(stats_select)}")
      c(
        alpha = res$para$alpha,
        ncell_select = ncell_select,
        ratio_select_vs_all = ncell_select / ncell_all,
        ratio_pos_vs_select = length(cell_pos) / ncell_select,
        unlist(stats_select),
        unlist(stats_id)
      )
    })
  stat_alphas <- dplyr::bind_rows(stat_alphas)
  stat_alphas <- .score_alpha_from_scissor(stat_alphas)
  list(stat_alphas = stat_alphas, res_scissor = res_scissor)
}

.score_alpha_from_scissor <- function(data,
  col_ratio = "ratio_select_vs_all",
  pattern_comp = "^RATIO_SELECT_VS_ALL_BY_",
  w_stability = 0.3, w_ratio = 0.5, w_dominance = 0.2)
{
  n <- nrow(data)
  if (n < 2) {
    stop('n < 2, too less "alpha" for scoring.')
  }
  # sort by alpha descending
  data <- data[order(data$alpha, decreasing = TRUE), ]
  comp_cols <- colnames(data)[ grpl(colnames(data), pattern_comp) ]
  # -------------------------
  # 1. stability score
  # -------------------------
  delta <- rep(NA_real_, n)
  for (i in 2:n) {
    delta[i] <- sum(abs(
      as.numeric(data[i, comp_cols]) -
      as.numeric(data[i - 1, comp_cols])
    ))
  }
  delta[1] <- max(delta, na.rm = TRUE)
  stability <- 1 - delta / max(delta, na.rm = TRUE)
  # -------------------------
  # 2. adaptive ratio score
  # -------------------------
  target_ratio <- median(data[[ col_ratio ]], na.rm = TRUE)
  ratio_score <- 1 - abs(data[[ col_ratio ]] - target_ratio)
  # -------------------------
  # 3. dominance score
  # -------------------------
  dom <- apply(data[, comp_cols], 1, function(x) sum(x^2))
  dominance <- dom / max(dom)
  # -------------------------
  # 4. final score
  # -------------------------
  score <- w_stability * stability +
           w_ratio * ratio_score +
           w_dominance * dominance
  data$score_stability <- stability
  data$score_ratio <- ratio_score
  data$score_dominance <- dominance
  data$score <- score
  dplyr::relocate(
    data, alpha, score, dplyr::starts_with("score")
  )
}

.replace_fun_diag_for_scissor <- function(family) {
  require(Matrix)
  fun_name <- switch(
    family, "binomial" = "LogL0", "cox" = "CoxL0", "gaussian" = "LmL0"
  )
  fun_internal <- get(
    fun_name, envir = asNamespace("Scissor")
  )
  body(fun_internal) <- as.call(
    append(
      as.list(body(fun_internal)),
      substitute(diag <- Matrix::diag),
      after = 1L
    )
  )
  replaceFunInPackage(fun_name, fun_internal, "Scissor")
}

.run_scissor_with_X_Y_network_KeepSparse <- function(X, Y, network, 
  alpha, family, cutoff = .2)
{
  if (!is(network, "sparseMatrix")) {
    stop('!is(network, "sparseMatrix").')
  }
  set.seed(123L)
  .replace_fun_diag_for_scissor(family)
  for (i in seq_along(alpha)) {
    cli::cli_alert_info("fit0: Scissor::APML1(...)")
    fit0 <- Scissor::APML1(
      X, Y, family = family, penalty = "Net",
      alpha = alpha[i], Omega = network, nlambda = 100, 
      nfolds = min(10L, nrow(X))
    )
    cli::cli_alert_info("fit1: Scissor::APML1(...)")
    fit1 <- Scissor::APML1(
      X, Y, family = family, penalty = "Net", 
      alpha = alpha[i], Omega = network, lambda = fit0$lambda.min
    )
    if (family == "binomial") {
      Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
    } else {
      Coefs <- as.numeric(fit1$Beta)
    }
    Cell1 <- colnames(X)[which(Coefs > 0)]
    Cell2 <- colnames(X)[which(Coefs < 0)]
    percentage <- (length(Cell1) + length(Cell2))/ncol(X)
    print(sprintf("alpha = %s", alpha[i]))
    print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
        length(Cell1), length(Cell2)))
    print(sprintf("The percentage of selected cell is: %s%%", 
        formatC(percentage * 100, format = "f", digits = 3)))
    if (percentage < cutoff) {
      break
    }
  }
  list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
      family = family), Coefs = Coefs, Scissor_pos = Cell1, 
    Scissor_neg = Cell2, cell_num = length(Cell1) + length(Cell2))
}

.prepare_scissor_X_Y_network_KeepSparse <- function(
  bulk_dataset, sc_dataset, phenotype,
  family = c("binomial", "cox", "gaussian"),
  Save_file = "scissor_inputs.qs", save = TRUE, qs_nthreads = 5L, ...)
{
  family <- match.arg(family)
  common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
  if (!length(common)) {
    stop("There is no common genes between the given single-cell and bulk samples.")
  }
  if (is(sc_dataset, "Seurat")) {
    message("Convert `RNA@data` as matrix...")
    sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
    message("For `RNA_snn`, Improved the Scissor code to keep 'RNA_stnn' as sparse matrix!")
    network <- as(sc_dataset@graphs$RNA_snn, "sparseMatrix")
  } else {
    stop("...")
  }
  # the same effect in `Scissor`, but with sparse matrix
  Matrix::diag(network) <- 0
  network@x[ network@x != 0 ] <- 1
  dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, ])
  dataset1 <- e(preprocessCore::normalize.quantiles(dataset0))
  rownames(dataset1) <- rownames(dataset0)
  colnames(dataset1) <- colnames(dataset0)
  Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
  Expression_cell <- dataset1[, (ncol(bulk_dataset) + 1):ncol(dataset1)]
  cli::cli_alert_info("stats::cor")
  X <- cor(Expression_bulk, Expression_cell)
  quality_check <- quantile(X)
  message("Performing quality-check for the correlations")
  message("The five-number summary of correlations:")
  print(quality_check)
  if (quality_check[3] < 0.01) {
    warning("The median correlation between the single-cell and bulk samples is relatively low.")
  }
  if (family == "binomial") {
    Y <- as.numeric(phenotype)
  }
  if (family == "gaussian") {
    Y <- as.numeric(phenotype)
  }
  if (family == "cox") {
    Y <- as.matrix(phenotype)
    if (ncol(Y) != 2) {
      stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
    }
  }
  if (save) {
    .qsave_multi(X, Y, network, file = Save_file, nthreads = qs_nthreads)
  }
  return(list(X = X, Y = Y, network = network))
}

.plot_scissor_alphas_selection <- function(data) {
  # -------------------------
  # Figure 1: score curves
  # -------------------------
  data_score <- tidyr::pivot_longer(
    dplyr::select(
      data, alpha, score, score_stability,
      score_ratio, score_dominance
      ),
    cols = -alpha,
    names_to = "metric",
    values_to = "value"
  )

  p1 <- ggplot(data_score, aes(x = -log(alpha), y = value, color = metric)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(x = NULL, y = "Score", color = "Score Type") +
    theme_bw()

  # -------------------------
  # Figure 2: selected / all
  # -------------------------
  data_select <- tidyr::pivot_longer(
    dplyr::select(
      data, alpha,
      dplyr::starts_with("RATIO_SELECT_VS_ALL_BY_")
      ),
    cols = -alpha,
    names_to = "celltype",
    values_to = "ratio"
  )

  data_select$celltype <- gsub(
    "^RATIO_SELECT_VS_ALL_BY_", "", data_select$celltype
  )

  p2 <- ggplot(data_select, aes(x = -log(alpha), y = ratio, fill = celltype)) +
    geom_area(position = "stack") +
    scale_fill_manual(values = color_set()) +
    labs(x = NULL, y = "Selected / All", fill = "Cell types") +
    theme_bw()

  # -------------------------
  # Figure 3: pos / selected
  # -------------------------
  data_pos <- tidyr::pivot_longer(
    dplyr::select(
      data,
      alpha,
      dplyr::starts_with("ratio_pos_vs_select_by_")
      ),
    cols = -alpha,
    names_to = "celltype",
    values_to = "ratio"
  )

  data_pos$celltype <- gsub(
    "^ratio_pos_vs_select_by_",
    "",
    data_pos$celltype
  )

  p3 <- ggplot(data_pos, aes(x = -log(alpha), y = ratio, fill = celltype)) +
    geom_area(position = "stack") +
    scale_fill_manual(values = color_set()) +
    labs(x = "-log(alpha)", y = "Positive / Selected", fill = "Cell types") +
    theme_bw()

  require(patchwork)
  p.alpha <- (p1 / p2 / p3) +
    patchwork::plot_layout(
      heights = c(1.1, 1.6, 1.6),
      guides = "collect"
    )
}

.note_legend_scissor_alpha <- "Scissor 参数 α 调节过程中综合评分及细胞选择结构的动态变化。横坐标均为参数 α，并采用对数尺度反向排列，表示从左至右模型约束逐渐减弱、被选中细胞数量逐渐增加。上图展示不同 α 条件下综合评分体系及其组成指标的变化趋势。黑线（score）表示综合评分，用于评估各候选 α 的整体表现；其余三条曲线分别表示组成稳定性得分（score_stability）、选择复杂度得分（score_ratio）及群体主导性得分（score_dominance）。其中，组成稳定性反映相邻 α 条件下细胞组成变化是否趋于收敛；选择复杂度反映被选中细胞比例是否处于适中范围；群体主导性反映筛选结果是否由有限优势细胞群主导。综合评分最高处对应推荐的最优 α 取值区间。中图展示各细胞类型被 Scissor 选中的比例占其原始总体细胞数量的变化。堆叠面积图中不同颜色代表不同细胞类型，纵坐标表示各类型细胞被纳入模型的比例。该图用于描述随着 α 逐渐减小，各类细胞被模型逐步释放并进入选择集合的动态过程，可用于识别对参数变化最敏感的细胞群体及主要贡献群体。下图展示各细胞类型在全部被选中细胞中属于 Scissor+ 的比例。纵坐标表示某细胞类型在被选中后，其正向关联细胞（Scissor+）所占比例。该图用于评估不同 α 条件下 phenotype 正相关信号的细胞来源构成，以及各细胞类型在正向生物学效应中的相对贡献。综合三幅图可同时评估参数 α 对模型性能、筛选规模及细胞组成结构的影响，从而实现对最优 α 的客观选择，并进一步解析 phenotype 相关细胞群体随参数变化的层级释放过程。"

