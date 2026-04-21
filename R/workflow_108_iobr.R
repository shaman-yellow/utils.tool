# ==========================================================================
# workflow of iobr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_iobr <- setClass("job_iobr",
  contains = c("job"),
  prototype = prototype(
    pg = "iobr",
    info = c("https://iobr.github.io/book/index.html"),
    cite = "[@Enhancing_immun_Zeng_2024]",
    method = "",
    tag = "iobr",
    analysis = "IOBR 免疫浸润分析"
    ))

job_iobr <- function(object, metadata)
{
  x <- .job_iobr(object = object)
  x$metadata <- metadata
  return(x)
}

setGeneric("asjob_iobr", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_iobr"))

setMethod("asjob_iobr", signature = c(x = "job_limma"),
  function(x, use = .guess_symbol(x), project = x$project, 
    levels = NULL, guess_which_level = 1L, ...)
  {
    if (x$rna) {
      stop('x$rna == TRUE, not support now. Use `job_deseq2` ...')
    }
    if (is.null(levels)) {
      levels <- .get_group_from_contrast_character(names(x@tables$step2$tops)[guess_which_level])
    }
    object <- x$normed_data
    if (FALSE) {
      object <- extract_unique_genes.job_limma(x)
    }
    if (range(object$E)[1] < 0) {
      stop('range(object$E)[1] < 0.')
    }
    mtx <- object$E
    rownames(mtx) <- gname(object$genes[[ use ]])
    x <- job_iobr(mtx, metadata = object$targets)
    x$levels <- levels
    return(x)
  })

setMethod("step0", signature = c(x = "job_iobr"),
  function(x){
    step_message("Prepare your data with function `job_iobr`.")
  })

setMethod("step1", signature = c(x = "job_iobr"),
  function(x, method = "cibersort",
    run_all = FALSE, workers = 1,
    cache = "tmp", tumor = FALSE, ..., rerun = FALSE)
  {
    step_message("Calculating ...")
    methods <- c(
      "cibersort",
      "xcell",
      "estimate",
      "timer", "quantiseq", "svr", "lsei",
      "mcpcounter",
      "epic",
      "ips",
      "quantiseq",
      "svr",
      "lsei"
    )
    method <- match.arg(method, methods)
    require(IOBR)
    dir.create(cache, FALSE)
    if (run_all) {
      methods <- methods
    } else {
      methods <- method
    }
    args <- list(...)
    args$tumor <- tumor
    args$eset <- object(x)
    set.seed(x$seed)
    x$allres <- pbapply::pbsapply(
      methods, simplify = FALSE, cl = workers,
      function(method) {
        args$method <- method
        try(expect_local_data(cache, "iobr", IOBR::deconvo_tme, args, rerun = rerun))
      }
    )
    x$res <- x$allres[[ method ]]
    x <- methodAdd(x, "为系统评估诊断基因与免疫微环境之间的潜在联系，进一步开展免疫浸润及相关性分析，筛选与免疫微环境密切相关的细胞亚群，可从“基因–细胞互作”角度深化对疾病发生发展机制的理解。")
    x <- methodAdd(x, "以 R 包 `IOBR` ⟦pkgInfo('IOBR')⟧ 选择算法 {method} 对数据集免疫浸润分析。")
    if (isNamespaceLoaded("IOBR")) {
      detach("package:IOBR")
    }
    x$method <- method
    return(x)
  })

setMethod("step2", signature = c(x = "job_iobr"),
  function(x, group.by = "group", cut.p = .05, cut.cor = .3, 
    add_noise = TRUE, keep_all = if (x$method == "xcell") FALSE else TRUE, method_cor = "spearman")
  {
    step_message("Significant test.")
    if (!is.null(x$allres)) {
      lst <- x$allres
    } else {
      lst <- setNames(list(data), x$method)
    }
    x$all_filter <- lst <- lapply(lst, 
      function(data) {
        data <- dplyr::rename_with(data, ~ sub("_[^_]+$", "", .x), -ID)
        data <- map(data, "ID", x$metadata, "sample", group.by, col = "group")
        if (x$method == "cibersort") {
          if (!is.null(cut.p)) {
            data <- trace_filter(data, `P-value` < !!cut.p)
            snap_filter <- snap(data)
            snap <- glue::glue("根据 cibersort 算法得到的 P-value，{snap_filter}{try_snap(data$group)}。")
          } else {
            snap <- NULL
          }
          data <- dplyr::select(data, -Correlation, -`P-value`, -RMSE)
        }
        namel(data, snap)
      })
    data <- mtx <- lst[[ x$method ]]$data
    if (x$method == "cibersort" && !is.null(cut.p)) {
      x <- snapAdd(x, "{lst[[ x$method ]]$snap}")
      x <- methodAdd(x, "剔除 P-value &gt; {cut.p} 的低可靠性样本。")
    }
    if (length(unique(data$group)) != 2) {
      stop('length(unique(data$group) != 2).')
    }
    data <- e(tidyr::pivot_longer(
        data, c(-ID, -group), names_to = "type", values_to = "level"
        ))
    p.stack <- ggplot(data, aes(x = ID, y = level, fill = type)) +
      geom_bar(stat = "identity", width = .7) +
      facet_grid(~ group, scales = "free_x", space = "free") +
      theme_light() +
      labs(y = "Relative Proportion", x = "Sample", fill = "Type") +
      theme(axis.text.x = element_blank()) +
      scale_fill_manual(values = color_set(TRUE))
    p.stack <- wrap_scale(
      p.stack, length(unique(data$ID)), 
      10, pre_width = 3, pre_height = 1, size = .3
    )
    p.stack <- set_lab_legend(
      p.stack,
      glue::glue("{x@sig} Infiltration Landscape"),
      glue::glue("免疫细胞渗透比例")
    )
    data <- e(dplyr::group_by(data, type))
    fun_pvalue <- function(level, group) {
      data <- data.frame(level = level, group = group)
      wilcox.test(level ~ group, data = data)$p.value
    }
    fun_measure <- function(level, group) {
      gp <- split(level, group)
      ms <- vapply(gp, function(x) fivenum(x)[3], double(1))
      names(gp)[which.max(ms)]
    }
    groupCor <- e(dplyr::summarise(data,
        pvalue = fun_pvalue(level, group), higher = fun_measure(level, group)))
    groupCor <- add_anno(groupCor)
    groupCor <- set_lab_legend(
      groupCor,
      glue::glue("{x@sig} {x$method} wilcox test data"),
      glue::glue("为 {x$method} 算法 wilcox test 组间比较附表。")
    )
    if (x$method == "xcell") {
      data <- dplyr::mutate(data, level = log2(level))
      ylab <- "log2(level)"
      exSnap <- " 对 xCell 算法得出的富集分数检验显著性，而图中的箱形图为了展示于同一尺度，做了 log2 变换"
    } else {
      ylab <- "level"
      exSnap <- ""
    }
    p.boxplot <- ggplot(data) +
      geom_jitter(aes(x = type, y = level, fill = group, group = group),
        stroke = 0, shape = 21, color = "transparent",
        position = position_jitterdodge(.2)) +
      geom_boxplot(
        aes(x = type, y = level, color = group),
        outlier.shape = NA, fill = "transparent", notchwidth = .7) +
      geom_text(data = groupCor, aes(x = type, label = sig, y = max(data$level))) +
      labs(y = ylab) +
      scale_color_manual(values = color_set()) +
      scale_fill_manual(values = color_set()) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_blank()
    p.boxplot <- set_lab_legend(
      wrap(p.boxplot, 11, 5),
      glue::glue("{x@sig} {x$method} Immune infiltration"),
      glue::glue("为 {x$method} 算法的免疫微环境分析箱形图|||使用 wilcox.test{exSnap} ({.md_p_significant}) 。")
    )
    dataSig <- dplyr::filter(groupCor, pvalue < cut.p)
    s.com <- dataSig$type
    feature(x) <- as_feature(s.com, "IOBR 免疫浸润分析差异细胞", nature = "cells")
    x <- snapAdd(x, "以 wilcox.test 组间差异分析{aref(p.boxplot)}，有显著区别 (⟦mark$blue('p &lt; {cut.p}')⟧) 的差异细胞有 {nrow(dataSig)} 类。")
    fun_snap <- function(higher) {
      ifelse(higher == x$levels[1], "升高", "下降")
    }
    snap_sig <- glue::glue("{s.com} 活性显著 {fun_snap(dataSig$higher)}")
    x <- snapAdd(x, "相比于 {x$levels[2]} 组，{x$levels[1]} 组的{bind(snap_sig)}。\n\n\n")
    mtx <- dplyr::select(mtx, -ID, -group)
    if (!keep_all) {
      mtx <- dplyr::select(mtx, dplyr::all_of(s.com))
    }
    if (add_noise) {
      mtx <- add_noise(mtx, seed = x$seed)
    }
    x$data_noise <- mtx
    # It's better to use it directly... oh my god!!!  `ggcor::fortify_cor`
    res_cor.test <- psych::corr.test(
      mtx, method = method_cor
    )
    x <- methodAdd(x, "以 R 包 `psych` ⟦pkgInfo('psych')⟧ 对免疫浸润细胞之间进行关联分析 ({method_cor}) 。⟦mark$blue('将|cor| &gt; 0.3 且 p &lt; {cut.p} 的分析结果判定为具有统计学意义')⟧。")
    t.cells_cor <- add_anno(.corp(as_data_long(
      res_cor.test$r, res_cor.test$p, "cells_x", "cells_y", 
      "cor", "pvalue"
    )))
    t.cells_cor <- set_lab_legend(
      t.cells_cor,
      glue::glue("{x@sig} cells correlation"),
      glue::glue("免疫浸润细胞相关性分析附表。")
    )
    colors <- ifelse(
      colnames(res_cor.test$p) %in% s.com, "red", "grey"
    )
    fun_plot <- function(x) {
      p <- ggcor::quickcor(ggcor::as_cor_tbl(x), type = "lower", cor.test = TRUE, method = method_cor)
      .ggcor_add_general_style(p, NULL)
    }
    p.cor <- wrap_scale(
      fun_plot(res_cor.test), 
      ncol(mtx), ncol(mtx), pre_height = 2, size = .5
    )
    if (keep_all) {
      exSnap <- "热图的细胞名若显示为红色，则为差异细胞；"
    } else {
      exSnap <- "展示差异细胞之间的关联性；"
    }
    p.cor <- set_lab_legend(
      p.cor,
      glue::glue("{x@sig} Correlation immune cells"),
      glue::glue("免疫细胞相关性分析热图|||{exSnap}热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。")
    )
    s.cc <- dplyr::filter(
      t.cells_cor, cells_x %in% s.com, cells_y %in% s.com,
      cells_x != cells_y,
      pvalue < cut.p, abs(cor) > cut.cor
    )
    if (nrow(s.cc)) {
      snap_cor <- .stat_correlation_table(
        s.cc, "cells_x", "cells_y"
      )
      x <- snapAdd(x, "免疫细胞之间的相关性分析结果{aref(p.cor)}，{bind(snap_cor)}。\n\n\n")
    } else {
      x <- snapAdd(x, "差异免疫浸润细胞之间未发现显著关联。 ")
    }
    x <- plotsAdd(x, p.boxplot, p.stack, p.cor)
    x <- tablesAdd(x, t.groupCor = groupCor, t.cells_cor)
    return(x)
  })

.ggcor_add_general_style <- function(p, type = "lower") {
  if (!is.null(type)) {
    p + ggcor::geom_square(data = ggcor::get_data(type = type, show.diag = FALSE)) +
      ggcor::geom_mark(data = ggcor::get_data(type = type, show.diag = FALSE), sep = "\n") +
      .scale_for_cor_palette()
  } else {
    p + ggcor::geom_square() +
      ggcor::geom_mark(sep = "\n") +
      .scale_for_cor_palette()
  }
}

.scale_for_cor_palette <- function(type = c("fill", "color")) {
  type <- match.arg(type)
  if (type == "fill") {
    scale_fill_distiller(
      palette = "RdBu",
      direction = -1,
      limits = c(-1, 1),
      name = "Correlation"
    )
  } else {
    scale_color_distiller(
      palette = "RdBu",
      direction = -1,
      limits = c(-1, 1),
      name = "Correlation"
    )
  }
}

add_noise <- function(mtx, level = 1e-6, seed = 12345) {
  var_values <- apply(mtx, 2, var, na.rm = TRUE)
  zero_var_cols <- names(var_values[var_values == 0 | is.na(var_values)])
  if (length(zero_var_cols)) {
    message("Column with variance of 0: ", bind(zero_var_cols))
    set.seed(seed)
    for (col in zero_var_cols) {
      # Add minimal random noise
      mtx[, col] <- mtx[, col] + rnorm(nrow(mtx), 0, level)
    }
  }
  return(mtx)
}

setMethod("step3", signature = c(x = "job_iobr"),
  function(x, ref, recode = NULL, cut.p = .05, cut.cor = .3, 
    add_noise = TRUE, keep_all = FALSE, 
    use_vst = FALSE, method_cor = "spearman")
  {
    step_message("Correlation for cells and genes.")
    mtx <- data.frame(dplyr::select(x$all_filter[[ x$method ]]$data, -group))
    rownames(mtx) <- mtx$ID
    data.cell <- t(as.matrix(dplyr::select(mtx, -ID)))
    if (!keep_all) {
      data.cell <- data.cell[ rownames(data.cell) %in% feature(x), ]
    }
    if (use_vst) {
      if (is.null(x$vst)) {
        stop('is.null(x$vst).')
      }
      data.expr <- x$vst
    } else {
      data.expr <- object(x)
    }
    if (!is.null(recode)) {
      if (!is(recode, "list")) {
        recode <- as.list(recode)
      }
      rownames(data.expr) <- dplyr::recode(
        rownames(data.expr), !!!recode, .default = rownames(data.expr)
      )
    }
    if (is(ref, "feature")) {
      ref.snap <- snap(ref)
    } else {
      ref.snap <- "基因"
    }
    if (is.list(ref)) {
      ref <- unlist(ref)
    }
    data.expr <- data.expr[rownames(data.expr) %in% ref, ]
    if (any(isNot <- !ref %in% rownames(data.expr))) {
      stop(glue::glue("Not cover: {bind(ref[isNot])}"))
    }
    data.expr <- data.expr[, match(colnames(data.cell), colnames(data.expr))]
    data.expr <- t(data.expr)
    data.cell <- t(data.cell)
    if (add_noise) {
      data.cell <- add_noise(data.cell)
    }
    # It's better to use `ggcor::fortify_cor` directly... oh my god!!!  
    dataCor <- psych::corr.test(data.expr, data.cell, method = method_cor)
    x <- methodAdd(x, "以 R 包 `psych` ⟦pkgInfo('psych')⟧ 对基因与免疫浸润细胞之间进行关联分析 ({method_cor}) 。将|cor| &gt; 0.3 且 p &lt; {cut.p} 的分析结果判定为具有统计学意义。")
    colors <- ifelse(
      colnames(dataCor$p) %in% feature(x), "red", "grey"
    )
    fun_plot <- function(x) {
      p <- ggcor::quickcor(ggcor::as_cor_tbl(x), type = "lower", cor.test = TRUE, method = method_cor)
      .ggcor_add_general_style(p)
    }
    p.GeneCellCor <- fun_plot(dataCor)
    p.GeneCellCor <- wrap_scale(
      p.GeneCellCor, ncol(data.cell), ncol(data.expr), min_height = 3
    )
    if (keep_all) {
      exSnap <- "热图的细胞名若显示为红色，则为差异细胞；"
    } else {
      exSnap <- "展示差异细胞之间的关联性；"
    }
    p.GeneCellCor <- set_lab_legend(
      p.GeneCellCor,
      glue::glue("{x@sig} correlation of Immune cells and selected genes"),
      glue::glue("{ref.snap}和免疫细胞相关性分析|||{exSnap}热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。")
    )
    x <- plotsAdd(x, p.GeneCellCor)
    t.geneCellCor <- add_anno(
      .corp(as_data_long(dataCor$r, dataCor$p, "genes", "cells", "cor", "pvalue"))
    )
    t.geneCellCor <- set_lab_legend(
      t.geneCellCor,
      glue::glue("{x@sig} correlation between cells and genes"),
      glue::glue("{ref.snap}与免疫细胞相关性分析")
    )
    s.gc <- dplyr::filter(
      t.geneCellCor, cells %in% feature(x),
      pvalue < cut.p, abs(cor) > cut.cor
    )
    x <- tablesAdd(x, t.geneCellCor)
    if (nrow(s.gc)) {
      snap_cor <- .stat_correlation_table(s.gc, "cells", "genes")
      x <- snapAdd(x, "免疫细胞与基因之间的关联性分析表明{aref(p.GeneCellCor)}，{bind(snap_cor)}。\n\n\n")
    } else {
      x <- snapAdd(x, "差异免疫浸润细胞与{ref.snap}之间未发现显著关联{aref(p.GeneCellCor)}。 ")
    }
    return(x)
  })

.stat_correlation_table <- function(data, x, y, cor = "cor", 
  pvalue = "pvalue", label.x = "", label.y = "", maxShow = 5)
{
  data <- dplyr::filter(data, !!rlang::sym(pvalue) < .05)
  if (!nrow(data)) {
    return(glue::glue("无显著关联"))
  }
  leader <- ""
  if (nrow(data) > maxShow) {
    leader <- glue::glue("共得到 {nrow(data)} 对显著关联，按相关系数从大到小排序，前 {maxShow} 的显著关联中，")
    data <- head(data, n = maxShow)
  }
  fun_measure <- function(x) ifelse(x > 0, "正", "负")
  fmt <- function(x) signif(x, 3)
  snap <- glue::glue(
    "{label.x} {data[[x]]} 和 {label.y} {data[[y]]} 之间显著{fun_measure(data[[cor]])}相关 (cor = {fmt(data[[cor]])}, p = {fmt(data[[pvalue]])})"
  )
  glue::glue("{leader}{bind(snap)}")
}

setMethod("add_anno", signature = c(x = "df"),
  function(x){
    if (!any(colnames(x) == "pvalue")) {
      stop('!any(colnames(x) == "pvalue").')
    }
    x <- dplyr::mutate(x,
      sig = ifelse(
        pvalue < .001, "***", ifelse(
          pvalue < .01, "**", ifelse(pvalue < .05, "*", "")
        )
      )
    )
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_iobr"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
