# ==========================================================================
# workflow of scMeta
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scMeta <- setClass("job_scMeta", 
  contains = c("job"),
  prototype = prototype(
    pg = "scMeta",
    info = c("https://github.com/wu-yc/scMetabolism"),
    cite = "",
    method = "",
    tag = "scMeta",
    analysis = "scMetabolism 代谢通路活性评估"
    ))

setGeneric("asjob_scMeta",
  function(x, ...) standardGeneric("asjob_scMeta"))

setMethod("asjob_scMeta", signature = c(x = "job_seurat"),
  function(x, ref = NULL, workers = 5, method = "AUCell", type = "KEGG")
  {
    if (is.null(ref)) {
      fea <- as_feature(levels(Seurat::Idents(object(x))), "关键细胞", "cell")
      ref <- snap(fea)
    }
    snapAdd_onExit("x", "将{snap(fea)}进行通路富集分析。")
    fun_compute <- function(...) {
      object <- e(scCustomize::Convert_Assay(
          object(x), "RNA", convert_to = "V3"
          ))
      object <- e(scMetabolism::sc.metabolism.Seurat(
          obj = object, method = method, imputation = FALSE,
          ncores = workers, metabolism.type = type
          ))
      object@assays$METABOLISM$score
    }
    data <- expect_local_data(
      "tmp", "scMeta", fun_compute, list(method, type, colnames(object(x)))
    )
    params <- params(x)
    x <- .job_scMeta(object = data)
    x@params <- append(x@params, params)
    x$method <- method
    x$type <- type
    return(x)
  })

setMethod("step0", signature = c(x = "job_scMeta"),
  function(x){
    step_message("Prepare your data with function `job_scMeta`.")
  })

setMethod("step1", signature = c(x = "job_scMeta"),
  function(x, use.p = c("p.value", "p.adjust"), show = 10L){
    step_message("Draw plot.")
    use.p <- match.arg(use.p)
    t.wilcox <- .get_test_list_by_compare_group_in_matrix(
      object(x), x$metadata, x$group.by, compare.by = "group"
    )
    t.wilcox <- set_lab_legend(
      t.wilcox,
      glue::glue("{x@sig} diff {x$type} pathway {x$method} activity scMetabolism"),
      glue::glue("{x$type} 活性差异分析")
    )
    p.hp <- .set_heatmap_with_test_list(t.wilcox, show, x$group.by)
    p.hp <- set_lab_legend(
      p.hp,
      glue::glue("{x@sig} Top {show} diff {x$type} heatmap scMetabolism"),
      glue::glue("scMetabolism 富集结果中各细胞类型按显著性排列靠前的通路|||横轴分布的是不同细胞类型，对应于各自组别；纵向为富集通路的名称。")
    )
    x <- plotsAdd(x, p.hp)
    snap <- .stat_table_by_pvalue(
      t.wilcox, n = show, use.p = use.p, colName = "Name"
    )
    x <- methodAdd(x, "scMetabolism 是一种基于单细胞转录组数据评估细胞代谢通路活性的计算方法，其主要目的是在单细胞分辨率下刻画不同细胞群体的代谢功能状态及其异质性。该方法通过整合基因表达数据与预定义的代谢通路数据库（如 {x$type}），对每个细胞或细胞群体进行通路活性打分，从而量化各类代谢过程的相对活跃程度。")
    x <- methodAdd(x, "使用 R 包 `scMetabolism` ⟦pkgInfo('scMetabolism')⟧ 基于基因表达量对 {x$type} 通路进行 {x$method} 活性评分。")
    x <- methodAdd(x, "以 wilcox.test 检验评估每个通路在组间的表达差异显著性，并以 {detail(use.p)} 对其排序。")
    x <- snapAdd(x, "通过对 scMetabolism 活性评分差异分析 {x$type} 富集通路 (⟦mark$blue('{detail(use.p)} &lt; 0.05')⟧)，一共富集到{snap}\n\n")
    x <- tablesAdd(x, t.wilcox)
    return(x)
  })

setMethod("step2", signature = c(x = "job_scMeta"),
  function(x){
    step_message("")
    
    return(x)
  })

.set_heatmap_with_test_list <- function(data, top, group.by, groups = NULL, col = "Name")
{
  vapply_mean <- function(x) {
    vapply(x, mean, double(1))
  }
  if (is.null(groups)) {
    groups <- colnames(data)[vapply(data, function(x) is(x, "list"), logical(1))]
  }
  fea <- head(unique(data[[col]]), n = top)
  data <- dplyr::filter(data, !!rlang::sym(col) %in% !!fea)
  data <- dplyr::select(data, -p.value, -p.adjust)
  data <- dplyr::mutate(
    data, dplyr::across(dplyr::all_of(groups), vapply_mean)
  )
  data <- tidyr::pivot_longer(
    data, dplyr::all_of(groups), 
    names_to = "group", values_to = "value"
  )
  data <- dplyr::mutate(
    data, Cell = paste0(group, "_", !!rlang::sym(group.by)),
    value = scale(log2(value + 1), center = TRUE)[, 1]
  )
  args <- list(
    .data = data, .row = rlang::quo(!!rlang::sym(col)), .column = rlang::quo(Cell),
    .value = rlang::quo(value), cluster_columns = FALSE, column_names_rot = 90,
    cluster_rows = TRUE, group_by = rlang::quo(group),
    row_names_max_width = grobWidth(textGrob(rownames(data), gpar(fontsize = 10, fontface = 1)))
  )
  wrap_scale_heatmap(
    # ComplexHeatmap::Heatmap
    funPlot(heatmap_with_group, args),
    data$Cell, data[[col]], pre_width = max(nchar(data[[col]])) * .1 + 2,
    pre_height = max(nchar(data[[col]])) * .05 + 1,
    w.size = .3
  )
}

.get_test_list_by_compare_group_in_matrix <- function(data, meta, group.by,
  compare.by = "group", sample.by = "cell", idcol = "Name")
{
  if (!is.null(idcol)) {
    data <- as_tibble(data, idcol = idcol)
  }
  data <- tidyr::pivot_longer(
    data, dplyr::where(is.double), names_to = "type", values_to = "value"
  )
  data <- map(data, "type", meta, sample.by, group.by, col = group.by)
  data <- map(data, "type", meta, sample.by, compare.by, col = compare.by)
  # drop the 'sample' columns, for pivot_wider gathered to generate 'list'
  data <- dplyr::select(data, -type)
  data <- tidyr::pivot_wider(data, names_from = !!rlang::sym(compare.by), values_from = value)
  groups <- unique(meta[[compare.by]])
  cli::cli_alert_info("wilcox.test")
  fun_test <- function(x, y) {
    pbapply::pbvapply(seq_along(x),
      function(n) {
        wilcox.test(x[[n]], y[[n]])$p.value
      }, double(1))
  }
  data <- dplyr::mutate(
    data, p.value = fun_test(!!rlang::sym(groups[1]), !!rlang::sym(groups[2]))
  )
  data <- dplyr::group_by(data, !!rlang::sym(group.by))
  data <- dplyr::mutate(
    data, p.adjust = p.adjust(p.value, method = "BH")
  )
  data <- dplyr::ungroup(data)
  dplyr::arrange(data, p.value)
}

setMethod("set_remote", signature = c(x = "job_scMeta"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
