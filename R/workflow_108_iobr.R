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
    analysis = "IOBR 肿瘤免疫微环境分析"
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
  function(x, use = .guess_symbol(x))
  {
    snap <- snap(x)
    object <- extract_unique_genes.job_limma(x)
    mtx <- object$E
    rownames(mtx) <- gname(object$genes[[ use ]])
    x <- snapAdd(
      job_iobr(mtx, metadata = object$targets), 
      "以数据集 ({x$project}, dataset: {x@sig}) 进行 IOBR 免疫微环境评估。"
    )
    x <- snapAdd(x, snap)
    return(x)
  })

setMethod("step0", signature = c(x = "job_iobr"),
  function(x){
    step_message("Prepare your data with function `job_iobr`.")
  })

setMethod("step1", signature = c(x = "job_iobr"),
  function(x, method = c("cibersort"), cache = "iobr")
  {
    step_message("Calculating ...")
    method <- match.arg(method)
    require(IOBR)
    dir.create(cache)
    hash <- digest::digest(list(object(x), method))
    file_cache <- add_filename_suffix(
      file.path(cache, "iobr.rds"), hash
    )
    if (is.null(x$res)) {
      if (file.exists(file_cache)) {
        x$res <- readRDS(file_cache)
      } else {
        x$res <- e(IOBR::deconvo_tme(
            eset = object(x), method = method, arrays = FALSE
            ))
      }
    }
    if (!file.exists(file_cache)) {
      saveRDS(x$res, file_cache)
    }
    x <- snapAdd(x, "使用算法为 {method}。")
    x$method <- method
    return(x)
  })

setMethod("step2", signature = c(x = "job_iobr"),
  function(x, group.by = "group")
  {
    step_message("Significant test.")
    data <- x$res
    data <- map(data, "ID", x$metadata, "sample", group.by, col = "group")
    if (length(unique(data$group)) != 2) {
      stop('length(unique(data$group) != 2).')
    }
    data <- e(tidyr::pivot_longer(
      data, c(-ID, -group), names_to = "type", values_to = "level"
    ))
    data <- dplyr::mutate(
      data, type = s(type, glue::glue("_{x$method}$"), "", TRUE)
    )
    data <- e(dplyr::filter(data, !grpl(type, "Correlation|P-value|RMSE")))
    data <- e(dplyr::group_by(data, type))
    fun <- function(level, group) {
      data <- data.frame(level = level, group = group)
      wilcox.test(level ~ group, data = data)$p.value
    }
    dataSig <- e(dplyr::summarise(data, pvalue = fun(level, group)))
    dataSig <- add_anno(dataSig)
    dataSig <- set_lab_legend(
      dataSig,
      glue::glue("{x@sig} {x$method} wilcox test data"),
      glue::glue("为 {x$method} 算法 wilcox test 组间比较附表。")
    )
    p.boxplot <- ggplot(data) +
      geom_jitter(aes(x = type, y = level, fill = group, group = group),
        stroke = 0, shape = 21, color = "transparent",
        position = position_jitterdodge(.2)) +
      geom_boxplot(
        aes(x = type, y = level, color = group),
        outlier.shape = NA, fill = "transparent", notchwidth = .7) +
      geom_text(data = dataSig, aes(x = type, label = sig, y = max(data$level))) +
      scale_color_manual(values = color_set()) +
      scale_fill_manual(values = color_set()) +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_blank()
    p.boxplot <- set_lab_legend(
      wrap(p.boxplot, 11, 5),
      glue::glue("{x@sig} {x$method} Immune infiltration"),
      glue::glue("为 {x$method} 算法的免疫微环境分析箱形图 (wilcox.test)。")
    )
    x <- plotsAdd(x, p.boxplot = p.boxplot)
    x <- tablesAdd(x, t.dataSig = dataSig)
    s.com <- dplyr::filter(dataSig, pvalue < .05)$type
    x <- snapAdd(x, "组间差异分析，有显著区别的是 {bind(s.com)}")
    return(x)
  })

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
