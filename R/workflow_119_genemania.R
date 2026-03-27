# ==========================================================================
# workflow of genemania
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_genemania <- setClass("job_genemania", 
  contains = c("job"),
  prototype = prototype(
    pg = "genemania",
    info = c("https://genemania.org/"),
    cite = "",
    method = "",
    tag = "genemania",
    analysis = "GeneMANIA 网络分析"
    ))


setGeneric("asjob_genemania",
  function(x, ...) standardGeneric("asjob_genemania"))

setMethod("asjob_genemania", signature = c(x = "feature"),
  function(x, ...){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_genemania(object = fea)
    return(x)
  })

setMethod("step0", signature = c(x = "job_genemania"),
  function(x){
    step_message("Prepare your data with function `job_genemania`.")
  })

setMethod("step1", signature = c(x = "job_genemania"),
  function(x, crawler = FALSE, dir_save = "geneMANIA"){
    step_message("Open webpage.")
    x$dir_save <- dir_save
    dir.create(x$dir_save, FALSE)
    if (crawler) {
      x$link <- start_drive(port = port, download.dir = dir_save)
      x$link$open()
      x$url_home <- "https://genemania.org/search"
      gett(object(x))
      x$link$navigate(x$url_home)
    }
    x <- methodAdd(x, "使用 GeneMANIA（https://genemania.org/）进行检索和创建功能相关的基因列表，并构建交互式网络，进行可视化 GeneMANIA 相互作用网络。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_genemania"),
  function(x, n_show = 7){
    step_message("Load file.")
    # files <- c("genemania-genes.txt", "genemania-interactions.txt", "genemania-functions.txt")
    files <- c("genemania-functions.txt")
    files <- file.path(x$dir_save, files)
    res_alls <- sapply(files, simplify = FALSE,
      function(file) {
        if (!file.exists(file)) {
          stop(glue::glue("Can not found: {file}"))
        }
        ftibble(file)
      })
    names(res_alls) <- s(tools::file_path_sans_ext(basename(names(res_alls))), 'genemania-', '')
    x$res_alls <- res_alls
    nOthers <- 20L
    x <- snapAdd(x, "关键基因与功能相关的 {nOthers} 个相关基因，它们之间通过 7 种不同的作用方式互作（共同表达、物理交互、遗传相互作用、共享蛋白质域、联合本地化、富集通路、预测）。")
    showFun <- head(res_alls$functions$Function, n = n_show)
    x <- snapAdd(x, "通过 GeneMANIA 计算、预测这些基因与多种功能相互关联，其中最显著 (FDR &lt; 0.05) 的 {n_show} 种为：{atrans(showFun)}。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_genemania"),
  function(x, rerun = FALSE){
    step_message("Render Picture.")
    file <- file.path(x$dir_save, "genemania-report.pdf")
    pngFile <- paste0(tools::file_path_sans_ext(file), ".png")
    png <- pdf_convert(file, filenames = pngFile, pages = 1, dpi = 300)
    fun_split_and_append <- function(png) {
      message("Click to split main figure with legends.")
      # top and bottom
      figure_and_legends <- split_image_by_click(png)
      figure <- split_image_by_click(figure_and_legends[[1]])[[2]]
      message("Click to split Legends with others (pagenumber).")
      legends <- split_image_by_click(figure_and_legends[[2]])[[1]]
      message("Click to split Legends.")
      legends <- split_image_by_click(legends)
      legend <- bind_image_add_spacer(legends[[1]], legends[[2]], stack = TRUE)
      p.network <- bind_image_add_spacer(figure, legend)
      grid::rasterGrob(as.raster(p.network))
    }
    p.network <- expect_local_data(
      "tmp", "genemania", fun_split_and_append, args = list(png = png), rerun = rerun
    )
    p.network <- set_lab_legend(
      wrap(p.network),
      glue::glue("{x@sig} geneMANIA network"),
      glue::glue("GeneMANIA 相互作用网络|||圆形代表不同基因，中间为所选基因 {bind(object(x))}，不同颜色的边代表不同相互作用，外围圆圈大小代表与关键基因的相关强度，圆形内部的不同颜色代表参与了不同功能。")
    )
    x <- plotsAdd(x, p.network)
    return(x)
  })
