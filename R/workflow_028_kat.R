# ==========================================================================
# workflow of kat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_kat <- setClass("job_kat", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/navinlabcode/copykat"),
    cite = "[@DelineatingCopGaoR2021]",
    method = "R package `copyKAT` used for aneuploid cell or cancer cell prediction",
    tag = "scrna:cancer",
    analysis = "CopyKAT 癌细胞鉴定"
    ))

job_kat <- function(x, use = names(x@object@assays)[[1]], layer = "counts")
{
  assay <- x@object@assays[[ use ]]
  if (is(assay, "Assay5")) {
    object <- do.call(`$`, list(assay, layer))
  } else {
    stop("`x@object@assays[[ use ]]` is not 'Assay5'.")
  }
  rownames(object) <- gs(rownames(object), "\\.[0-9]*$", "")
  .job_kat(object = object)
}

setGeneric("asjob_kat", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_kat"))

setMethod("asjob_kat", signature = c(x = "job_seurat"),
  function(x, use = names(x@object@assays)[[1]], layer = "counts"){
    job_kat(x, use, layer)
  })

setMethod("step0", signature = c(x = "job_kat"),
  function(x){
    step_message("Prepare your data with function `job_kat`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_kat"),
  function(x, workers = 5, path = timeName("copykat"), test = FALSE)
  {
    step_message("Run copyKAT.")
    x$savepath <- path
    if (is.null(x$res_copykat)) {
      wd <- getwd()
      if (dir.exists(path)) {
        path <- paste0("copykat", gs(Sys.time(), " |:", "_"))
      }
      dir.create(path, FALSE)
      setwd(path)
      x@params$wd <- path
      if (isNamespaceLoaded("copykat") && !test) {
        pkgload::unload("copykat")
      } else {
        e(pkgload::load_all("~/copykat"))
      }
      full.anno <- dplyr::mutate(
        copykat::full.anno, hgnc_symbol = s(hgnc_symbol, "\\.[0-9]*", "")
      )
      full.anno <<- full.anno
      cyclegenes <<- copykat::cyclegenes
      DNA.hg20 <<- copykat::DNA.hg20
      message("According to `full.anno` to filter `object(x)`, and duplicated.")
      object(x) <- object(x)[rownames(object(x)) %in% full.anno$hgnc_symbol, ]
      object(x) <- object(x)[!duplicated(rownames(object(x))), ]
      res_copykat <- tryCatch(
        e(copykat::copykat(rawmat = object(x),
            genome = "hg20", n.cores = workers)), finally = setwd(wd)
      )
      try(rm(list = c("full.anno", "cyclegenes", "DNA.hg20"), envir = .GlobalEnv), silent = TRUE)
      x$res_copykat <- res_copykat
      x <- methodAdd(x, "R 包 `CopyKAT` 用于鉴定恶性细胞 {cite_show('DelineatingCopGaoR2021')}。`CopyKAT` 可以区分整倍体与非整倍体，其中非整倍体被认为是肿瘤细胞，而整倍体是正常细胞 {cite_show('CausesAndConsGordon2012')}。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_kat"),
  function(x){
    step_message("Results visualization.")
    p.copykat <- plot_heatmap.copyKAT(x$res_copykat)
    x@plots[[ 2 ]] <- namel(p.copykat)
    res_copykat <- dplyr::as_tibble(x$res_copykat$prediction)
    res_copykat <- dplyr::mutate(res_copykat,
      copykat.pred = gs(copykat.pred, ".*:([a-z]+):.*", "\\1"),
      copykat_cell = ifelse(copykat.pred == "aneuploid", "Cancer cell", "Normal cell")
    )
    x@tables[[ 2 ]] <- namel(res_copykat)
    return(x)
  })

setMethod("step3", signature = c(x = "job_kat"),
  function(x, name = x@params$wd) {
    step_message("save heatmap.")
    if (!is.null(x@params$res_copykat)) {
      x@params$res_copykat <- NULL
      x@object <- NULL
      file <- select_savefun(x@plots$step2$p.copykat)(x@plots$step2$p.copykat,
        name = name)
      newfile <- gs(file, "\\.pdf$", ".png")
      pdf_convert(file, filenames = newfile, dpi = 300)
      fig <- .file_fig(newfile)
      fig <- .set_lab(fig, sig(x), "malignant cells prediction heatmap")
      x@plots$step3$p.copykat <- fig
    }
    return(x)
  })

setMethod("merge", 
  signature = c(x = "job_seurat", y = "job_kat"),
  function(x, y, ref = "ChatGPT_cell", to = "copykat_cell") {
    if (!any(colnames(x@object@meta.data) == ref)) {
      stop("`ref` not found in meta.data.")
    }
    res <- y@tables$step2$res_copykat
    fun_cell <- function(x) gs(colnames(x@object), "_[0-9]+", "")
    res <- res$copykat_cell[ match(fun_cell(x), res$cell.names) ]
    res <- ifelse(is.na(res) | vapply(res, function(x) identical(x, "Normal cell"), logical(1)),
      x@object@meta.data[[ ref ]], res)
    x@object@meta.data[[ to ]] <- res
    p.map_cancer <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = FALSE, pt.size = .7,
        group.by = "copykat_cell", cols = color_set()
        ))
    p.map_cancer <- wrap(as_grob(p.map_cancer), 7, 4)
    p.map_cancer <- .set_lab(p.map_cancer, sig(x), "Cancer", "Cell type annotation")
    x@params$p.map_cancer <- p.map_cancer
    meth(x)$merge <- glue::glue("将 `CopyKAT` 的预测结果映射细胞注释中。")
    .append_heading(glue::glue("Seurat 癌细胞注释 ({x@sig})"), TRUE)
    return(x)
  })

setMethod("regroup", signature = c(x = "job_seurat", ref = "job_kat"),
  function(x, ref, k){
    if (ref@step < 2) {
      stop("ref@step < 2")
    }
    ka.tree <- readRDS(paste0(ref$savepath, "/_copykat_clustering_results.rds"))
    ka.tree$labels <- gs(ka.tree$labels, "\\.", "-")
    x <- regroup(x, ka.tree, k, TRUE)
    x$ka.tree <- ka.tree
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_kat"),
  function(x, ref, merge = x@params$group.by, cutree = NULL, pt.size = 1.5)
  {
    ref <- ref@tables$step2$res_copykat
    object(x)@meta.data$copykat_cell <-
      ref$copykat_cell[match(rownames(object(x)@meta.data), ref$cell.names)]
    if (!is.null(object(x)@meta.data[[ merge ]])) {
      anno_name <- paste0(gs(merge, "_cell$", "_"), "copykat")
      object(x)@meta.data[[ anno_name ]] <- as.factor(ifelse(
          object(x)@meta.data$copykat_cell == "Cancer cell",
          object(x)@meta.data$copykat_cell,
          as.character(object(x)@meta.data$scsa_cell))
      )
    }
    if (!is.null(cutree)) {
      if (is.logical(cutree)) {
        cutree <- seq(3, 36, 3)
      }
      ## visualize the marked cells in ...
      if (length(cutree) == 1) {
        obj <- regroup(x, kat, cutree)
        cells <- rownames(obj@object@meta.data)
        tree <- obj$ka.tree
        data <- tibble::tibble(name = tree$labels, value = order(tree$order), var = "Position")
        data <- dplyr::filter(data, name %in% dplyr::all_of(cells))
        data <- dplyr::mutate(data, group = obj@object@meta.data$regroup.hclust[match(name, cells)])
        levels <- names(sort(vapply(split(data$value, data$group), mean, double(1))))
        data <- dplyr::mutate(data, x = factor(group, levels = !!levels))
        p <- .map_boxplot2(data, FALSE, x = "x") +
          vis(obj, "regroup.hclust", pt.size)@data
        return(wrap(p, 12, 5))
      }
      p.list <- lapply(cutree,
        function(k) {
          obj <- regroup(x, kat, k)
          vis(obj, "regroup.hclust", pt.size)@data
        })
      p <- do.call(patchwork::wrap_plots, p.list)
      message("`k` for `cutree` set as ", paste0(cutree, collapse = ", "))
      return(wrap(p, 16, 12))
    }
    return(x)
  })

plot_heatmap.copyKAT <- function(copykat.obj) {
  require(copykat)
  ## the following code was wrriten by author of copykat, so terrible.
  my_palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
  col_breaks = c(seq(-1, -0.4, length = 50), seq(-0.4, -0.2, length = 150),
    seq(-0.2, 0.2, length = 600), seq(0.2, 0.4, length = 150),
    seq(0.4, 1, length = 50))
  pred.data <- data.frame(copykat.obj$prediction)
  pred.data <- pred.data[pred.data$copykat.pred != "not.defined",]
  CNA.data <- data.frame(copykat.obj$CNAmat)
  chr <- as.numeric(CNA.data$chrom) %% 2+1
  rbPal1 <- colorRampPalette(c('black', 'grey'))
  CHR <- rbPal1(2)[as.numeric(chr)]
  chr1 <- cbind(CHR, CHR)
  rbPal5 <- colorRampPalette(pal_group <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:2])
  com.preN <- pred.data$copykat.pred
  if (com.preN[1] == "diploid") {
    rbPal5 <- rev(rbPal5(2))
  } else {
    rbPal5 <- rbPal5(2)
  }
  pred <- rbPal5[as.numeric(factor(com.preN))]
  cells <- rbind(pred, pred)
  heatmap.3(t(CNA.data[, 4:ncol(CNA.data)]), dendrogram = "r",
    distfun = function(x) parallelDist::parDist(x, threads = 4, method = "euclidean"),
    hclustfun = function(x) hclust(x, method = "ward.D2"), ColSideColors = chr1,
    RowSideColors = cells, Colv = NA, Rowv = TRUE, notecol = "black", col = my_palette,
    breaks = col_breaks, key = FALSE, keysize = .5,
    density.info = "none", trace = "none", cexRow = 0.1, cexCol = 0.1, cex.main = 1,
    cex.lab = 0.1, symm = FALSE, symkey = FALSE, symbreaks = TRUE, cex = 1, cex.main = 4,
    margins = c(7, 5))
  data <- data.frame(x = c(-1, 1), y = c(-1, 1),
    Cell = factor(that <- c("diploid", "aneuploid"), levels = that))
  p <- ggplot(data) +
    geom_col(aes(x = x, y = y, fill = Cell)) +
    geom_col(aes(x = x, y = y, color = x)) +
    labs(fill = "", color = "") +
    scale_fill_manual(values = pal_group) +
    scale_color_gradient2(low = my_palette[1], high = tail(my_palette, 1),
      breaks = c(-1, 0, 1), labels = c("Loss", "", "Gain")) +
    theme(legend.position = "bottom")
  legend <- .get_legend(p)
  draw(ggather(legend, vp = viewport(unit(.5, "npc"), unit(1, "cm"))))
  plot <- recordPlot()
  return(plot)
  ## not used
  if (FALSE) {
    tumor.cells <- pred.data$cell.names[grepl("aneuploid", pred.data$copykat.pred)]
    tumor.mat <- CNA.data[, colnames(CNA.data) %in% make.names(tumor.cells)]
    hcc <- hclust(parallelDist::parDist(t(tumor.mat), threads = 4, method = "euclidean"), method = "ward.D2")
    hc.umap <- cutree(hcc, 2)
    rbPal6 <- colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4])
    subpop <- rbPal6(2)[as.numeric(factor(hc.umap))]
    cells <- rbind(subpop, subpop)
    heatmap.3(t(tumor.mat), dendrogram = "r", distfun = function(x)
      parallelDist::parDist(x, threads = 4, method = "euclidean"), hclustfun =
        function(x) hclust(x, method = "ward.D2"),
      ColSideColors = chr1, RowSideColors = cells, Colv = NA, Rowv = TRUE, 
      notecol = "black", col = my_palette, breaks = col_breaks, key = TRUE, keysize = 1,
      density.info = "none", trace = "none",
      cexRow = 0.1, cexCol = 0.1, cex.main = 1, cex.lab = 0.1, 
      symm = FALSE, symkey = FALSE, symbreaks = TRUE, cex = 1, cex.main = 4, margins = c(10, 10))
    legend("topright", c("c1", "c2"), pch = 15, 
      col = RColorBrewer::brewer.pal(n = 8, name = "Dark2")[3:4], cex = 0.9, bty = 'n')
  }
}
