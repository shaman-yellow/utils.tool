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

job_kat <- function(x, refs = "")
{
  if (!is(x, "dgCMatrix")) {
    stop('!is(x, "dgCMatrix").')
  }
  x <- .job_kat(object = x)
  x$refs <- refs
  return(x)
}

setGeneric("asjob_kat", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_kat"))

setMethod("asjob_kat", signature = c(x = "job_seurat"),
  function(x, refs = NULL, group.by = x$group.by, use = names(x@object@assays)[[1]], layer = "counts")
  {
    assay <- x@object@assays[[ use ]]
    if (is(assay, "Assay5")) {
      object <- do.call(`$`, list(assay, layer))
    } else {
      stop("`x@object@assays[[ use ]]` is not 'Assay5'.")
    }
    rownames(object) <- gname(rownames(object))
    if (is.null(refs)) {
      stop('is.null(refs).')
    }
    metadata <- object(x)@meta.data
    refs <- dplyr::filter(metadata, !!rlang::sym(group.by) %in% !!refs)
    if (!nrow(refs)) {
      stop('!nrow(refs).')
    }
    refs <- rownames(refs)
    if (is.null(refs)) {
      stop('is.null(refs), no rownames.')
    }
    x <- job_kat(object, refs)
    return(x)
  })

setMethod("step0", signature = c(x = "job_kat"),
  function(x){
    step_message("Prepare your data with function `job_kat`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_kat"),
  function(x, workers = 5, path = glue::glue("copykat_{x@sig}"), test = FALSE, ...)
  {
    step_message("Run copyKAT.")
    x$savepath <- path
    x$workers <- workers
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L, ...,
        {
          full.anno <- cyclegenes <- DNA.hg20 <- NULL
          if (packageVersion("Matrix") < "1.7-3") {
            x@object <- as.matrix(x@object)
          }
          x <- step1(x, workers = "{workers}", path = "{path}")
          x@object <- NULL
        }
      )
      return(x)
    }
    if (is.null(x$res_copykat)) {
      savedir <- getwd()
      if (dir.exists(path)) {
        path <- paste0("copykat", gs(Sys.time(), " |:", "_"))
      }
      dir.create(path, FALSE)
      setwd(path)
      x$savedir <- path
      if (isNamespaceLoaded("copykat")) {
        pkgload::unload("copykat")
      }
      full.anno <- dplyr::mutate(
        copykat::full.anno, hgnc_symbol = s(hgnc_symbol, "\\.[0-9]*", "")
      )
      full.anno <<- full.anno
      cyclegenes <<- copykat::cyclegenes
      DNA.hg20 <<- copykat::DNA.hg20
      message("According to `full.anno` to filter `object(x)`, and duplicated.")
      require(Matrix)
      object(x) <- object(x)[rownames(object(x)) %in% full.anno$hgnc_symbol, ]
      object(x) <- object(x)[!duplicated(rownames(object(x))), ]
      refs <- x$ref_cells
      if (is.null(refs)) {
        refs <- ""
      }
      res_copykat <- tryCatch(
        e(copykat::copykat(rawmat = object(x),
            genome = "hg20", n.cores = workers, norm.cell.names = refs)), finally = setwd(savedir)
      )
      # try(rm(list = c("full.anno", "cyclegenes", "DNA.hg20"), envir = .GlobalEnv), silent = TRUE)
      x$res_copykat <- res_copykat
      x <- methodAdd(x, "R 包 `CopyKAT` 用于鉴定恶性细胞 {cite_show('DelineatingCopGaoR2021')}。`CopyKAT` 可以区分整倍体与非整倍体，其中非整倍体被认为是肿瘤细胞，而整倍体是正常细胞 {cite_show('CausesAndConsGordon2012')}。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_kat"),
  function(x, workers = x$workers, inherits = TRUE, 
    local = FALSE, ignore = FALSE, ...)
  {
    step_message("Results visualization.")
    if (is.remote(x) && !local) {
      x <- run_job_remote(x, wait = 1L, inherit_last_result = inherits,
        ignore_local_cache = ignore, ignore_remote_cache = ignore, ...,
        {
          if ("{ignore}") {
            x@step <- as.integer(1L)
          }
          x <- step2(x, "{workers}")
        })
      x <- transmute_remote_figs(x)
    } else {
      p.copykat <- plot_heatmap_copyKAT(
        x$res_copykat, workers, name = x$savedir
      )
      p.copykat <- .set_lab(p.copykat, sig(x), "malignant cells prediction heatmap")
      p.copykat <- setLegend(p.copykat, "为 copyKAT 预测恶质细胞的热图。")
      x@plots[[ 2 ]] <- namel(p.copykat)
      res_copykat <- dplyr::as_tibble(x$res_copykat$prediction)
      res_copykat <- dplyr::mutate(res_copykat,
        copykat.pred = gs(copykat.pred, ".*:([a-z]+):.*", "\\1"),
        copykat_cell = ifelse(copykat.pred == "aneuploid", "Malignant cell", "Benign cell")
      )
      res_copykat <- set_lab_legend(
        res_copykat,
        glue::glue("{sig(x)} copyKAT prediction data"),
        glue::glue("为 copyKAT 预测恶质细胞结果附表。")
      )
      x@tables[[ 2 ]] <- namel(res_copykat)
    }
    return(x)
  })

transmute_remote_figs <- function(x, mode = c("plots", "params"), step = x@step)
{
  mode <- match.arg(mode)
  if (mode == "plots") {
    x@plots[[step]] <- lapply(x@plots[[step]], 
      function(plot) {
        if (is(plot, "fig")) {
          transmute(x, plot)
        } else {
          plot
        }
      })
  }
  return(x)
}

setMethod("transmute", signature = c(x = "job", ref = "fig"),
  function(x, ref){
    if (!is.remote(x)) {
      stop('!is.remote(x).')
    }
    if (length(ref) != 1) {
      stop('length(ref) != 1.')
    }
    file_remote <- ref@.Data
    if (file.exists(file_remote)) {
      stop('file.exists(file_remote), has been transmuted?')
    }
    file_local <- file.path(x$map_local, basename(file_remote))
    get_file_from_remote(file_remote, x$wd, file_local, x$remote)
    ref@.Data <- file_local
    return(ref)
  })

# setMethod("step3", signature = c(x = "job_kat"),
#   function(x, name = x$savedir, force = FALSE) {
#     step_message("save heatmap.")
#     if (!is.null(x@params$res_copykat) || force) {
#       x@object <- NULL
#       file <- select_savefun(x@plots$step2$p.copykat)(x@plots$step2$p.copykat,
#         name = name)
#       newfile <- s(file, "\\.pdf$", ".png")
#       pdf_convert(file, filenames = newfile, dpi = 300)
#       fig <- .file_fig(newfile)
#       fig <- .set_lab(fig, sig(x), "malignant cells prediction heatmap")
#       x@plots$step3$p.copykat <- fig
#     }
#     return(x)
#   })

setMethod("map", signature = c(x = "job_seurat", ref = "job_kat"),
  function(x, ref, from = "scsa_cell", to = "copykat_cell")
  {
    if (!any(colnames(x@object@meta.data) == from)) {
      stop("`from` not found in meta.data.")
    }
    # all used in ref
    res <- ref@tables$step2$res_copykat
    fun_cell <- function(x) colnames(x@object)
    res <- res$copykat_cell[ match(fun_cell(x), res$cell.names) ]
    res <- ifelse(is.na(res) | vapply(res, function(x) identical(x, "Benign cell"), logical(1)),
      as.character(x@object@meta.data[[ from ]]), res)
    x@object@meta.data[[ to ]] <- res
    x <- mutate(
      x, isCancer = ifelse(
        !!rlang::sym(to) == "Malignant cell", "Malignant cell", "Benign cell"
      )
    )
    Raws <- x@object@meta.data[[from]]
    palette <- .set_palette_in_ending(Raws, "Malignant cell")
    x@object@meta.data[[to]] <- factor(
      x@object@meta.data[[to]], levels = names(palette)
    )
    p.map_cancer <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = FALSE, pt.size = .7,
        group.by = to, cols = palette
        ))
    p.map_cancer <- wrap(as_grob(p.map_cancer), 7, 4)
    p.map_cancer <- .set_lab(p.map_cancer, sig(x), "Cancer", "Cell type annotation")
    p.props_cancer <- plot_cells_proportion(
      object(x)@meta.data, "isCancer", from, relative = FALSE
    )
    x$p.props_cancer <- set_lab_legend(
      p.props_cancer,
      glue::glue("{sig(x)} cancer cell proportions"),
      glue::glue("为 copyKAT 注释的恶质细胞在各个细胞类型中的占比。")
    )
    x@params$p.map_cancer <- p.map_cancer
    x <- snapAdd(x, "将 `CopyKAT` 的预测结果映射细胞注释中。", step = class(ref))
    x$.map_heading <- glue::glue("Seurat-copyKAT 癌细胞注释")
    return(x)
  })

.set_palette_in_ending <- function(Raws, extras = "Malignant cell", 
  extra_colors = "black", colors = color_set())
{
  if (is.character(Raws)) {
    allRawCells <- sort(unique(Raws))
  } else if (is.factor(Raws)) {
    allRawCells <- levels(Raws)
  } else {
    stop("Not either factor or character.")
  }
  palette <- nl(
    levels <- c(allRawCells, extras), 
    c(colors[seq_along(allRawCells)], extra_colors), FALSE
  )
  palette
}

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

setMethod("merge", signature = c(x = "job_seurat", y = "job_kat"),
  function(x, y, merge = x@params$group.by, cutree = NULL, pt.size = 1.5)
  {
    y <- y@tables$step2$res_copykat
    object(x)@meta.data$copykat_cell <-
      y$copykat_cell[match(rownames(object(x)@meta.data), y$cell.names)]
    if (!is.null(object(x)@meta.data[[ merge ]])) {
      anno_name <- paste0(gs(merge, "_cell$", "_"), "copykat")
      object(x)@meta.data[[ anno_name ]] <- as.factor(ifelse(
          object(x)@meta.data$copykat_cell == "Malignant cell",
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

plot_heatmap_copyKAT <- function(copykat.obj, workers = 4L,
  name = "copyKAT_heatmap")
{
  require(copykat)
  pred.data <- data.frame(copykat.obj$prediction)
  pred.data <- pred.data[pred.data$copykat.pred != "not.defined", ]
  CNA.data <- data.frame(copykat.obj$CNAmat)
  # set colors for values
  color_gainOrLoss <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 3, name = "RdBu")))(n = 999)
  breaks_gainOrLoss = c(seq(-1, -0.4, length = 50), seq(-0.4, -0.2, length = 150),
    seq(-0.2, 0.2, length = 600), seq(0.2, 0.4, length = 150),
    seq(0.4, 1, length = 50))
  # set colors for chroms
  color_chr <- c('black', 'grey')
  ColSideColors_chrs <- color_chr[ as.numeric(CNA.data$chrom) %% 2+1 ]
  ColSideColors_chrs <- cbind(ColSideColors_chrs, ColSideColors_chrs)
  colnames(ColSideColors_chrs) <- rep("CHR", 2)
  # set colors for cell groups, green and red (aneuploid, diploid)
  color_cellGroup <- c("#D95F02", "#1B9E77")
  RowSideColors_cells <- color_cellGroup[ as.integer(factor(pred.data$copykat.pred)) ]
  RowSideColors_cells <- rbind(RowSideColors_cells, RowSideColors_cells)
  rownames(RowSideColors_cells) <- rep("Cells", 2)
  distfun <- function(x) {
    file_cache <- paste0("dist_", digest::digest(CNA.data))
    if (file.exists(file_cache)) {
      message("File exists: ", file_cache)
      res <- readRDS(file_cache)
    } else {
      res <- parallelDist::parDist(x, threads = workers, method = "euclidean")
      saveRDS(res, file_cache)
    }
    res
  }
  # draw...
  pngDpi(pngfile <- paste0(name, ".png"))
  heatmap.3(t(CNA.data[, 4:ncol(CNA.data)]), dendrogram = "r",
    distfun = distfun, hclustfun = function(x) hclust(x, method = "ward.D2"),
    ColSideColors = ColSideColors_chrs, RowSideColors = RowSideColors_cells,
    Colv = NA, Rowv = TRUE, notecol = "black", col = color_gainOrLoss,
    breaks = breaks_gainOrLoss, key = FALSE, keysize = .5,
    density.info = "none", trace = "none", cexRow = 0.1, cexCol = 0.1, cex.main = 1,
    cex.lab = 0.1, symm = FALSE, symkey = FALSE, symbreaks = TRUE, cex = 1, cex.main = 4,
    margins = c(7, 5))
  data <- data.frame(
    x = c(-1, 1), y = c(-1, 1),
    Cell = as.factor(c("aneuploid", "diploid"))
  )
  p <- ggplot(data) +
    geom_col(aes(x = x, y = y, fill = Cell)) +
    geom_col(aes(x = x, y = y, color = x)) +
    labs(fill = "", color = "") +
    scale_fill_manual(values = color_cellGroup) +
    scale_color_gradient2(
      low = head(color_gainOrLoss, n = 1), high = tail(color_gainOrLoss, 1),
      breaks = c(-1, 0, 1), labels = c("Loss", "", "Gain")
    ) +
    theme(legend.position = "bottom")
  legend <- .get_legend(p)
  draw(ggather(legend, vp = viewport(unit(.5, "npc"), unit(1, "cm"))))
  dev.off()
  fig <- .file_fig(pngfile)
  return(fig)
}

setMethod("set_remote", signature = c(x = "job_kat"),
  function(x, wd = glue::glue("~/kat_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })



setupCopykat <- function() {
  require("copykat")
  data <- dplyr::mutate(
    copykat::full.anno, hgnc_symbol = gname(hgnc_symbol)
  )
  # full.anno <<- full.anno
  # cyclegenes <<- copykat::cyclegenes
  # DNA.hg20 <<- copykat::DNA.hg20
  replaceFunInPackage("full.anno", data, "copykat")
}
