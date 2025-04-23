# ==========================================================================
# workflow of infercnv
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_infercnv <- setClass("job_infercnv",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://bioconductor.org/packages/release/bioc/html/infercnv.html"),
    method = "Package inferCNV used for CNV anlysis and cancer cell prediction",
    tag = "scrna:cancer",
    analysis = "InferCNV 变异拷贝数分析"
    ))

setGeneric("asjob_infercnv", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_infercnv"))

setMethod("asjob_infercnv", signature = c(x = "job_seurat"),
  function(x, ref, groups = NULL, ..., recluster = TRUE, subset = "seurat_clusters",
    group.by = x$group.by, outdir = glue::glue("infercnv_{x@sig}"))
  {
    if (missing(ref)) {
      stop('missing(ref).')
    }
    metadata <- dplyr::select(
      as_tibble(object(x)@meta.data), rownames, !!rlang::sym(group.by), !!rlang::sym(subset)
    )
    if (any(!ref %in% metadata[[ group.by ]])) {
      stop('any(!ref %in% metadata[[ group.by ]]).')
    }
    if (is.null(groups)) {
      groups <- unique(as.character(metadata[[group.by]]))
      groups <- groups[ !groups %in% ref ]
    }
    allGroups <- unique(c(groups, ref))
    metadata <- dplyr::filter(metadata, !!rlang::sym(group.by) %in% allGroups)
    dir.create(outdir, FALSE)
    if (recluster && subset == "seurat_clusters" && !is.null(groups)) {
      hash <- digest::digest(
        list(x@sig, rownames(object(x)), groups)
      )
      file_submeta <- add_filename_suffix(
        file.path(outdir, "submetadata.rds"), hash
      )
      if (!file.exists(file_submeta)) {
        sub <- asjob_seurat_sub(x, !!rlang::sym(group.by) %in% groups)
        sub <- step1(sub)
        sub <- step2(sub)
        sub <- step3(sub, ...)
        metadata_sub <- as_tibble(object(sub)@meta.data)
        snap <- snap(metadata_sub) <- snap(sub)
        saveRDS(metadata_sub, file_submeta)
      } else {
        metadata_sub <- readRDS(file_submeta)
        snap <- snap(metadata_sub)
      }
      metadata <- map(
        metadata, "rownames", metadata_sub, "rownames", subset, col = subset
      )
    }
    if (!is.null(subset)) {
      metadata <- dplyr::mutate(
        metadata, group_subset = paste0(
          !!rlang::sym(group.by), "_", !!rlang::sym(subset)
          ),
        group_subset = ifelse(
          !!rlang::sym(group.by) %in% !!ref, 
          as.character(!!rlang::sym(group.by)), group_subset
        )
      )
      message(glue::glue("\n{showStrings(metadata$group_subset, trunc = FALSE)}"))
      metadata <- dplyr::select(metadata, rownames, group_subset)
      group.by <- "group_subset"
    }
    counts <- e(SeuratObject::LayerData(object(x), "count"))
    rownames(counts) <- gname(rownames(counts))
    counts <- counts[ !duplicated(rownames(counts)), ]
    message(glue::glue("Before Cells filter: {bind(dim(counts))}"))
    counts <- counts[, colnames(counts) %in% metadata$rownames]
    message(glue::glue("After Cells filter: {bind(dim(counts))}"))
    ranges <- get_gene_ranges(rownames(counts))
    counts <- counts[rownames(counts) %in% ranges$symbols, ]
    ranges <- ranges[match(rownames(counts), ranges$symbols), ]
    genes <- dplyr::select(
      ranges, symbols, seqnames, start, end
    )
    tmp.metadata <- file.path(outdir, "metadata.tsv")
    write_tsv(metadata, tmp.metadata, col.names = FALSE)
    tmp.genes <- file.path(outdir, "genes.tsv")
    write_tsv(genes, tmp.genes, col.names = FALSE)
    cell_groups <- unique(as.character(metadata[[ group.by ]]))
    obj.cnv <- e(infercnv::CreateInfercnvObject(counts,
        tmp.genes, tmp.metadata, ref))
    x <- .job_infercnv(object = obj.cnv)
    x$outdir <- outdir
    x$tmp.metadata <- tmp.metadata
    x$tmp.genes <- tmp.genes
    x <- methodAdd(x, "以 R 包 `infercnv` ({packageVersion('infercnv')}) 探索肿瘤单细胞 RNA 测序数据，以识别体细胞大规模染色体拷贝数的变异。")
    if (exists("snap") && !is.function(snap)) {
      message(glue::glue("Add seurat recluster snap."))
      x <- snapAdd(x, snap)
    }
    x <- snapAdd(x, "使用 `inferCNV` 识别肿瘤细胞的染色体拷贝数变异，选择 {bind(ref)} 为参考细胞 (正常细胞)，识别 {bind(groups)} 中的拷贝数变异。分析中，{bind(groups)} 以 {subset} 标识次级聚类。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Prepare your data with function `asjob_infercnv`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_infercnv"),
  function(x, workers = 4, cutoff = .1, hmm = FALSE, ...)
  {
    step_message("Run inferCNV.")
    message(crayon::yellow(glue::glue("In linux, too many cells may result in blank heatmap.")))
    # https://github.com/broadinstitute/infercnv/issues/362
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L, ...,
        {
          x <- step1(
            x, workers = "{workers}",
            cutoff = "{cutoff}", hmm = "{hmm}"
          )
        }
      )
      return(x)
    }
    if (dir.exists(x$outdir)) {
      unlink(x$outdir, TRUE, TRUE)
    }
    options(scipen = 100)
    object(x) <- e(infercnv::run(
        object(x),
        # cutoff = 1 works well for Smart-seq2
        # and cutoff = 0.1 works well for 10x Genomics
        num_threads = workers, cluster_by_groups = TRUE,
        cutoff = cutoff, denoise = TRUE, HMM = hmm, 
        analysis_mode = "subclusters",
        out_dir = x$outdir, save_rds = FALSE, save_final_rds = FALSE
        ))
    return(x)
  })

setMethod("step2", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Got results.")
    if (is.remote(x)) {
      dir <- file.path(x$map_local, x$outdir)
      if (!dir.exists(dir)) {
        get_file_from_remote(
          x$outdir, x$wd, x$map_local, recursive = TRUE
        )
      }
    } else {
      dir <- x$outdir
    }
    p.infer <- .file_fig(.cut_png_blank_space(file.path(dir, "infercnv.png")))
    p.infer <- set_lab_legend(
      p.infer,
      glue::glue("{x@sig} infercnv heatmap"),
      glue::glue("为 CNV 层次聚类热图。{.infercnv_heatmap_method()}")
    )
    x <- plotsAdd(x, p.infer = p.infer)
    return(x)
  })

setMethod("step3", signature = c(x = "job_infercnv"),
  function(x, k = 10, clear = TRUE)
  {
    step_message("Kmean and scoring CNV...")
    if (is.null(object(x))) {
      stop('is.null(object(x)).')
    }
    expr <- object(x)@expr.data
    obs <- unlist(
      object(x)@observation_grouped_cell_indices, use.names = FALSE
    )
    refs <- unlist(
      object(x)@reference_grouped_cell_indices, use.names = FALSE
    )
    clusters <- x$clusters <- kmeanMiniBatch(t(expr[, obs]), k)
    x$kmean_group_indices <- split(obs, paste0("C", clusters))
    x$data_group <- tibble::tibble(
      cells = colnames(expr)[obs], clusters = paste0("C", unname(clusters))
    )
    expr <- colMeans((expr - 1) ^ 2)
    expr_obs <- expr[obs]
    expr_refs <- expr[refs]
    data <- tibble::tibble(
      group = c(paste0("C", clusters), rep("Ref", length(expr_refs))),
      expr = c(expr_obs, expr_refs)
    )
    data <- dplyr::mutate(
      data, type = ifelse(group == "Ref", "Reference", "Observation"),
      type = factor(type, c("Reference", "Observation"))
    )
    dataPvalue <- test_unbalance_groups(
      data, "Ref", "group", "expr", alternative = "greater"
    )
    x$dataPvalue <- dataPvalue <- dplyr::mutate(
      dataPvalue, type = factor("Observation", levels(data$type))
    )
    x$dataPvalue <- set_lab_legend(
      x$dataPvalue,
      glue::glue("{x@sig} kmean cluster CNV score wilcox test"),
      glue::glue("`kmean` 聚类后各组细胞的 CNV score 与 Reference 比较的 wilcox.test 统计结果。")
    )
    p.violin <- ggplot(data) +
      geom_violin(aes(x = reorder(group, expr), y = expr, fill = group)) +
      geom_text(data = dataPvalue, aes(x = group, label = sig, y = max(data$expr))) +
      ggforce::facet_row(~ type, space = "free", scales = "free_x") +
      scale_fill_manual(values = color_set()) +
      theme_minimal()
    p.violin <- wrap(p.violin, 8, 3.5)
    p.violin <- set_lab_legend(
      p.violin, glue::glue("{x@sig} kmean cluster CNV score violin plot"), ""
    )
    x <- snapAdd(
      x, "以 `kmean` 聚类，将 Observation 细胞分为 {k} 组，对各组的 CNV score 与 Reference 细胞检验，显著差异的组为恶质细胞。"
    )
    Legend(p.violin) <- "为 `kmean` 聚类后，对各组细胞的 CNV score (细胞基因表达量的均方值 $\\text{RMS} = \\sqrt{\\frac{1}{n} \\sum_{i=1}^{n} x_i^2}$) 以 `wilcox.test` 与 Reference 组比较，单尾检验是否显著大于 Reference 组。P-value 结果以 * 标记。(&lt; 0.001, ***; &lt; 0.01, **; &lt; 0.05, *)"
    if (clear) {
      object(x) <- NULL
    }
    x <- plotsAdd(x, p.violin = p.violin)
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_infercnv"),
  function(x, ref, from = x$group.by, to = "infercnv_cell"){
    if (ref@step < 3) {
      stop('ref@step < 3.')
    }
    res <- ref$data_group
    group.sig <- dplyr::filter(ref$dataPvalue, pvalue < .05)$group
    res <- dplyr::mutate(
      res, isCancer = ifelse(clusters %in% group.sig, "Malignant cell", "Benign cell")
    )
    message(glue::glue("Is Malignant: {try_snap(res$isCancer == 'Malignant cell')}"))
    isCancer <- res$isCancer[ match(rownames(object(x)@meta.data), res$cells) ]
    message(glue::glue("Got annotation: {try_snap(!is.na(isCancer))}"))
    object(x)@meta.data[[ to ]] <- ifelse(
      is.na(isCancer) | isCancer == "Benign cell",
      as.character(object(x)@meta.data[[ from ]]),
      "Malignant cell"
    )
    palette <- .set_palette_in_ending(
      object(x)@meta.data[[ from ]], "Malignant cell"
    )
    object(x)@meta.data[[ to ]] <- factor(
      object(x)@meta.data[[ to ]], names(palette)
    )
    p.map_cancer <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = FALSE,
        group.by = to, cols = palette
        ))
    p.map_cancer <- wrap(as_grob(p.map_cancer), 7, 4)
    p.map_cancer <- .set_lab(p.map_cancer, sig(x), "Cancer", "Cell type annotation")
    x$p.map_cancer <- p.map_cancer
    x <- snapAdd(x, "将 `infercnv` 的预测结果映射细胞注释中。", step = class(ref))
    x$.map_heading <- glue::glue("Seurat-InferCNV 癌细胞注释")
    return(x)
  })

test_unbalance_groups <- function(data, ref, col.group = "group",
  col.value = "value", fun = function(...) wilcox.test(...)$p.value, ...)
{
  fun <- match.fun(fun)
  data <- split(data, data[[col.group]] == ref)
  if (any(vapply(data, nrow, integer(1)) <= 1)) {
    stop('any(vapply(data, nrow, integer(1)) <= 1).')
  }
  fun_test <- function(x, ...) {
    fun(x, data[["TRUE"]][[col.value]], ...)
  }
  res <- dplyr::reframe(
    data[["FALSE"]], pvalue = fun_test(
      !!rlang::sym(col.value), ...
      ),
    max = max(!!rlang::sym(col.value), na.rm = TRUE),
    min = min(!!rlang::sym(col.value), na.rm = TRUE),
    .by = !!rlang::sym(col.group)
  )
  res <- dplyr::mutate(res,
    sig = ifelse(
      pvalue < .001, "***", ifelse(
        pvalue < .01, "**", ifelse(pvalue < .05, "*", "")
      )
    )
    # sig = paste0(sig, " (", signif(pvalue, 4), ")")
  )
  res
}

kmeanMiniBatch <- function(mtx, k = 10, batch = 100, force = FALSE, ...) {
  hash <- digest::digest(list(mtx, k, batch, force))
  file_cache <- add_filename_suffix("kmean.rds", hash)
  if (file.exists(file_cache)) {
    message(glue::glue("Read from cache: {file_cache}"))
    clusters <- readRDS(file_cache)
  } else {
    if (nrow(mtx) < 1e6 && !force) {
      message("Use `kmeans` for small size data.")
      clusters <- kmeans(mtx, k)
    } else {
      if (!requireNamespace("ClusterR", quietly = TRUE)) {
        install.packages("ClusterR")
      }
      object <- e(ClusterR::MiniBatchKmeans(
        mtx, k, batch_size = batch, num_init = 5
      ))
      clusters <- e(ClusterR::predict_MBatchKMeans(mtx, object$centroids))
    }
    saveRDS(clusters, file_cache)
  }
  if (is(clusters, "kmeans")) {
    clusters <- clusters$cluster
  }
  clusters
}

.cut_png_blank_space <- function(file_png, threshold = .99, 
  max = .01)
{
  img <- png::readPNG(file_png)
  dims <- dim(img)
  is_rgb <- length(dims) >= 3 && (dims[3] %in% c(3, 4))
  if (is_rgb) {
    rgb <- img[,,1:3]
    is_white <- (rgb[,,1] >= threshold) & (rgb[,,2] >= threshold) & (rgb[,,3] >= threshold)
  } else {
    gray <- img[,,1]
    is_white <- (gray >= threshold)
  }
  white_rows <- apply(is_white, 1, all)
  n <- 1L
  status <- white_rows[1]
  group <- rep(0L, length(white_rows))
  for (i in seq_along(white_rows)) {
    if (status != white_rows[i]) {
      n <- n + 1L
      status <- !status
    }
    group[i] <- n
  }
  groups <- split(seq_along(white_rows), group)
  ratio <- lengths(groups) / length(white_rows)
  allStatus <- as.logical(seq_along(ratio) %% 2)
  if (!white_rows[1]) {
    allStatus <- !allStatus
  }
  areThat <- ratio > max & allStatus
  maxNum <- floor(length(white_rows) * max)
  groups[areThat] <- lapply(groups[areThat], 
    function(x) {
      head(x, maxNum)
    })
  keep <- unlist(groups)
  if (is_rgb) {
    cropped_img <- img[keep, , , drop = FALSE]
  } else {
    cropped_img <- img[keep, , drop = FALSE]
  }
  newfile <- add_filename_suffix(file_png, "crop")
  png::writePNG(cropped_img, newfile)
  return(newfile)
}

.infercnv_heatmap_method <- function() {
  "正常细胞的表达值绘制在顶部热图中，肿瘤细胞的表达值绘制在底部热图中，基因在染色体上从左到右排列。正常细胞表达数据实际上从肿瘤细胞表达数据中减去，从而得出差异，其中染色体区域扩增显示为红色块，染色体区域缺失显示为蓝色块。参考<https://github.com/broadinstitute/inferCNV/wiki/Interpreting-the-figure>"
}

setMethod("set_remote", signature = c(x = "job_infercnv"),
  function(x, wd = glue::glue("~/infercnv_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })


get_gene_ranges <- function(symbols, version = c("hg38", 
  "hg19"), gname = TRUE)
{
  raw <- symbols
  if (gname) {
    symbols <- gname(symbols)
  }
  version <- match.arg(version)
  name_fun <- glue::glue("TxDb.Hsapiens.UCSC.{version}.knownGene")
  if (!requireNamespace(name_fun, quietly = TRUE)) {
    BiocManager::install(name_fun)
  }
  db <- get_fun(name_fun, asNamespace(name_fun), "S4")
  ranges <- e(GenomicFeatures::genes(db))
  # if (!requireNamespace("EnsDb.Hsapiens.v86")) {
  #   BiocManager::install("EnsDb.Hsapiens.v86")
  #   # EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  # }
  entrez <- e(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "ALIAS",
    keys = symbols,
    column = c("ENTREZID")
  ))
  hasThats <- !is.na(entrez) & (entrez %in% ranges$gene_id)
  entrez <- entrez[ hasThats ]
  ranges <- ranges[ match(entrez, ranges$gene_id) ]
  ranges <- data.frame(ranges)
  ranges$symbols <- raw[ hasThats ]
  return(ranges)
}

