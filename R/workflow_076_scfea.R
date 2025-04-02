# ==========================================================================
# workflow of scfea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scfea <- setClass("job_scfea", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "scfea",
    info = c("https://github.com/changwn/scFEA/tree/master"),
    cite = "[@AGraphNeuralAlgham2021]",
    method = "The `scFEA` (python) was used to estimate cell-wise metabolic via single cell RNA-seq data",
    tag = "scrna:flux",
    analysis = "scFEA 单细胞数据的代谢通量预测"
    ))

setGeneric("asjob_scfea", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_scfea"))

setMethod("asjob_scfea", signature = c(x = "job_seurat"),
  function(x, cells = NULL, groups = NULL, group.by = x$group.by,
    org = c("human", "mouse"), assay = "RNA", dir = "scfea", ...)
  {
    org <- match.arg(org)
    message("Use org: ", org)
    if (org == "mouse") {
      moduleGene_file <- "module_gene_complete_mouse_m168.csv"
    } else if (org == "human") {
      moduleGene_file <- "module_gene_m168.csv"
    }
    if (nrow(object(x)@meta.data) > 50000) {
      stop('nrow(object(x)@meta.data) > 50000, too many cells for prediction.')
    }
    if (TRUE) {
      file <- file.path(pg("scfea_db"), moduleGene_file)
      fun_test <- function(file) {
        x <- readLines(file)[-1]
        x <- lapply(strsplit(x, ","), function(x) x[-1])
        unique(unlist(x))
      }
      refs <- fun_test(file)
    }
    objAssay <- object(x)@assays[[ assay ]]
    if (isMultiAssays(objAssay)) {
      message("Combined data from multiple layers of the assay.")
      data <- objAssay$counts
    } else {
      data <- objAssay@counts
      if (identical(dim(data), c(0L, 0L))) {
        data <- objAssay@data
      }
    }
    # actually used genes
    data <- data[ rownames(data) %in% refs, ]
    metadata <- meta(x)
    message("Subset by `cells`.")
    if (!is.null(groups)) {
      if (is.null(group.by)) {
        stop('is.null(group.by).')
      }
      cells <- metadata[[ group.by ]] %in% groups
      snapAdd_onExit("x", "将 `Seurat` ({bind(groups)} 细胞) 以 `scFEA` 预测代谢通量。")
    }
    if (!is.null(cells)) {
      data <- data[, cells]
      metadata <- metadata[cells, ]
    } else {
      snapAdd_onExit("x", "将 `Seurat` (所有细胞) 以 `scFEA` 预测代谢通量。")
    }
    message("Dim: ", paste0(dim(data), collapse = ", "))
    message("Max: ", max(apply(data, 2, max)))
    # test expression
    noAnyExpre <- which(apply(data, 2, sum) == 0)
    if (length(noAnyExpre)) {
      message(
        glue::glue("Some cells no any expression (of actually used genes),
          randomly mutate them (+1).\n{bind(colnames(data)[noAnyExpre])}")
      )
      set.seed(143254345L)
      whichGene <- sample(seq_len(nrow(data)), 10)
      data[whichGene, noAnyExpre] <- data[whichGene, noAnyExpre] + 1
    }
    if (dir.exists(dir)) {
      if (sureThat("Directory of `dir` exists, remove that?")) {
        unlink(dir, recursive = TRUE)
      }
    }
    dir.create(dir, FALSE)
    expr_file <- file.path(dir, "data.csv")
    data.table::fwrite(data.frame(data, check.names = FALSE), expr_file, row.names = TRUE)
    x <- job_scfea(expr_file, org = org)
    x@params <- c(x@params, list(metadata = metadata, from_seurat = TRUE))
    x <- methodAdd(x, "将 Seurat 的 `{assay}` Assay ('counts') 作为输入数据，以 `scFEA` 预测细胞的代谢通量 {cite_show('AGraphNeuralAlgham2021')}。参考 <https://github.com/changwn/scFEA/blob/master/scFEA_tutorial1.ipynb> 和 <https://github.com/changwn/scFEA/blob/master/scFEA_tutorial2.ipynb>。")
    return(x)
  })

isMultiAssays <- function(assay) {
  if (is(assay, "Assay5")) {
    return(length(assay@layers) > 1)
  } else {
    FALSE
  }
}

integrateLayers <- function(assay) {
  names <- names(assay@layers)
  objs <- sapply(names, simplify = FALSE,
    function(name) {
      do.call(`$`, list(assay, name))
    })
  commonGenes <- ins(lst = lapply(objs, function(x) rownames(x)))
  objs <- lapply(objs, function(x) x[ commonGenes, ])
  obj <- SeuratObject::CreateSeuratObject(counts = do.call(cbind, objs))
  obj@meta.data$.layer <- unlist(lapply(names, function(name) rep(name, ncol(objs[[ name ]]))))
  return(obj)
}

job_scfea <- function(expr_file, org = c("mouse", "human"), test = FALSE)
{
  org <- match.arg(org)
  x <- .job_scfea()
  x$wd <- dirname(expr_file)
  x$dir <- "."
  x$input_file <- basename(expr_file)
  if (org == "mouse") {
    x$moduleGene_file <- "module_gene_complete_mouse_m168.csv"
    x$stoichiometry_matrix <- "cmMat_complete_mouse_c70_m168.csv"
    x$cName_file <- "cName_complete_mouse_c70_m168.csv"
    x$module_annotation <- "Human_M168_information.symbols.csv"
  } else if (org == "human") {
    x$moduleGene_file <- "module_gene_m168.csv"
    x$stoichiometry_matrix <- "cmMat_c70_m168.csv"
    x$cName_file <- "cName_c70_m168.csv"
    x$module_annotation <- "Human_M168_information.symbols.csv"
  }
  if (test) {
    file <- paste0(pg("scfea_db"), "/", x$moduleGene_file)
    fun_test <- function(file) {
      x <- readLines(file)[-1]
      x <- lapply(strsplit(x, ","), function(x) x[-1])
      unique(unlist(x))
    }
    refs <- fun_test(file)
    data <- ftibble(expr_file)
    print(table(refs %in% data[[1]]))
  }
  return(x)
}

setMethod("step0", signature = c(x = "job_scfea"),
  function(x)
  {
    step_message("Prepare your data with function `asjob_scfea`.")
  })

setMethod("step1", signature = c(x = "job_scfea"),
  function(x)
  {
    step_message("Run prediction.")
    if (is.remote(x)) {
      dir <- x$map_local
    } else {
      dir <- x$wd
      dir.create(file.path(dir, "output"))
    }
    if (!file.exists(file.path(dir, "flux.csv"))) {
      rem_run(
        scriptPrefix(x), pg("scfeaPython", is.remote(x)), " ", pg("scfea", is.remote(x)),
        " --data_dir ", pg("scfea_db", is.remote(x)),
        " --input_dir ", x$dir,
        " --test_file ", x$input_file,
        " --res_dir ", x$dir,
        " --moduleGene_file ", x$moduleGene_file,
        " --stoichiometry_matrix ", x$stoichiometry_matrix,
        " --cName_file ", x$cName_file,
        " --output_flux_file ", "flux.csv",
        " --output_balance_file ", "balance.csv"
      )
      if (is.remote(x)) {
        testRem_file.exists(x, "flux.csv", x$wait)
      }
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_scfea"),
  function(x){
    step_message("Collate results.")
    if (is.remote(x)) {
      cdRun("scp -r ", x$remote, ":", x$wd, "/* ", x$map_local)
      dir <- x$map_local
    } else {
      dir <- x$wd
    }
    t.balance <- ftibble(file.path(dir, "balance.csv"))
    t.flux <- ftibble(file.path(dir, "flux.csv"))
    t.flux <- .set_lab(t.flux, sig(x), "metabolic flux matrix")
    t.flux <- setLegend(t.flux, "为细胞代谢通量矩阵 (各 `M_` 为代谢模块)。")
    fun <- function(file) {
      lst <- strsplit(readLines(file)[-1], ",")
      names(lst) <- lapply(lst, function(x) x[1])
      lst <- lapply(lst, function(x) x[-1])
      as_df.lst(lst, "module", "gene")
    }
    dir_db <- pg("scfea_db")
    t.moduleGenes <- fun(paste0(dir_db, "/", x$moduleGene_file))
    t.moduleGenes <- .set_lab(t.moduleGenes, sig(x), "genes related to metabolic flux")
    t.moduleGenes <- setLegend(t.moduleGenes, "代谢通量对应的模块的基因。")
    maybePlot <- list.files(dir, "loss_[0-9]+-[0-9]+\\.png", full.names = TRUE)
    if (length(maybePlot)) {
      p.loss <- .file_fig(maybePlot)
      p.loss <- .set_lab(p.loss, sig(x), "Convergency of the loss terms during training")
      p.loss <- setLegend(p.loss, "为 `scFEA` 训练过程的收敛曲线。")
    } else {
      p.loss <- NULL
    }
    t.anno <- ftibble(paste0(dir_db, "/", x$module_annotation))
    t.anno <- dplyr::mutate(t.anno, name = paste0(Compound_IN_name, " -> ", Compound_OUT_name))
    t.anno <- .set_lab(t.anno, sig(x), "annotation of metabolic flux")
    t.anno <- setLegend(t.anno, "各代谢模块的注释。")
    t.compounds <- tibble::tibble(
      name = c(t.anno$Compound_IN_name, t.anno$Compound_OUT_name),
      kegg = c(t.anno$Compound_IN_ID, t.anno$Compound_OUT_ID)
    )
    t.compounds <- dplyr::distinct(t.compounds)
    x@tables[[ 2 ]] <- namel(t.flux, t.anno, t.balance, t.moduleGenes, t.compounds)
    x@plots[[ 2 ]] <- namel(p.loss)
    return(x)
  })

setMethod("vis", signature = c(x = "job_scfea"),
  function(x, group.by, feature, groups = NULL, order = TRUE,
    cutoff = "tail", mag = 100, return_type = c("job", "plot"), name = NULL)
  {
    if (missing(group.by)) {
      stop('missing(group.by), group.by of Cell Type can not be missing.')
    }
    if (missing(feature)) {
      stop(
        'missing(feature), name of feature (e.g., Glucose -> G6P) can not be missing.'
      )
    }
    data <- x@tables$step2$t.flux
    anno <- x@tables$step2$t.anno
    data <- e(tidyr::pivot_longer(data, -V1, names_to = "Module", values_to = "Flux"))
    data <- map(data, "Module", anno, "V1", "name", col = "name")
    data <- dplyr::filter(data, name %in% !!feature)
    if (!nrow(data)) {
      stop('!nrow(data), can not match any feature.')
    }
    data <- map(data, "V1", x$metadata, "rownames", group.by, col = "group.by")
    if (order) {
      data <- dplyr::mutate(data, name = factor(name, levels = !!rev(feature)))
    }
    if (identical(cutoff, "tail")) {
      message("Cut the tailing low value.")
      density <- density(data$Flux)
      seqs <- seq(min(density$x), max(density$x), length.out = 10)
      valueCutOff <- max(density$y) / mag
      for (i in seqs) {
        if (all(density$y[ density$x > i ] < valueCutOff)) {
          data <- dplyr::filter(data, Flux < i)
          message(glue::glue("Cut the Flux by {i}"))
          break
        }
      }
    }
    if (!is.null(groups)) {
      data <- dplyr::filter(data, group.by %in% !!groups)
    }
    p.flux <- ggplot(data, aes(x = Flux, y = name, fill = group.by)) +
      e(ggridges::geom_density_ridges()) +
      facet_wrap(~ group.by, nrow = 1) +
      guides(fill = "none") +
      ggridges::theme_ridges()
    message("Automated scaling.")
    p.flux <- wrap(
      p.flux, min(20, length(unique(data$group.by)) * 3 + max(nchar(as.character(unique(data$name)))) * .1), 
      min(10, length(feature) * .4 + 1)
    )
    p.flux <- set_lab_legend(
      p.flux, glue::glue("{x@sig} {name} Cell flux ridge plot"),
      glue::glue("为 {name} 的代谢通量山脊图。")
    )
    return_type <- match.arg(return_type)
    message("Finished.")
    if (return_type == "job") {
      x$p.flux <- p.flux
      return(x)
    } else {
      return(p.flux)
    }
  })

setMethod("map", signature = c(x = "job_scfea", ref = "job_limma"),
  function(x, ref, group.by = ref$group.by, groups = NULL, feature = NULL, 
    which = NULL, ...)
  {
    if (ref@step < 2L) {
      stop('ref@step < 2L.')
    }
    if (x@step < 2L) {
      stop('x@step < 2L.')
    }
    if (is.null(group.by)) {
      stop('is.null(group.by).')
    }
    if (is.null(which)) {
      which <- names(ref@tables$step2$tops)
    } else if (is.numeric(which)) {
      which <- names(ref@tables$step2$tops)[which]
    }
    p.fluxs <- lapply(which, 
      function(which) {
        if (is.null(feature)) {
          top <- ref@tables$step2$tops[[which]]
          if (!nrow(top)) {
            return(NULL)
          }
          feature <- top$name
          groups <- .get_versus_cell(which, NULL)
        }
        message(glue::glue("The top has: {nrow(top)}."))
        vis(
          x, group.by = group.by, feature = feature, groups = groups, 
          return_type = "plot", name = .get_versus_cell(which)
        )
      })
    names(p.fluxs) <- .get_versus_cell(which)
    x$p.fluxs <- p.fluxs
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_scfea"),
  function(x, ref, dims = 1:5)
  {
    data_flux <- ref@tables$step2$t.flux
    if (nrow(data_flux) != nrow(x@object@meta.data)) {
      stop('nrow(data_flux) != nrow(x@object@meta.data).')
    }
    data_flux <- data.frame(
      data_flux[, -1], row.names = data_flux[[1]], check.names = FALSE
    )
    data_flux <- t(as.matrix(data_flux))
    object(x)[["FLUX"]] <- e(SeuratObject::CreateAssayObject(counts = data_flux))
    oAssay <- SeuratObject::DefaultAssay(object(x))
    SeuratObject::DefaultAssay(object(x)) <- 'FLUX'
    object(x) <- e(
      Seurat::FindVariableFeatures(object(x), selection.method = "vst", nfeatures = 2000)
    )
    object(x) <- e(Seurat::ScaleData(object(x), features = rownames(object(x)), assay = 'FLUX'))
    object(x) <- e(Seurat::RunPCA(
      object(x), features = SeuratObject::VariableFeatures(object = object(x)), 
      npcs = 10, reduction.name = 'pca.flux', assay = "FLUX"
    ))
    p.flux_pca_rank <- pretty_elbowplot(
      e(Seurat::ElbowPlot(object(x), reduction = "pca.flux"))
    )
    p.flux_pca_rank <- wrap(p.flux_pca_rank, 4, 4)
    p.flux_pca_rank <- .set_lab(p.flux_pca_rank, sig(x), "Ranking of principle components of metabolic flux")
    x$p.flux_pca_rank <- p.flux_pca_rank
    object(x) <- e(Seurat::FindNeighbors(
      object(x), dims = dims, assay = "FLUX", reduction = "pca.flux"
    ))
    x <- mutate(
      x, seurat_clusters_RNA = seurat_clusters, .after = seurat_clusters
    )
    object(x) <- e(Seurat::FindClusters(object(x), resolution = 0.5))
    object(x) <- e(Seurat::RunUMAP(
      object(x), dims = dims, assay = 'FLUX', reduction.name = "umap.flux"
    ))
    p.map_flux <- vis(x, reduction = "umap.flux")
    p.map_flux <- .set_lab(p.map_flux, sig(x), "cells metabolic flux")
    x$p.map_flux <- setLegend(p.map_flux, "为细胞代谢通量 (`scFEA` 预测，输入 `Seurat`) 的 UMAP 聚类。")
    SeuratObject::DefaultAssay(object(x)) <- oAssay
    x <- mutate(
      x, seurat_clusters_FLUX = seurat_clusters,
      seurat_clusters = seurat_clusters_RNA,
      .after = seurat_clusters
    )
    x <- snapAdd(
      x, "将 `scFEA` 的代谢通量预测，输入 `Seurat` 数据对象中，按标准工作流分析，以细胞群聚类。",
      step = "job_scfea"
    )
    x$.map_snap <- "job_scfea"
    return(x)
  })

setMethod("cal_corp", signature = c(x = "job_limma", y = "job_seurat"),
  function(x, y, from = tops$name, to = tops$gene, names = NULL,
    tops = y@tables$step2$tops[[ 1 ]])
  {
    if (is.null(x$from_scfea)) {
      stop("The 'job_limma' is not converted from job_scfea.")
    }
    if (is(to, "list")) {
      groupTo <- to
      to <- unique(unlist(to))
    }
    stop("...")
  })

setMethod("regroup", signature = c(x = "job_scfea", "job_seurat"),
  function(x, ref){
    if (is.null(x$metadata$rownames)) {
      stop('is.null(x$metadata$rownames), no cell names?')
    }
    newMeta <- as_tibble(ref@object@meta.data)
    x$metadata <- newMeta[match(x$metadata$rownames, newMeta$rownames), ]
    return(x)
  })

setMethod("asjob_limma", signature = c(x = "job_scfea"),
  function(x, metadata, group, scale_sample = TRUE, scale_var = FALSE, ...)
  {
    if (missing(metadata) && !is.null(x$from_seurat)) {
      message("Use metadata from `x$metadata`")
      metadata <- x$metadata
      if (missing(group)) {
        group <- "group"
        guess_group <- tail(colnames(metadata), n = 1)
        message(
          glue::glue("Guess group based on '{guess_group}', trunc the numbering in 'orig.ident'.")
        )
        metadata <- dplyr::mutate(
          metadata, sample = rownames, group = paste0(
            !!rlang::sym(guess_group), ".", gs(orig.ident, "[0-9]*$", "")
          ), group = gs(make.names(group), "[.]+", "_")
        )
        message("Use group:\n", showStrings(metadata[[ group ]]))
      }
    }
    counts <- data.frame(x@tables$step2$t.flux)
    rownames(counts) <- counts[[1]]
    counts <- counts[, -1]
    if (scale_sample) {
      counts <- scale(counts, ...)
    }
    counts <- t(counts)
    if (scale_var) {
      counts <- scale(counts, ...)
    }
    counts <- as_tibble(counts, .name_repair = "minimal")
    genes <- x@tables$step2$t.anno
    annoExtra <- reframe_col(as_tibble(x@tables$step2$t.moduleGenes), "gene", list)
    genes <- map(genes, colnames(genes)[1], annoExtra, "module", "gene", col = "gene")
    if (!any(metadata[[1]] %in% colnames(counts))) {
      stop("Use first column of `metadata` as sample name (cell name), but not match any.")
    } else {
      message("Use first column of `metadata` as sample name (cell name).")
      metadata <- dplyr::filter(metadata, !!rlang::sym(colnames(metadata)[1]) %in% !!colnames(counts))
      if (!all(c("sample", "group") %in% colnames(metadata))) {
        if (any(colnames(metadata) == "group")) {
          metadata <- metadata[, colnames(metadata) != "group"]
        }
        metadata <- dplyr::rename(
          metadata, sample = !!rlang::sym(colnames(metadata)[1]),
          group = !!rlang::sym(group)
        )
      }
      metadata <- dplyr::relocate(metadata, sample, group)
    }
    cli::cli_alert_info("job_limma_normed")
    compounds_annotation <- x@tables$step2$t.compounds
    x <- job_limma_normed(counts, metadata, genes = genes)
    x <- snapAdd(x, "以 `limma` 的线形分析策略，对细胞的代谢通量差异分析。")
    x@analysis <- "Limma 代谢通量差异分析"
    x$from_scfea <- TRUE
    x$compounds_annotation <- compounds_annotation
    x$group.by <- group
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_scfea"),
  function(x, wd = "scfea")
  {
    local_wd <- x$wd
    x$wd <- "."
    if (rem_file.exists(wd)) {
      isThat <- sureThat("Dir exists, remove all files ?")
      if (isThat) {
        cdRun("ssh ", x$remote, " 'rm -r ", wd, "'")
      }
    }
    rem_dir.create(wd)
    x$wd <- wd
    rem_dir.create("output")
    cdRun("scp ", file.path(local_wd, "data.csv"), " ", x$remote, ":", wd)
    x$map_local <- local_wd
    x$dir <- "."
    return(x)
  })

setMethod("intersect", signature = c(x = "job_limma", y = "character"),
  function(x, y, key = 1, use = "kegg", split = if (use == "kegg") "\\+" else NULL,
    x.name = "Diff_flux", y.name = "Diff_meta")
  {
    if (is.null(x$from_scfea)) {
      stop("The value got TRUE: `is.null(x$from_scfea)`")
    }
    if (is.null(x@tables$step2[[ key ]])) {
      stop("The value got TRUE: `is.null(x@tables$step2[[ key ]])`")
    }
    data <- x@tables$step2[[ key ]]
    if (!is.null(split) && use == "kegg") {
      data <- dplyr::mutate(data,
        kegg_split = lapply(kegg,
          function(x) {
            unlist(strsplit(x, split))
          }))
      alls <- unlist(data$kegg_split)
      p.venn <- new_venn(lst = nl(c(x.name, y.name), list(alls, y)))
      data <- dplyr::filter(data,
        vapply(kegg_split,
          function(x) {
            any(x %in% y)
          }, logical(1))
      )
      data <- .set_lab(data, sig(x), "intersection-flux-data")
      .append_heading("交集：差异代谢物+单细胞差异代谢通量相关代谢物")
      return(namel(data, p.venn))
    } else {
      stop("...")
    }
  })

