# ==========================================================================
# workflow of monocle
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_monocle <- setClass("job_monocle", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/"),
    cite = "[@ReversedGraphQiuX2017; @TheDynamicsAnTrapne2014]",
    method = "R package `Monocle3` used for cell pseudotime analysis",
    tag = "scrna:pseudo",
    analysis = "Monocle3 拟时分析"
    ))

setGeneric("do_monocle", 
  function(x, ref, ...) standardGeneric("do_monocle"))

setMethod("do_monocle", signature = c(x = "job_seurat", ref = "character"),
  function(x, ref, dims = 1:15, resolution = 1.2, group.by = x@params$group.by)
  {
    x <- getsub(x, cells = grp(x@object@meta.data[[group.by]], ref))
    snapAdd_onExit("x", "从 `Seurat` 数据对象 {group.by} 中提取 {ref} 类型的细胞，对其重新聚类分析。")
    x@step <- 2L
    sr_sub <- step3(x, dims, resolution)
    methodAdd_onExit("x", "从 `Seurat` 数据对象 {group.by} 中提取 {ref} 类型的细胞，对其重新聚类分析。")
    methodAdd_onExit("x", sr_sub@meth[[ "step3" ]])
    methodAdd_onExit(
      "x", "使用 `SeuratWrappers` (`SeuratWrappers::as.cell_data_set`, 参考 <http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html>) 将 `Seurat` 转化为 `Monocle3` 的 `cell_data_set` 数据。该转化将继承 `Seurat` 前期分析的 PCA、UMAP 等聚类结果，用于 `Monocle3` 的拟时分析，使前后分析一致。"
    )
    snapAdd_onExit("x", "将 `Seurat` 数据对象转化为 `Monocle3` 数据对象 (详见方法章节)。")
    x <- asjob_monocle(sr_sub, "seurat_clusters")
    x$sr_sub <- sr_sub
    return(x)
  })

setGeneric("asjob_monocle", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_monocle"))

setMethod("asjob_monocle", signature = c(x = "job_seurat"),
  function(x, group.by = x@params$group.by, ..., use = names(x@object@assays)[[1]]){
    step_message("
      Other parameters would be passed to `SeuratWrappers::as.cell_data_set`.
      <http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html>
      "
    )
    if (x@step < 3) {
      stop("x@step < 3. At least step 3 was needed with `Seurat::FindClusters`.")
    }
    if (is.null(group.by))
      stop("is.null(group.by) == TRUE")
    if (FALSE) {
      # https://github.com/cole-trapnell-lab/monocle3/issues/438
      actAssay <- SeuratObject::DefaultAssay(object(x))
      if (actAssay != "RNA") {
        names <- names(object(x)@assays)
        names(object(x)@assays)[match(c("RNA", actAssay), names)] <- c(actAssay, "RNA")
        object(x)@active.assay <- 'RNA'
        message(crayon::red("As issues in <https://github.com/cole-trapnell-lab/monocle3/issues/438>,",
            "herein, the assays would be changed name before obtain."
            ))
      }
    }
    if (TRUE) {
      # https://github.com/cole-trapnell-lab/monocle3/issues/438
      message(crayon::red("As issues in <https://github.com/cole-trapnell-lab/monocle3/issues/438>,",
          "herein, set default assays as 'RNA' before obtain.
          NOTE: this results in all genes but not genes of integrated stored in monocle object;
          however, the 'UMAP' data were inherits from 'Seurat' object, so the plots were consistent.
          "
          ))
      object(x)@active.assay <- use
    }
    palette <- x$palette
    object <- e(SeuratWrappers::as.cell_data_set(object(x), group.by = group.by, ...))
    mn <- .job_monocle(object = object)
    object(mn) <- e(monocle3::estimate_size_factors(object(mn)))
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- object(x)@reductions[["pca"]]@feature.loadings
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- object(x)@reductions$pca@stdev
    mn@params$group.by <- group.by
    mn@params$palette <- palette
    mn 
  })

setMethod("step0", signature = c(x = "job_monocle"),
  function(x){
    step_message("Prepare your data with methods `asjob_monocle`.
      A processed 'Seurat' object is needed. (job_seurat at least in step 3).
      "
    )
  })

setMethod("step1", signature = c(x = "job_monocle"),
  function(x, groups = x@params$group.by, pt.size = 1.5, pre = FALSE, norm_method = "none"){
    step_message("Constructing single-cell trajectories.
      red{{`groups`}} would passed to `monocle3::plot_cells` for
      annotation in plot. Mutilple group could be given.
      "
    )
    x$pt.size <- pt.size
    palette <- x$palette
    if (is.null(palette)) {
      palette <- color_set()
    }
    if (!all(groups %in% colnames(object(x)@colData)))
      stop("Some of `groups` not found in `colData` of `object(x)`")
    if (pre) {
      object(x) <- e(monocle3::preprocess_cds(object(x), norm_method = norm_method))
    }
    object(x) <- e(monocle3::cluster_cells(object(x)))
    x <- methodAdd(
      x, "以 `monocle3::cluster_cells` 计算细胞群的 'clusters' 和 'partitions'。"
    )
    object(x) <- e(monocle3::learn_graph(object(x)))
    x <- methodAdd(x, "以 `monocle3::learn_graph` 从高维空间 (high-dimensional space) 中构建 'trajectory'。")
    x <- snapAdd(x, "构建轨迹图 (Trajectory)。")
    p.traj <- e(sapply(groups, simplify = FALSE,
        function(group) {
          p <- monocle3::plot_cells(object(x), color_cells_by = group,
            label_cell_groups = TRUE, label_branch_points = TRUE,
            group_label_size = 4, graph_label_size = 2,
            cell_size = x$pt.size, cell_stroke = 0, alpha = .7
          )
          p <- wrap(p + scale_color_manual(values = palette))
          p <- .set_lab(p, sig(x), "trajectories of", group)
        }))
    p.prin <- monocle3::plot_cells(
      object(x), color_cells_by = groups[1],
      label_cell_groups = FALSE, label_principal_points = TRUE,
      graph_label_size = 3, cell_size = x$pt.size, cell_stroke = 0,
      alpha = .7
    )
    p.prin <- p.prin + scale_color_manual(values = palette)
    p.prin <- wrap(p.prin, 10, 7)
    p.prin <- .set_lab(p.prin, sig(x), "principal points")
    p.prin <- setLegend(p.prin, "为拟时轨迹与 principal point 示意。")
    x@plots[[ 1 ]] <- namel(p.traj, p.prin)
    return(x)
  })

setMethod("step2", signature = c(x = "job_monocle"),
  function(x, roots){
    step_message("
      red{{`roots`}} would passed to `monocle3::order_cells` for setting
      as roots of `pseudotime`. If `roots` is character with names, processed
      by red{{`get_earliest_principal_nodes`}} and then passed to parameter `root_pr_nodes`;
      else, directly passed to param `root_pr_nodes`.
      "
    )
    x <- methodAdd(x, "选择 {bind(roots)} (principle points) 为拟时起点，以 `monocle3::order_cells` 将细胞排序，随后构建细胞拟时变化图。")
    x <- snapAdd(x, "选择 {bind(roots)} (principle points) 为拟时起点。")
    if (!is.null(names(roots))) {
      if (length(roots) > 1)
        stop("With name of `roots`, but length(roots) > 1")
      roots <- get_earliest_principal_nodes(object(x), names(roots), unname(roots))
    }
    object(x) <- e(monocle3::order_cells(object(x), root_pr_nodes = roots))
    x <- snapAdd(x, "构建细胞的拟时变化图 (Pseudotime)。")
    p.pseu <- e(monocle3::plot_cells(
        object(x), color_cells_by = "pseudotime",
        label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = FALSE, graph_label_size = 3,
        cell_size = x$pt.size, cell_stroke = 0, alpha = .7
        ))
    p.pseu <- wrap(p.pseu, 6, 5)
    p.pseu <- setLegend(p.pseu, "为细胞拟时图。")
    p.pseu <- .set_lab(p.pseu, sig(x), "pseudotime")
    x@plots[[ 2 ]] <- namel(p.pseu)
    return(x)
  })

setMethod("step3", signature = c(x = "job_monocle"),
  function(x, formula_string = NULL, group.by = NULL, cores = 4){
    step_message("This step do:
      1. Regression analysis;
      2. Graph-autocorrelation analysis;
      3. Finding modules of co-regulated genes (Significant genes).
      4. Significant genes in co-regulated modules.
      "
    )
    if (!is.null(group.by)) {
      x$group.by <- group.by
    }
    if (length(unique(object(x)@colData[[ x$group.by ]])) <= 2) {
      stop("Too few clusters for gene modules finding...")
    }
    if (length(x@tables) < 3)
      x@tables[[ 3 ]] <- list()
    if (!is.null(formula_string)) {
      if (!is.character(formula_string)) {
        stop("is.character(formula_string) == FALSE")
      }
      if (is.null(x@tables[[ 3 ]]$fit_coefs)) {
        if (FALSE) {
          if (isNamespaceLoaded("Seurat")) {
            ## Oh, this is useless, still error
            message(crayon::red("As <https://github.com/cole-trapnell-lab/monocle3/issues/522>,",
                "herein, unload the package `Seurat`"))
            pkgload::unload("Seurat")
          }
        }
        fits <- e(monocle3::fit_models(object(x), formula_string, cores = cores, verbose = TRUE))
        fit_coefs <- e(monocle3::coefficient_table(fits))
        fit_coefs.sig <- dplyr::filter(fit_coefs, term != "(Intercept)", q_value < .05)
        fit_coefs.sig <- dplyr::arrange(fit_coefs.sig, q_value)
        fit_goodness <- e(monocle3::evaluate_fits(fits))
      } else {
        fit_coefs <- x@tables[[ 3 ]]$fit_coefs
        fit_coefs.sig <- x@tables[[ 3 ]]$fit_coefs.sig
        fit_goodness <- x@tables[[ 3 ]]$fit_goodness
      }
    }
    x <- methodAdd(x, "以 `monocle3::graph_test` 寻找单细胞拟时轨迹中差异表达的基因。")
    x <- snapAdd(x, "寻找单细胞拟时轨迹 (`monocle3::graph_test`) 中差异表达的基因。")
    if (is.null(x@tables[[ 3 ]]$graph_test)) {
      graph_test <- e(monocle3::graph_test(object(x), neighbor_graph = "knn", cores = cores))
      graph_test <- as_tibble(graph_test)
      graph_test.sig <- dplyr::filter(graph_test, q_value < .05)
      graph_test.sig <- dplyr::arrange(graph_test.sig, dplyr::desc(morans_I), q_value)
      graph_test.sig <- dplyr::rename(graph_test.sig, gene_id = 1)
      graph_test.sig <- .set_lab(graph_test.sig, sig(x), "Graph Test Significant genes")
      graph_test.sig <- setLegend(graph_test.sig, "为 Moran’s I test (Graph Test) 筛选的差异表达基因 (Q-value cutoff: 0.05)。")
    } else {
      graph_test <- x@tables[[ 3 ]]$graph_test
      graph_test.sig <- x@tables[[ 3 ]]$graph_test.sig
    }
    if (!is.null(formula_string)) {
      cross.sig <- graph_test.sig$gene_id[ graph_test.sig$gene_id %in% fit_coefs.sig$gene_id ]
      gene_sigs <- list(
        fit_coefs.sig = fit_coefs.sig$gene_id,
        graph_test.sig = graph_test.sig$gene_id,
        cross.sig = cross.sig
      )
    } else {
      gene_sigs <- list(graph_test.sig = graph_test.sig$gene_id)
    }
    if (!is.factor(object(x)@colData[[ x$group.by ]])) {
      object(x)@colData[[ x$group.by ]] %<>% as.factor()
    }
    cell_group <- tibble::tibble(
      cell = row.names(SummarizedExperiment::colData(object(x))), 
      group = SummarizedExperiment::colData(object(x))[[ x@params$group.by ]]
    )
    x@params$cell_group <- cell_group
    if (is.null(x@tables[[ 3 ]]$gene_module)) {
      gene_module <- try(cal_modules.cds(object(x), gene_sigs, cell_group), TRUE)
      gene_module_heatdata <- lapply(gene_module,
        function(lst) {
          if (!is.null(lst$aggregate)) {
            if (nrow(lst$aggregate) > 0)
              callheatmap(new_heatdata(as_data_long(lst$aggregate, , "module", "cell")))
          }
        })
      gene_module_heatdata$graph_test.sig <- .set_lab(gene_module_heatdata$graph_test.sig,
        sig(x), "gene module heatmap")
      x@plots[[ 3 ]] <- namel(gene_module_heatdata)
    } else {
      gene_module <- x@tables[[ 3 ]]$gene_module
    }
    if (!is.null(formula_string)) {
      x@tables[[ 3 ]] <- namel(fit_coefs, fit_coefs.sig, fit_goodness,
        graph_test, graph_test.sig, gene_module)
    } else {
      x@tables[[ 3 ]] <- namel(graph_test, graph_test.sig, gene_module)
    }
    x$cellClass_tree.gene_module <- try(hclust(dist(t(gene_module$graph_test.sig$aggregate))), silent = TRUE)
    return(x)
  })

setMethod("step4", signature = c(x = "job_monocle"),
  function(x, groups = ids(x), genes, group.by = NULL, cutoff = .5,
    cutoff.den = 1, group.den = "orig.ident")
  {
    step_message("Plot genes (in branch) that change as a function of pseudotime.
      red{{`groups`}} and red{{`genes`}} were used to subset the `object(x)`."
    )
    cds <- object(x)
    fun_sub <- selectMethod("[", class(cds))
    gene.groups <- grouping_vec2list(genes, 10, TRUE)
    pblapply <- pbapply::pblapply
    colData <- SummarizedExperiment::colData
    if (is.null(group.by)) {
      group.by <- x@params$group.by
    }
    genes_in_pseudotime <- e(pblapply(gene.groups,
        function(genes) {
          if (!is.null(groups)) {
            cds <- fun_sub(cds,
              rownames(cds) %in% genes,
              colData(cds)[[ group.by ]] %in% groups
            )
          } else {
            cds <- fun_sub(cds, rownames(cds) %in% genes, )
          }
          p <- monocle3::plot_genes_in_pseudotime(cds,
            label_by_short_name = FALSE,
            color_cells_by = group.by,
            min_expr = cutoff
          )
          data <- p$data
          p <- wrap(p, 6, length(genes) * 1.6)
          p <- setLegend(p, "为 {bind(genes)} 在拟时 (Pseudotime) 中的表达变化。")
          namel(p, data)
        }))
    names(genes_in_pseudotime) <- paste0("pseudo", seq_along(genes_in_pseudotime))
    dat_genes_in_pseudotime <- lapply(genes_in_pseudotime, function(x) x$data)
    theme <- rstyle("theme")
    plot_density <- lapply(dat_genes_in_pseudotime,
      function(data) {
        p <- ggplot(dplyr::filter(data, expression > cutoff.den)) +
          ggplot2::geom_density(aes(pseudotime, color = !!rlang::sym(group.den))) +
          labs(y = paste0("density (expression > ", cutoff.den, ")")) +
          theme
        wrap(p, 7, 2.5)
      })
    genes_in_pseudotime <- lapply(genes_in_pseudotime, function(x) x$p)
    genes_in_pseudotime <- .set_lab(
      genes_in_pseudotime, sig(x), 
      paste0("Set", seq_along(genes_in_pseudotime)), "genes in pseudotime"
    )
    x@plots[[ 4 ]] <- namel(genes_in_pseudotime, plot_density)
    x@tables[[ 4 ]] <- namel(dat_genes_in_pseudotime)
    return(x)
  })


setMethod("regroup", signature = c(x = "job_seurat", ref = "hclust"),
  function(x, ref, k, by.name = FALSE, rename = NULL){
    ref <- cutree(ref, k)
    x$cutree <- ref
    regroup(x, ref, k, by.name, rename)
  })

setMethod("regroup", signature = c(x = "job_seurat", ref = "integer"),
  function(x, ref, k, by.name = FALSE, rename = NULL){
    if (!is.null(rename)) {
      if (is.character(rename)) {
        ref[] <- paste0(rename, "_", ref[])
      } else {
        ref[] <- dplyr::recode(ref, !!!(rename),
          .default = as.character(unname(ref)))
      }
    }
    idents <- SeuratObject::Idents(object(x))
    if (by.name) {
      re.idents <- nl(names(idents), dplyr::recode(names(idents), !!!ref), FALSE)
    } else {
      re.idents <- nl(names(idents), dplyr::recode(unname(idents), !!!ref), FALSE)
    }
    re.idents <- factor(re.idents, levels = sort(as.character(unique(re.idents))))
    e(SeuratObject::Idents(object(x)) <- re.idents)
    object(x)@meta.data$regroup.hclust <- unname(re.idents)
    return(x)
  })

setMethod("ids", signature = c(x = "job_monocle"),
  function(x, id = x@params$group.by, unique = TRUE){
    ids <- SummarizedExperiment::colData(object(x))[[ id ]]
    if (unique)
      unique(ids)
    else
      ids
  })

setMethod("add_anno", signature = c(x = "job_monocle"),
  function(x, branches = NULL)
  {
    metaPrin <- igraph::V(monocle3::principal_graph(object(x))[[ "UMAP" ]])
    cellPrin <- object(x)@principal_graph_aux[[ "UMAP" ]]$pr_graph_cell_proj_closest_vertex[, ]
    cellPrin[] <- metaPrin$name[ cellPrin ]
    object(x)@colData[[ "principal_node" ]] <-
      unname(cellPrin)[ match(rownames(object(x)@colData), names(cellPrin)) ]
    object(x)@colData[[ "pseudotime" ]] <- monocle3::pseudotime(object(x))
    if (!is.null(branches)) {
      branches <- get_branches.mn(x, branches)
      branches <- as_df.lst(branches)
      branches <- dplyr::distinct(branches, name, .keep_all = TRUE)
      branches <- dplyr::arrange(branches, type)
      branches <- dplyr::mutate(branches,
        branch = paste0("time_", unlist(lapply(table(type), seq))),
        branch = paste0(gs(type, "Branch ", "B"), ":", branch)
      )
      which <- match(object(x)@colData[[ "principal_node" ]], branches$name)
      object(x)@colData[[ "branch" ]] <- branches$branch[ which ]
    }
    return(x)
  })

get_branches.mn <- function(x, branches) {
  # branches: list(c("Y_start", "Y_end"))
  if (is.null(names(branches))) {
    names(branches) <- paste0("Branch ", seq_along(branches))
  }
  linkPrin <- monocle3::principal_graph(object(x))[["UMAP"]]
  branches <- lapply(branches,
    function(x) {
      names(igraph::shortest_paths(linkPrin, x[1], x[2])$vpath[[1]])
    })
  branches
}

setMethod("map", signature = c(x = "job_monocle", ref = "job_seurat"),
  function(x, ref, use.x, use.ref, name = "cell_mapped"){
    matched <- match(rownames(object(x)@colData), rownames(object(ref)@meta.data))
    object(x)@colData[[name]] <- as.character(object(ref)@meta.data[[use.ref]])[matched]
    object(x)@colData[[name]] <- ifelse(
      is.na(object(x)@colData[[name]]),
      as.character(object(x)@colData[[use.x]]),
      object(x)@colData[[name]]
    )
    object(x)@colData[[name]] <- factor(object(x)@colData[[name]],
      levels = sort(unique(object(x)@colData[[name]]))
    )
    x$group.by <- name
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_monocle"),
  function(x, ref, cols = c("pseudotime", "branch"))
  {
    if (any(cols == "pseudotime") && !any(colnames(meta) == "pseudotime")) {
      ref <- add_anno(ref)
    }
    meta <- object(ref)@colData
    cols <- cols[ cols %in% colnames(meta) ]
    if (!length(cols)) {
      stop("No available columns in metadata of `ref`.")
    }
    for (i in cols) {
      lst <- nl(rownames(meta), meta[[ i ]])
      if (is.character(meta[[ i ]])) {
        default <- NA_character_
      } else if (is.numeric(meta[[ i ]])) {
        default <- as.numeric("")
      }
      object(x)@meta.data[[ i ]] <- dplyr::recode(
        rownames(object(x)@meta.data), !!!lst, .default = default
      )
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_monocle", ref = "character"),
  function(x, ref, branches = NULL, enrich = NULL, seurat = NULL, HLs = NULL,
    assay = NULL, enrichExtra = NULL, group_by = NULL, use.enrich = c("go", "kegg"), ...)
  {
    message("Plot pseudotime heatmap.")
    use.enrich <- match.arg(use.enrich)
    if (is.null(seurat)) {
      if (is.null(x$sr_sub)) {
        stop("Object Seurat (`x$sr_sub`) should provided.")
      } else {
        seurat <- x$sr_sub@object
      }
    }
    if (!is.null(assay)) {
      seurat@active.assay <- assay
    }
    seurat <- seurat[ rownames(seurat) %in% ref, ]
    pseudotime <- monocle3::pseudotime(x@object)
    title <- ""
    if (!is.null(branches)) {
      # branches: list(c("Y_start", "Y_end"))
      title <- bind(
        vapply(branches, bind, character(1), co = " -> "), co = ";"
      )
      title <- paste0("在 Branch {title}")
      branches <- get_branches.mn(x, branches)
      hpMode <- "normal"
      if (length(branches) == 2) {
        if (branches[[1]][1] == branches[[2]][1]) {
          hpMode <- "212"
        } else if (tail(branches[[1]], 1) == tail(branches[[2]], 1)) {
          hpMode <- "121"
        }
      }
      if (is.null(object(x)@colData$principal_node)) {
        x <- add_anno(x)
      }
      groupCells <- nl(rownames(object(x)@colData), object(x)@colData$principal_node, FALSE)
      groupCells <- lapply(branches,
        function(yn) {
          names(groupCells[ groupCells %in% yn ])
        })
      n <- 0L
      lst <- lapply(groupCells,
        function(cells) {
          n <<- n + 1L
          pseudotime <- pseudotime[ names(pseudotime) %in% cells ]
          if (hpMode == "212" && n == 1L) {
            rev <- TRUE
          } else if (hpMode == "121" && n == 2L) {
            rev <- TRUE
          } else {
            rev <- FALSE
          }
          prepare_pseudo_heatmap_tidydata(seurat[, cells], pseudotime, rev)
        })
      if (length(lst) > 1) {
        dat <- frbind(lst, idcol = "Branch")
        n <- 0L
        levels <- unlist(lapply(lst,
            function(x) {
              n <<- n + 1L
              paste0("Branch ", n, ":", levels(x[[ ".Pseudo_Time" ]]))
            }))
        dat <- dplyr::mutate(dat,
          .Pseudo_Time = factor(paste0(Branch, ":", .Pseudo_Time), levels = levels)
        )
        dat <- dplyr::group_by(dat, Branch)
      } else if (length(lst)){
        dat <- lst[[1]]
      }
      p.hp <- plot_pseudo_heatmap.seurat(dat, enrich = enrich, HLs = HLs, enrichExtra = enrichExtra,
        group_by = group_by, use.enrich = use.enrich, ...)
      p.hp <- wrap(p.hp, 5 + length(lst) * 3, 8)
    } else {
      p.hp <- plot_pseudo_heatmap.seurat(
        prepare_pseudo_heatmap_tidydata(seurat, pseudotime),
        enrich = enrich, HLs = HLs, enrichExtra = enrichExtra,
        group_by = group_by, use.enrich = use.enrich, ...
      )
      p.hp <- wrap(p.hp, 5, 8)
    }
    # ComplexHeatmap::Heatmap
    # ComplexHeatmap::top_annotation
    # ComplexHeatmap::HeatmapAnnotation
    p.hp <- .set_lab(p.hp, sig(x), "Pseudotime heatmap of genes")
    p.hp <- setLegend(p.hp, "为基因 ({less(ref)}) 表达{title}拟时热图。")
    # .append_heading("Pseudotime Heatmap")
    p.hp
  })

map.fluxGene <- function(x, ref, branches = NULL, enrich = NULL, seurat = NULL, HLs = NULL,
  assay = NULL, enrichExtra = NULL, group_by = NULL, ...)
{
  if (!is(ref, 'df')) {
    stop("The value got TRUE: `!is(ref, 'df')`")
  }
  belong.flux <- reframe_col(dplyr::select(ref, gene, name), "gene",
    function(x) unlist(strsplit(unlist(x), "-")))
  belong.flux <- dplyr::relocate(belong.flux, gene, Metabolic_flux = name)
  if (!is.null(enrichExtra)) {
    message("`enrichExtra` overwrite by `belong.flux`.")
  }
  enrichExtra <- belong.flux
  ref <- belong.flux$gene
  map(x, ref, branches = branches, enrich = enrich, seurat = seurat, HLs = HLs,
    assay = assay, enrichExtra = enrichExtra, group_by = group_by, ...)
}

prepare_pseudo_heatmap_tidydata <- function(seurat, pseudotime, rev.pseudotime = FALSE)
{
  dat <- pseudotime_heatmap(seurat,
    show_rownames = TRUE, pseudotime = pseudotime
  )
  dat <- as_tibble(dat)
  dat <- tidyr::pivot_longer(dat, -rownames, names_to = "Pseudo_Time", values_to = "Levels")
  dat <- dplyr::mutate(dat, Pseudo_Time = as.integer(Pseudo_Time),
    .Pseudo_Time = factor(Pseudo_Time,
      levels = sort(unique(Pseudo_Time), decreasing = rev.pseudotime))
  )
  dat
}

plot_pseudo_heatmap.seurat <- function(dat, enrich = NULL,
  use.enrich = c("go", "kegg"),
  top.enrich = 10, cutoff.enrich = .05, split = 3, HLs = NULL,
  enrichExtra = NULL, group_by = NULL, ...)
{
  .check_columns(dat, c("rownames", ".Pseudo_Time", "Pseudo_Time", "Levels"))
  use.enrich <- match.arg(use.enrich)
  if (!is.null(enrich)) {
    if (is(enrich, "job_enrich")) {
      use.enrich <- paste0("res.", use.enrich)
      enrich <- enrich@tables$step1[[ use.enrich ]][[1]]
      if (use.enrich == "res.go") {
        enrich <- split(enrich, ~ ont)
        enrich <- lapply(enrich,
          function(x) {
            x <- dplyr::filter(x, p.adjust < !!cutoff.enrich)
            x <- head(x, top.enrich)
            x <- dplyr::reframe(dplyr::group_by(x, Description), genes = unlist(geneName_list))
          })
        names(enrich) <- paste0("GO_", names(enrich))
      } else if (use.enrich == "res.kegg") {
        .enrich <- enrich
        enrich <- dplyr::filter(enrich, p.adjust < !!cutoff.enrich)
        if (!nrow(enrich)) {
          message("No results for `cutoff.enrich`: ", cutoff.enrich, "\nCancel filter.")
          enrich <- .enrich
        }
        enrich <- head(enrich, top.enrich)
        enrich <- dplyr::reframe(dplyr::group_by(enrich, Description), genes = unlist(geneName_list))
        enrich <- list(KEGG = enrich)
      }
      if (!is.null(enrichExtra)) {
        message("Use first column of `enrichExtra` as ID.")
        enrichExtra <- dplyr::rename(enrichExtra, genes = 1)
        enrichExtra <- tidyr::pivot_longer(enrichExtra, -genes, names_to = "ont", values_to = "Description")
        enrichExtra <- split(enrichExtra, ~ ont)
        enrich <- c(enrichExtra, enrich)
      }
      for (col in names(enrich)) {
        ## Here, the higher rank pathway will be priority to match
        dat <- map(dat, "rownames", enrich[[ col ]], "genes", "Description", col = col)
        dat[[ col ]] <- ifelse(is.na(dat[[ col ]]), "No match", dat[[ col ]])
        dat[[ col ]] <- factor(dat[[ col ]], levels = c(unique(enrich[[ col ]]$Description), "No match"))
        dat[[ col ]] <- droplevels(dat[[ col ]])
      }
      anno_enrich <- names(enrich)
    }
  } else {
    use.enrich <- ""
  }
  if (!is.null(HLs)) {
    if (!is(HLs, "list") || is.null(names(HLs)) || !all(vapply(HLs, is.character, logical(1)))) {
      stop("`HLs` should be 'list' with names.")
    }
    HLs <- split(as_df.lst(HLs, "type", "genes"), ~ type)
    for (i in seq_along(HLs)) {
      col <- names(HLs)[[ i ]]
      dat <- map(dat, "rownames", HLs[[i]], "genes", "type", col = col)
      dat[[ col ]] <- ifelse(is.na(dat[[ col ]]), "No match", "Match")
      dat[[ col ]] <- factor(dat[[ col ]], levels = c("Match", "No match"))
      dat[[ col ]] <- droplevels(dat[[ col ]])
    }
  }
  maxBreak <- max(ceiling(abs(range(dat$Levels))))
  if (!is.null(group_by)) {
    if (length(group_by) == 1 && any(names(HLs) == group_by)) {
      dat <- dplyr::mutate(dat,
        dplyr::across(!!rlang::sym(group_by),
          function(x) {
            x <- as.character(x)
            ifelse(x == "Match", group_by, "Others")
          }))
    }
    HLs <- HLs[ !names(HLs) %in% group_by ]
    if (is.character(group_by)) {
      group_by <- list(group_by)
    }
    if (!any(dplyr::groups(dat) %in% unlist(group_by))) {
      group_by <- c(group_by, dplyr::groups(dat))
    }
    dat <- dplyr::group_by(dat, !!!rlang::syms(group_by))
  }
  p.hp <- tidyHeatmap::heatmap(dat, rownames, .Pseudo_Time, Levels,
    cluster_columns = FALSE, cluster_rows = TRUE,
    row_title = "Genes", column_title = character(0),
    row_km = split, palette_value = fun_color(-maxBreak, maxBreak),
    show_column_names = FALSE,
    ...
  )
  maxBreak <- max(ceiling(abs(range(dat$Pseudo_Time))))
  # p.hp <- tidyHeatmap::annotation_line(p.hp, Pseudo_Time, show_annotation_name = FALSE)
  p.hp <- tidyHeatmap::annotation_tile(p.hp, Pseudo_Time, show_annotation_name = FALSE)
  # p.hp@top_annotation$color[[ which(p.hp@top_annotation$col_name == "Pseudo_Time") ]] <- fun_color(0, maxBreak, TRUE, "seq")
  if (!is.null(enrich)) {
    dat <- dplyr::ungroup(dat)
    allAnno <- lapply(dplyr::select(dat, dplyr::all_of(anno_enrich)), levels)
    palAnno <- head(color_set(TRUE), length(unlist(allAnno)))
    palAnno <- split(palAnno, rep(seq_along(allAnno), lengths(allAnno)))
    palAnno <- lapply(seq_along(palAnno),
      function(n) {
        nl(allAnno[[n]], palAnno[[n]], FALSE)
      })
    names(palAnno) <- anno_enrich
    pals <- list()
    for (i in anno_enrich) {
      pal <- palAnno[[ i ]]
      pal[[ "No match" ]] <- "white"
      p.hp <- tidyHeatmap::annotation_tile(p.hp, !!rlang::sym(i), palette = pal)
      pals[[ i ]] <- pal
    }
    if (TRUE) {
      # as 'tidyHeatmap::annotation_tile' not support (or bug?) for named vector for color map
      p.hp@left_annotation$color[p.hp@left_annotation$col_name %in% anno_enrich] <- pals
    }
  }
  if (!is.null(HLs)) {
    pal <- c("Match" = "black", "No match" = "white")
    for (i in names(HLs)) {
      p.hp <- tidyHeatmap::annotation_tile(p.hp, !!rlang::sym(i), palette = pal)
    }
    p.hp@left_annotation$color[p.hp@left_annotation$col_name %in% names(HLs)] <- rep(list(pal), length(HLs))
  }
  p.hp
}

get_earliest_principal_nodes <- function(cds, col, target) {
  closest_vertex <- cds@principal_graph_aux[[ "UMAP" ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[ colnames(cds), ])
  if (!missing(col)) {
    cell_ids <- which(SummarizedExperiment::colData(cds)[, col] == target)
    pos <- as.numeric(names(which.max(table(closest_vertex[ cell_ids, ]))))
  } else {
    pos <- as.numeric(names(which.max(table(closest_vertex[,]))))
  }
  igraph <- monocle3::principal_graph(cds)[["UMAP"]]
  vertical <- igraph::V(igraph)
  root_pr_nodes <- vertical$name[ pos ]
  root_pr_nodes
}


cal_modules.cds <- function(cds, gene_sigs, cell_group) {
  n <- 0
  e(lapply(gene_sigs,
    function(genes) {
      n <<- n + 1
      if (length(genes) < 10)
        return(data.frame())
      anno <- paste0("with ", names(gene_sigs)[n])
      res <- try(
        monocle3::find_gene_modules(cds[genes, ], cores = 4, resolution = 10 ^ seq(-6,-1)),
        silent = FALSE
      )
      if (!inherits(res, "try-error")) {
        aggregate <- monocle3::aggregate_gene_expression(cds, res, cell_group)
        row.names(aggregate) <- paste0("Module ", row.names(aggregate))
        module <- as_tibble(res)
      } else {
        module <- tibble::tibble()
        aggregate <- data.frame()
      }
      namel(module, aggregate)
    }))
}

.get_data.plot_genes_in_pseudotime <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
  ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
  trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
  vertical_jitter = NULL, horizontal_jitter = NULL) 
{
  colData <- SummarizedExperiment::colData
  rowData <- SummarizedExperiment::rowData
  fit_models <- monocle3::fit_models
  colData(cds_subset)$pseudotime <- monocle3::pseudotime(cds_subset)
  f_id <- NA
  Cell <- NA
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  if (is.null(min_expr)) {
    min_expr <- 0
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- colData(cds_subset)
  cds_rowData <- rowData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", 
    by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", 
    by.y = "row.names")
  cds_exprs$adjusted_expression <- cds_exprs$expression
  if (label_by_short_name) {
    if (!is.null(cds_exprs$gene_short_name)) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  new_data <- as.data.frame(colData(cds_subset))
  new_data$Size_Factor <- 1
  model_tbl <- fit_models(cds_subset, model_formula_str = trend_formula)
  model_tbl
}


setMethod("asjob_seurat", signature = c(x = "job_monocle"),
  function(x, k, rename = NULL, reset_palette = TRUE){
    x <- regroup(x$sr_sub, x$cellClass_tree.gene_module, k, rename = rename)
    if (reset_palette)
      x$palette <- NULL
    show(vis(x, "regroup.hclust", 1.5))
    return(x)
  })

setMethod("vis", signature = c(x = "job_monocle"),
  function(x, refs, group.by = x@param$group.by,
    use = "logcounts", rownames = TRUE, rownames.size = 5, smooth = TRUE)
  {
    if (!is(refs, "list")) {
      refs <- list(Value = refs)
    }
    if (is.null(names(refs))) {
      names(refs) <- paste0("Hm ", seq_along(refs))
    }
    cell_order <- order(e(monocle3::pseudotime(object(x))))
    hp.anno <- ComplexHeatmap::HeatmapAnnotation(
      Group = object(x)@colData[[ group.by ]][ cell_order ],
      col = list(Group = .setPaletteForMN(x, group.by, ggsci::pal_npg()(10)))
    )
    n <- 0L
    h.lst <- lapply(refs,
      function(ref) {
        n <<- n + 1L
        expr <- object(x)@assays@data[[ use ]][ ref, cell_order ]
        if (smooth) {
          expr <- e(apply(expr, 1, function(x) stats::smooth.spline(x, df = 3)$y))
          expr <- t(scale(expr))
        } else {
          expr <- t(scale(Matrix::t(expr)))
        }
        params <- list(
          expr, name = names(refs)[n], km = 2,
          col = .get_col_fun(expr),
          cluster_columns = FALSE,
          row_names_gp = gpar(fontsize = rownames.size),
          show_column_names = FALSE,
          show_row_names = rownames,
          clustering_method_rows = "ward.D2",
          heatmap_legend_param = list(title = names(refs)[n], ncol = 1)
        )
        if (n == 1) {
          params <- c(params, list(top_annotation = hp.anno))
        }
        do.call(ComplexHeatmap::Heatmap, params)
      })
    hs <- h.lst[[1]]
    if (length(h.lst) >= 2) {
      for (i in 2:length(h.lst)) {
        hs <- ComplexHeatmap::`%v%`(hs, h.lst[[ i ]])
      }
    }
    return(hs)
  })

.setPaletteForMN <- function(x, group.by = x@params$group.by, palette = color_set()) {
  groups <- x@object@colData[[ group.by ]]
  if (is.factor(groups)) {
    levels <- levels(groups)
  } else {
    levels <- sort(unique(groups))
  }
  nl(levels, palette[seq_along(levels)], FALSE)
}

.get_col_fun <- function(data) {
  range <- range(data)
  circlize::colorRamp2(seq(from = floor(range[1]), to = ceiling(range[2]), length = 11),
    rev(RColorBrewer::brewer.pal(11, "RdYlGn")))
}

setMethod("skel", signature = c(x = "job_monocle"),
  function(x, suffix, pattern, sig.mn = paste0("mn.", suffix),
    sig.sr = paste0("sr.", suffix),
    sig.sr_sub = paste0("sr_sub.", suffix))
  {
    code <- c('',
      paste0('mn <- do_monocle(sr, "', pattern, '")'),
      '',
      'mn <- step1(mn, "cell_type", pre = TRUE)',
      'mn@plots$step1$p.prin',
      'mn <- step2(mn, "Y_2")',
      'mn@plots$step2$p.pseu',
      'mn <- step3(mn, group.by = "seurat_clusters")',
      'mn@plots$step3$gene_module_heatdata$graph_test.sig',
      '',
      paste0('sr_sub <- asjob_seurat(mn, 5, rename = "', pattern, "_", suffix, '")'),
      'vis(sr_sub, "regroup.hclust")',
      '',
      'mn <- clear(mn)',
      '',
      'sr <- map(sr, sr_sub, "scsa_cell", "regroup.hclust")',
      'vis(sr, "cell_mapped")',
      'vis(sr, "scsa_cell")',
      '',
      'sr <- clear(sr)'
    )
    code <- gs(code, "\\bsr\\b", sig.sr)
    code <- gs(code, "\\bsr_sub\\b", sig.sr_sub)
    code <- gs(code, "\\bmn\\b", sig.mn)
    writeLines(code)
  })



