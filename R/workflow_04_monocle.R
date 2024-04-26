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
    method = "R package `Monocle3` used for cell pseudotime analysis"
    ))

setGeneric("do_monocle", 
  function(x, ref, ...) standardGeneric("do_monocle"))

setMethod("do_monocle", signature = c(x = "job_seurat", ref = "character"),
  function(x, ref, dims = 1:15, resolution = 1.2, group.by = x@params$group.by)
  {
    x <- getsub(x, cells = grp(x@object@meta.data[[group.by]], ref))
    x@step <- 2L
    sr_sub <- step3(x, dims, resolution)
    x <- asjob_monocle(sr_sub, "seurat_clusters")
    x$sr_sub <- sr_sub
    return(x)
  })

setGeneric("asjob_monocle", 
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
      stop("is.null(group.by) == T")
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
  function(x, groups = x@params$group.by, pt.size = 1.5, pre = F, norm_method = "none"){
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
    object(x) <- e(monocle3::learn_graph(object(x)))
    p.traj <- e(sapply(groups, simplify = F,
        function(group) {
          p <- monocle3::plot_cells(object(x), color_cells_by = group,
            label_cell_groups = T, label_branch_points = T,
            group_label_size = 4, graph_label_size = 2,
            cell_size = x$pt.size, cell_stroke = 0, alpha = .7
          )
          p <- wrap(p + scale_color_manual(values = palette))
          p <- .set_lab(p, sig(x), "trajectories of", group)
        }))
    p.prin <- monocle3::plot_cells(
      object(x), color_cells_by = groups[1],
      label_cell_groups = F, label_principal_points = T,
      graph_label_size = 3, cell_size = x$pt.size, cell_stroke = 0,
      alpha = .7
    )
    p.prin <- p.prin + scale_color_manual(values = palette)
    p.prin <- wrap(p.prin, 10, 7)
    p.prin <- .set_lab(p.prin, sig(x), "principal points")
    x@plots[[ 1 ]] <- namel(p.traj, p.prin)
    return(x)
  })

setMethod("step2", signature = c(x = "job_monocle"),
  function(x, roots){
    step_message("
      red{{`roots`}} would passed to `monocle3::order_cells` for setting
      as roots of `pseudotime`. If `roots` is character with names, processed
      by red{{`get_principal_nodes`}} and then passed to parameter `root_pr_nodes`;
      else, directly passed to param `root_pr_nodes`.
      "
    )
    if (!is.null(names(roots))) {
      if (length(roots) > 1)
        stop("With name of `roots`, but length(roots) > 1")
      roots <- get_principal_nodes(object(x), names(roots), unname(roots))
    }
    object(x) <- e(monocle3::order_cells(object(x), root_pr_nodes = roots))
    p.pseu <- e(monocle3::plot_cells(
        object(x), color_cells_by = "pseudotime",
        label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = FALSE, graph_label_size = 3,
        cell_size = x$pt.size, cell_stroke = 0, alpha = .7
        ))
    p.pseu <- wrap(p.pseu, 6, 5)
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
        stop("is.character(formula_string) == F")
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
        fits <- e(monocle3::fit_models(object(x), formula_string, cores = cores, verbose = T))
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
    if (is.null(x@tables[[ 3 ]]$graph_test)) {
      graph_test <- e(monocle3::graph_test(object(x), neighbor_graph = "knn", cores = cores))
      graph_test <- as_tibble(graph_test)
      graph_test.sig <- dplyr::filter(graph_test, q_value < .05)
      graph_test.sig <- dplyr::arrange(graph_test.sig, dplyr::desc(morans_I), q_value)
      graph_test.sig <- dplyr::rename(graph_test.sig, gene_id = 1)
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
      gene_module <- try(cal_modules.cds(object(x), gene_sigs, cell_group), T)
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
    x$cellClass_tree.gene_module <- try(hclust(dist(t(gene_module$graph_test.sig$aggregate))), silent = T)
    return(x)
  })

setMethod("step4", signature = c(x = "job_monocle"),
  function(x, groups = ids(x), genes, group.by = NULL)
  {
    step_message("Plot genes (in branch) that change as a function of pseudotime.
      red{{`groups`}} and red{{`genes`}} were used to subset the `object(x)`."
    )
    cds <- object(x)
    fun_sub <- selectMethod("[", class(cds))
    gene.groups <- grouping_vec2list(genes, 10, T)
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
              colData(cds)[[ x@params$group.by ]] %in% groups
            )
          }
          p <- monocle3::plot_genes_in_pseudotime(cds,
            label_by_short_name = F,
            color_cells_by = group.by,
            min_expr = 0.5
          )
          wrap(p, 6, length(genes) * 1.6)
        }))
    names(genes_in_pseudotime) <- paste0("pseudo", 1:length(genes_in_pseudotime))
    x@plots[[ 4 ]] <- namel(genes_in_pseudotime)
    return(x)
  })


setMethod("regroup", signature = c(x = "job_seurat", ref = "hclust"),
  function(x, ref, k, by.name = F, rename = NULL){
    ref <- cutree(ref, k)
    x$cutree <- ref
    regroup(x, ref, k, by.name, rename)
  })

setMethod("regroup", signature = c(x = "job_seurat", ref = "integer"),
  function(x, ref, k, by.name = F, rename = NULL){
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
      re.idents <- nl(names(idents), dplyr::recode(names(idents), !!!ref), F)
    } else {
      re.idents <- nl(names(idents), dplyr::recode(unname(idents), !!!ref), F)
    }
    re.idents <- factor(re.idents, levels = sort(as.character(unique(re.idents))))
    e(SeuratObject::Idents(object(x)) <- re.idents)
    object(x)@meta.data$regroup.hclust <- unname(re.idents)
    return(x)
  })

setMethod("ids", signature = c(x = "job_monocle"),
  function(x, id = x@params$group.by, unique = T){
    ids <- SummarizedExperiment::colData(object(x))[[ id ]]
    if (unique)
      unique(ids)
    else
      ids
  })

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

setMethod("map", signature = c(x = "job_monocle", ref = "character"),
  function(x, ref, cells = NULL, seurat = NULL, ...)
  {
    message("Plot pseudotime heatmap.")
    if (is.null(seurat)) {
      if (is.null(x$sr_sub)) {
        stop("Object Seurat should provided.")
      } else {
        seurat <- x$sr_sub@object
      }
    }
    seurat <- seurat[ rownames(seurat) %in% ref, ]
    if (is.null(cells)) {
      seurat <- seurat[, cells]
    }
    p.hp <- plot_pseudotime_heatmap(seurat,
      show_rownames = T,
      pseudotime = monocle3::pseudotime(x@object)
    )
    p.hp <- .set_lab(p.hp, sig(x), "Pseudotime heatmap of genes")
    p.hp
  })

get_principal_nodes <- function(cds, col, target) {
  cell_ids <- which(SummarizedExperiment::colData(cds)[, col] == target)
  closest_vertex <- cds@principal_graph_aux[[ "UMAP" ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[ colnames(cds), ])
  pos <- as.numeric(names(which.max(table(closest_vertex[ cell_ids, ]))))
  root_pr_nodes <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[ pos ]
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
        silent = F
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
  function(x, k, rename = NULL, reset_palette = T){
    x <- regroup(x$sr_sub, x$cellClass_tree.gene_module, k, rename = rename)
    if (reset_palette)
      x$palette <- NULL
    show(vis(x, "regroup.hclust", 1.5))
    return(x)
  })

setMethod("vis", signature = c(x = "job_monocle"),
  function(x, refs, group.by = x@param$group.by,
    use = "logcounts", rownames = T, rownames.size = 5, smooth = T)
  {
    if (!is(refs, "list")) {
      refs <- list(Value = refs)
    }
    if (is.null(names(refs))) {
      names(refs) <- paste0("Hm ", 1:length(refs))
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
          cluster_columns = F,
          row_names_gp = gpar(fontsize = rownames.size),
          show_column_names = F,
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
  nl(levels, palette[1:length(levels)], F)
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
      'mn <- step1(mn, "cell_type", pre = T)',
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


