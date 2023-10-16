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
    info = c("Tutorial: https://cole-trapnell-lab.github.io/monocle3/docs/getting_started/")
    ))

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
    object <- e(SeuratWrappers::as.cell_data_set(object(x), group.by = group.by, ...))
    mn <- .job_monocle(object = object)
    object(mn) <- e(monocle3::estimate_size_factors(object(mn)))
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- object(x)@reductions[["pca"]]@feature.loadings
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- object(x)@reductions$pca@stdev
    mn@params$group.by <- group.by
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
  function(x, groups = x@params$group.by, pt.size = .7){
    step_message("Constructing single-cell trajectories.
      red{{`groups`}} would passed to `monocle3::plot_cells` for
      annotation in plot. Mutilple group could be given.
      "
    )
    x$pt.size <- pt.size
    if (!all(groups %in% colnames(object(x)@colData)))
      stop("Some of `groups` not found in `colData` of `object(x)`")
    object(x) <- e(monocle3::cluster_cells(object(x)))
    object(x) <- e(monocle3::learn_graph(object(x)))
    p.traj <- e(sapply(groups, simplify = F,
        function(group) {
          p <- monocle3::plot_cells(object(x), color_cells_by = group,
            label_cell_groups = T, label_branch_points = T,
            group_label_size = 4, graph_label_size = 2,
            cell_size = x$pt.size, cell_stroke = 0, alpha = .7
          )
          p + scale_color_manual(values = color_set())
        }))
    p.prin <- monocle3::plot_cells(
      object(x), color_cells_by = groups[1],
      label_cell_groups = F, label_principal_points = T,
      graph_label_size = 3, cell_size = x$pt.size, cell_stroke = 0,
      alpha = .7
    )
    p.prin <- p.prin + scale_color_manual(values = color_set())
    p.prin <- wrap(p.prin, 10, 7)
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
    x@plots[[ 2 ]] <- namel(p.pseu)
    return(x)
  })

setMethod("step3", signature = c(x = "job_monocle"),
  function(x, formula_string = NULL, cores = 4){
    step_message("This step do:
      1. Regression analysis;
      2. Graph-autocorrelation analysis;
      3. Finding modules of co-regulated genes (Significant genes).
      4. Significant genes in co-regulated modules.
      "
    )
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
    x$cellClass_tree.gene_module <- hclust(dist(t(gene_module$graph_test.sig$aggregate)))
    return(x)
  })

setMethod("step4", signature = c(x = "job_monocle"),
  function(x, groups, genes){
    step_message("Plot genes (in branch) that change as a function of pseudotime.
      red{{`groups`}} and red{{`genes`}} were used to subset the `object(x)`."
    )
    cds <- object(x)
    fun_sub <- selectMethod("[", class(cds))
    gene.groups <- grouping_vec2list(genes, 10, T)
    pblapply <- pbapply::pblapply
    colData <- SummarizedExperiment::colData
    genes_in_pseudotime <- e(pblapply(gene.groups,
        function(genes) {
          cds <- fun_sub(cds,
            rownames(cds) %in% genes,
            colData(cds)[[ x@params$group.by ]] %in% groups
          )
          p <- monocle3::plot_genes_in_pseudotime(cds,
            label_by_short_name = F,
            color_cells_by = x@params$group.by,
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
      ref[] <- dplyr::recode(ref, !!!(rename),
        .default = as.character(unname(ref)))
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
        silent = T
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
