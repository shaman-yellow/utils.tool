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
  function(x, group.by = "SingleR_cell", ...){
    step_message("
      Other parameters would be passed to `SeuratWrappers::as.cell_data_set`.
      <http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html>
      "
    )
    if (x@step < 3) {
      stop("x@step < 3. At least step 3 was needed with `Seurat::FindClusters`.")
    }
    if (group.by == "SingleR_cell") {
      if (x@step < 4L)
        stop("`step4(x)` has not been performed.")
    }
    object <- e(SeuratWrappers::as.cell_data_set(object(x), group.by = group.by, ...))
    mn <- .job_monocle(object = object)
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- object(x)@reductions[["pca"]]@feature.loadings
    if (is.null(object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']]))
      object(mn)@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- object(x)@reductions$pca@stdev
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
  function(x, groups = "SingleR_cell"){
    step_message("Constructing single-cell trajectories.
      red{{`groups`}} would passed to `monocle3::plot_cells` for
      annotation in plot. Mutilple group could be given.
      "
    )
    if (!all(groups %in% colnames(object(x)@colData)))
      stop("Some of `groups` not found in `colData` of `object(x)`")
    object(x) <- e(monocle3::cluster_cells(object(x)))
    object(x) <- e(monocle3::learn_graph(object(x)))
    p.traj <- e(sapply(groups, simplify = F,
        function(group) {
          p <- monocle3::plot_cells(object(x), color_cells_by = group,
            label_cell_groups = T, label_branch_points = T,
            group_label_size = 4, graph_label_size = 2
          )
          p + scale_color_manual(values = color_set())
        }))
    p.prin <- monocle3::plot_cells(
      object(x), color_cells_by = groups[1],
      label_cell_groups = F, label_principal_points = T,
      graph_label_size = 3
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
        cell_size = .5
        ))
    x@plots[[ 2 ]] <- namel(p.pseu)
    return(x)
  })

setMethod("step3", signature = c(x = "job_monocle"),
  function(x, formula_string){
    step_message("This step do:
      1. Regression analysis;
      2. Graph-autocorrelation analysis grey{{(see monocle3::graph_test for help)}};
      3. Finding modules of co-regulated genes.
      "
    )
    if (!is.character(formula_string)) {
      stop("is.character(formula_string) == F")
    }
    fits <- e(monocle3::fit_models(object(x), formula_string, cores = 4))
    fit_coefs <- e(monocle3::coefficient_table(fits))
    fit_coefs.sig <- dplyr::filter(fit_coefs, term != "(Intercept)", q_value < .05)
      # goodness <- e(monocle3::evaluate_fits(gene_fits))
      # fit_coefs.sig <- dplyr::select(fit_coefs.si, gene = gene_short_name, term, q_value, estimate)
    graph_test <- e(monocle3::graph_test(object(x), neighbor_graph = "knn", cores = 4))
    graph_test <- as_tibble(graph_test)
    gene_module <- try(gene_module <- e(monocle3::find_gene_modules(object(x), cores = 4)), T)
    if (!inherits(gene_module, "try-error")) {
      gene_module <- as_tibble(gene_module)
    } else {
      gene_module <- tibble::tibble()
    }
    # x@tables[[ 3 ]] <- namel(fit_coefs, fit_coefs.sig, graph_test, gene_module)
    return(x)
  })

setMethod("step4", signature = c(x = "job_monocle"),
  function(x){
    step_message("Finding genes that change as a function of pseudotime.")
  })

setMethod("ids", signature = c(x = "job_monocle", id = "missing"),
  function(x, ...){
    ids(x, id = "SingleR_cell", ...)
  })

setMethod("ids", signature = c(x = "job_monocle", id = "character"),
  function(x, id, unique = T){
    ids <- SummarizedExperiment::colData(object(x))[[ id ]]
    if (unique)
      unique(ids)
  })

get_principal_nodes <- function(cds, col, target) {
  cell_ids <- which(SummarizedExperiment::colData(cds)[, col] == target)
  closest_vertex <- cds@principal_graph_aux[[ "UMAP" ]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[ colnames(cds), ])
  pos <- as.numeric(names(which.max(table(closest_vertex[ cell_ids, ]))))
  root_pr_nodes <- igraph::V(monocle3::principal_graph(cds)[["UMAP"]])$name[ pos ]
  root_pr_nodes
}
