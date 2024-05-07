# ==========================================================================
# all these function were extract from monocle (not Monocle 3)
# but revised for Seurat (5) and Monocle3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pseudotime_heatmap <- function(
    object,
    cluster_rows = TRUE,
    hclust_method = "ward.D2",
    num_clusters = 6,
    hmcols = NULL,
    add_annotation_row = NULL,
    add_annotation_col = NULL,
    show_rownames = FALSE,
    use_gene_short_name = TRUE,
    scale_max = 3,
    scale_min = -3,
    trend_formula = "~ sm.ns(Pseudotime, df = 3)",
    cores = 1,
    nCol = 100,
    pseudotime = NULL,
    plot_heatmap = F)
{
  ## NOTE: object of cell_data_set converted by 'SeuratWrapper' may caused error.' 
  if (is(object, "cell_data_set")) {
    stop("`object` can not be cell_data_set.")
  } else if (is(object, "Seurat")) {
    if (!is.null(pseudotime)) {
      object@meta.data[[ "Pseudotime" ]] <- pseudotime
    }
  }
  if (is.null(object@meta.data[[ "Pseudotime" ]]) && is.null(pseudotime)) {
    stop("`Pseudotime` should be provided.")
  }
  if (exists("sm.ns", envir = .GlobalEnv)) {
    fun <- get("sm.ns", envir = .GlobalEnv)
    if (!identical(fun, VGAM::sm.ns, ignore.environment = T)) {
      stop("`sm.ns` in .GlobalEnv should be VGAM::sm.ns")
    }
  } else {
    sm.ns <<- VGAM::sm.ns
  }
  num_clusters <- min(num_clusters, nrow(object))
  if (is.null(.metaData(object)$Pseudotime)) {
    stop("is.null(.metaData(object)$Pseudotime)")
  }
  newdata <- data.frame(
    Pseudotime = seq(min(.metaData(object)$Pseudotime),
      max(.metaData(object)$Pseudotime),
      length.out = nCol
    )
  )
  m <- genSmoothCurves(object, new_data = newdata,
    trend_formula = trend_formula, relative_expr = T)
  # remove genes with no expression in any condition
  m <- m[apply(m, 1, sum) != 0, ]
  # Row-center the data.
  m <- m[apply(m, 1, sd) != 0, ]
  m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
  m <- m[!is.na(row.names(m)), ]
  m[is.nan(m)] <- 0
  m[m > scale_max] <- scale_max
  m[m < scale_min] <- scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix))) / 2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- blue2green2red(length(bks) - 1)
  } else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap::pheatmap(heatmap_matrix,
    useRaster = T,
    cluster_cols = FALSE,
    cluster_rows = cluster_rows,
    show_rownames = F,
    show_colnames = F,
    clustering_distance_rows = row_dist,
    clustering_method = hclust_method,
    cutree_rows = num_clusters,
    silent = TRUE,
    filename = NA,
    breaks = bks,
    border_color = NA,
    color = hmcols
  )
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, num_clusters)))
  } else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    if (nrow(add_annotation_col) != nCol) {
      stop("add_annotation_col should have only ", nCol, " rows!")
    }
    annotation_col <- add_annotation_col
  } else {
    annotation_col <- NA
  }
  feature_label <- row.names(heatmap_matrix)
  if (!is.null(annotation_row)) {
    row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row)) {
    row.names(annotation_row) <- row_ann_labels
  }
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  if (plot_heatmap) {
    ph_res <- pheatmap::pheatmap(heatmap_matrix[, ], # ph$tree_row$order
      useRaster = T,
      cluster_cols = FALSE, cluster_rows = cluster_rows,
      show_rownames = show_rownames, show_colnames = F,
      clustering_distance_rows = row_dist, # row_dist
      clustering_method = hclust_method, # ward.D2
      cutree_rows = num_clusters,
      # cutree_cols = 2,
      annotation_row = annotation_row, annotation_col = annotation_col,
      treeheight_row = 20, breaks = bks, fontsize = 6,
      color = hmcols, border_color = NA, silent = TRUE, filename = NA
    )
  } else {
    return(heatmap_matrix)
  }
}


genSmoothCurves <- function(
    object, new_data, trend_formula = "~ sm.ns(Pseudotime, df = 3)",
    relative_expr = T, response_type = "response")
{
  if (is(object, "Seurat") || is(object, "cell_data_set")) {
    expressionFamily <- VGAM::uninormal()
  } else {
    stop("...")
  }
  expression_curve_matrix <- smartEsApply(object, 1,
    function(x, trend_formula, expressionFamily, relative_expr, new_data) {
      environment(fit_model_helper) <- environment()
      environment(responseMatrix) <- environment()
      model_fits <- fit_model_helper(x,
        modelFormulaStr = trend_formula, expressionFamily = expressionFamily,
        relative_expr = relative_expr, disp_func = NULL
      )
      if (is.null(model_fits)) {
        expression_curve <- as.data.frame(matrix(rep(NA, nrow(new_data)), nrow = 1))
      } else {
        res <- responseMatrix(list(model_fits), new_data, response_type = response_type)
        expression_curve <- as.data.frame(res)
      }
      colnames(expression_curve) <- row.names(new_data)
      expression_curve
    },
    convert_to_dense = TRUE,
    trend_formula = trend_formula, expressionFamily = expressionFamily, relative_expr = relative_expr, new_data = new_data
  )
  expression_curve_matrix <- as.matrix(do.call(rbind, expression_curve_matrix))
  return(expression_curve_matrix)
}

fit_model_helper <- function(
    x, modelFormulaStr, expressionFamily,
    relative_expr, disp_func = NULL, verbose = FALSE, ...)
{
  modelFormulaStr <- paste("f_expression", modelFormulaStr, sep = "")
  orig_x <- x
  if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    if (relative_expr) {
      x <- x / Size_Factor
    }
    f_expression <- round(x)
  } else if (expressionFamily@vfamily %in% c("uninormal", "binomialff")) {
    f_expression <- x
  } else {
    f_expression <- log10(x)
  }
  tryCatch(
    {
      if (verbose) {
        FM_fit <- VGAM::vglm(as.formula(modelFormulaStr),
          family = expressionFamily, epsilon = 1e-1
        )
      } else {
        FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr),
          family = expressionFamily, epsilon = 1e-1
        ))
      }
      FM_fit
    },
    error = function(e) {
      print(e)
      # If we threw an exception, re-try with a simpler model.  Which one depends on
      # what the user has specified for expression family
      # print(disp_guess)
      backup_expression_family <- NULL
      if (expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
        stop("The `disp_func` is NULL.")
      } else if (expressionFamily@vfamily %in% c("uninormal")) {
        backup_expression_family <- NULL
      } else if (expressionFamily@vfamily %in% c("binomialff")) {
        backup_expression_family <- NULL
      } else {
        backup_expression_family <- NULL
      }
      if (!is.null(backup_expression_family)) {
        # FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family=backup_expression_family, trace=T, epsilon=1e-1, checkwz=F)
        test_res <- tryCatch(
          {
            if (verbose) {
              FM_fit <- VGAM::vglm(as.formula(modelFormulaStr), family = backup_expression_family, epsilon = 1e-1, checkwz = TRUE)
            } else {
              FM_fit <- suppressWarnings(VGAM::vglm(as.formula(modelFormulaStr), family = backup_expression_family, epsilon = 1e-1, checkwz = TRUE))
            }
            FM_fit
          },
          # warning = function(w) { FM_fit },
          error = function(e) {
            # print (e);
            NULL
          }
        )
        # print(test_res)
        test_res
      } else {
        # print(e);
        NULL
      }
    }
  )
}

responseMatrix <- function(models, newdata = NULL, response_type = "response", cores = 1)
{
  res_list <- parallel::mclapply(models, function(x) {
    if (is.null(x)) {
      NA
    } else {
      if (x@family@vfamily %in% c("negbinomial", "negbinomial.size")) {
        predict(x, newdata = newdata, type = response_type)[, "mean"]
      } else if (x@family@vfamily %in% c("uninormal")) {
        predict(x, newdata = newdata, type = response_type)[, "mean"]
      } else {
        10^predict(x, newdata = newdata, type = response_type)[, "mean"]
      }
    }
  }, mc.cores = cores)
  res_list_lengths <- lapply(
    res_list[!is.na(res_list)],
    length
  )
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0) {
    na_matrix <- matrix(rep(
      rep(NA, res_list_lengths[[1]]),
      num_na_fits
    ), nrow = num_na_fits)
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    non_na_matrix <- Matrix::t(do.call(cbind,
        lapply(res_list[!is.na(res_list)], unlist)))
    row.names(non_na_matrix) <- names(res_list[!is.na(res_list)])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list), ]
  } else {
    res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[!is.na(res_list)])
  }
  res_matrix
}

smartEsApply <- function(X, MARGIN, FUN, convert_to_dense, ...) {
  parent <- environment(FUN)
  if (is.null(parent)) {
    parent <- emptyenv()
  }
  e1 <- new.env(parent = parent)
  Biobase::multiassign(names(.metaData(X)), .metaData(X), envir = e1)
  environment(FUN) <- e1
  if (is(X, "Seurat")) {
    exprs <- function(x) {
      x@assays[[ x@active.assay ]]@scale.data
    }
  } else if (is(X, "cell_data_set")) {
    exprs <- function(x) {
      x@assays@data$logcounts
    }
  } else {
    stop("`X` should be 'Seurat' or 'cell_data_set'.")
  }
  dat <- exprs(X)
  if (isSparseMatrix(dat)) {
    dat <- as.data.frame(dat)
  }
  res <- apply(dat, MARGIN, FUN, ...)
  if (MARGIN) {
    names(res) <- row.names(X)
  } else {
    names(res) <- colnames(X)
  }
  res
}

isSparseMatrix <- function(x) {
  is(x, "dgCMatrix") || is(x, "dgTMatrix")
}

.metaData <- function(x) {
  if (is(x, "Seurat")) {
    x@meta.data
  } else if (is(x, "cell_data_set")) {
    x@colData
  } else {
    stop("`x` is not 'Seurat' or 'cell_data_set'")
  }
}

blue2green2red <- function(n) {
  rgb.tables(n,
    red = c(0.8, 0.2, 1),
    green = c(0.5, 0.4, 0.8),
    blue = c(0.2, 0.2, 1))
}

rgb.tables <- function(n,
  red = c(0.75, 0.25, 1),
  green = c(0.5, 0.25, 1),
  blue = c(0.25, 0.25, 1))
{
  rr <- do.call("table.ramp", as.list(c(n, red)))
  gr <- do.call("table.ramp", as.list(c(n, green)))
  br <- do.call("table.ramp", as.list(c(n, blue)))
  rgb(rr, gr, br)
}

table.ramp <- function(n, mid = 0.5, sill = 0.5, base = 1, height = 1)
{
    x <- seq(0, 1, length.out = n)
    y <- rep(0, length(x))
    sill.min <- max(c(1, round((n - 1) * (mid - sill / 2)) + 1))
    sill.max <- min(c(n, round((n - 1) * (mid + sill / 2)) + 1))
    y[sill.min:sill.max] <- 1
    base.min <- round((n - 1) * (mid - base / 2)) + 1
    base.max <- round((n - 1) * (mid + base / 2)) + 1
    xi <- base.min:sill.min
    yi <- seq(0, 1, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    xi <- sill.max:base.max
    yi <- seq(1, 0, length.out = length(xi))
    i <- which(xi > 0 & xi <= n)
    y[xi[i]] <- yi[i]
    height * y
}
