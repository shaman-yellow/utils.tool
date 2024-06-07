# ==========================================================================
# workflow of lasso
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_lasso <- setClass("job_lasso", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@EfsAnEnsemblNeuman2017]",
    method = "R package `glmnet` used for LASSO analysis and EFS used for feature selection",
    tag = "filter:lasso",
    analysis = "Lasso 回归"
    ))

setGeneric("asjob_lasso", 
  function(x, ...) standardGeneric("asjob_lasso"))

setMethod("asjob_lasso", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = "gene_name", from_normed = T, fun_scale = function(x) scale(x + 1)){
    step_message("The default, 'job_limma' from 'job_tcga' were adapted to
      convertion.
      "
    )
    if (from_normed) {
      message("Ignore param `fun_scale`")
      if (!is(x@params$normed_data, 'EList'))
        stop("is(x@params$normed_data, 'EList') == F")
      if (!is.null(use.filter)) {
        pos <- x@params$normed_data$genes[[use]] %in% use.filter
        object <- e(limma::`[.EList`(x@params$normed_data, pos, ))
      } else {
        object <- x@params$normed_data
      }
      mtx <- t(object$E)
      colnames(mtx) <- object$genes[[ use ]]
      x <- .job_lasso(object = mtx)
      x@params$metadata <- as_tibble(object$targets)
    } else {
      message("Ignore param `use.filter`")
      if (x@step != 0L) {
        stop("x@step != 0L")
      }
      mtx <- object(x)$counts
      if (!is.null(fun_scale)) {
        mtx <- fun_scale(mtx)
      }
      mtx <- t(mtx)
      colnames(mtx) <- object(x)$genes[[ use ]]
      metadata <- as_tibble(object(x)$samples)
      x <- .job_lasso(object = mtx)
      x@params$metadata <- metadata
    }
    return(x)
  })

setMethod("step0", signature = c(x = "job_lasso"),
  function(x){
    step_message("Prepare your data with function `asjob_lasso`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_lasso"),
  function(x, target = x@params$metadata$group, levels = c("Alive", "Dead"),
    n.train = .8, seed = 555)
  {
    step_message("Data preparation.")
    len <- nrow(object(x))
    all <- 1:len
    set.seed(seed)
    wh.train <- sample(all, round(len * n.train))
    wh.valid <- all[!all %in% wh.train]
    x@params$train <- object(x)[wh.train, ]
    x@params$valid <- object(x)[wh.valid, ]
    x@params$train.target <- target[wh.train]
    x@params$valid.target <- target[wh.valid]
    x@params$levels <- levels
    p.all <- gtitle(new_pie(target)@data, "All")
    p.train <- gtitle(new_pie(x@params$train.target)@data, "Train")
    p.valid <- gtitle(new_pie(x@params$valid.target)@data, "Validate")
    p.consist <- wrap(xf(p.all = 1, p.train = 1, p.valid = 1), 6, 2)
    x@plots[[ 1 ]] <- namel(p.consist)
    return(x)
  })

setMethod("step2", signature = c(x = "job_lasso"),
  function(x, top = 30){
    step_message("Evaluate variable (genes) importance.")
    target <- x@params$train.target
    target <- ifelse(target == x@params$levels[1], 0L, 1L)
    data <- data.frame(cbind(x@params$train, target))
    efs <- e(EFS::ensemble_fs(data, ncol(data), runs = 10))
    colnames(efs) <- colnames(object(x))
    p.efs <- plot_colStack.efs(efs, top = top)
    tops <- attr(p.efs, "tops")
    x@params$tops <- tops
    x@tables[[ 2 ]] <- namel(efs)
    x@plots[[ 2 ]] <- namel(p.efs)
    return(x)
  })

setMethod("step3", signature = c(x = "job_lasso"),
  function(x, family = "binomial", nfold = 10, use_tops = T){
    if (use_tops) {
      fun <- function(mtx) mtx[, colnames(mtx) %in% x@params$tops ]
      x@params$train %<>% fun()
      x@params$valid %<>% fun()
    }
    x@params$model <- e(glmnet::cv.glmnet(
        x@params$train, x@params$train.target,
        alpha = 1, family = family, nfold = nfold
        ))
    model <- x@params$model
    p.model <- wrap(as_grob(expression(plot(model)), environment()), 6, 6)
    levels <- x@params$levels
    res <- sapply(c("lambda.min", "lambda.1se"), simplify = F,
      function(lam) {
        preds <- e(stats::predict(model, newx = x@params$valid, s = model[[ lam ]]))
        roc <- e(pROC::roc(x@params$valid.target, as.numeric(preds), plot = F, levels = levels))
        namel(preds, roc)
      })
    x@params$preds <- lapply(res, function(lst) lst$preds)
    x@params$roc <- lapply(res, function(lst) lst$roc)
    x@params$coef <- coefficients(model, s = c(lambda.min = model$lambda.min,
        lambda.1se = model$lambda.1se))
    p.roc <- lapply(x@params$roc,
      function(roc) {
        plot_roc(roc)
      })
    coef <- as_tibble(Matrix::as.matrix(x@params$coef))
    coef <- arrange(coef, dplyr::desc(abs(lambda.min)))
    coef <- rename(coef, variable = rownames)
    p.coef <- plot_sig(filter(coef, !grepl("Inter", variable)))
    x@plots[[ 3 ]] <- namel(p.model, p.roc, p.coef)
    x@tables[[ 3 ]] <- namel(coef)
    return(x)
  })

plot_roc <- function(roc) {
  wrap(as_grob(
      expression({
        plot(roc)
        text(.1, .1, paste0("AUC: ", round(roc$auc[[1]], 2)), cex = 1.2)
      }),
    environment()), 6, 6)
}

plot_sig <- function(data, x = colnames(data)[1], y = colnames(data)[2],
  lab_x = "Variables", lab_y = "Coefficients") 
{
  data <- filter(data, abs(!!rlang::sym(y)) != 0)
  p <- ggplot(data, aes(x = stats::reorder(!!rlang::sym(x),
        !!rlang::sym(y)), y = !!rlang::sym(y))) +
    geom_point(aes(color = !!rlang::sym(y))) +
    labs(x = lab_x, y = lab_y) +
    theme_minimal() +
    coord_flip() +
    guides(color = "none") +
    theme()
  wrap(p, 6, 3 + nrow(data) * .2)
}

plot_colStack.efs <- function(raw, top = 30,
  lab_x = "Variable", lab_y = "Importance", lab_fill = "Method")
{
  data <- as_tibble(raw)
  data <- rename(data, group = rownames)
  data <- tidyr::gather(data, var, value, -group)
  tops <- sort(vapply(split(data$value, data$var), sum, double(1)), decreasing = T)
  tops <- head(names(tops), top)
  data <- filter(data, var %in% dplyr::all_of(tops))
  if (top > 50) {
    data_tops <- data.frame(var = tops, .rank = 1:length(tops))
    groups <- grouping_vec2list(data_tops$.rank, 50, T)
    groups <- rep(paste0("Rank ", 1:length(groups)), lengths(groups))
    data_tops$.groups <- groups
    data <- tbmerge(data, data_tops, by = "var", all.x = T)
    data <- mutate(data, var = stringr::str_trunc(var, 30))
  }
  p <- ggplot(data, aes(x = reorder(var, value), y = value)) +
    geom_col(aes(fill = group), position = "stack", width = .6) +
    coord_flip() +
    scale_fill_manual(values = ggsci::pal_npg()(8)) +
    labs(x = lab_x, y = lab_y, fill = lab_fill)+
    theme_minimal()
  if (top > 50) {
    p <- p + facet_wrap(~ .groups, scales = "free")
    p <- wrap(p, 2 + length(unique(data$.groups)) * 5.5, 9.5)
  } else {
    p <- wrap(p, 6, 3 + length(tops) * .2)
  }
  attr(p, "tops") <- tops
  p
}

setMethod("map", signature = c(x = "job_lasso"),
  function(x, ref, pvalue = T){
    data <- select(object(x), dplyr::all_of(ref))
    data$group <- x$metadata$group
    data <- tidyr::gather(data, var, value, -group)
    p <- .map_boxplot(data, pvalue)
    p
  })

.map_boxplot <- function(data, pvalue) {
  p <- ggplot(data, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(stroke = .1, shape = 21, width = .1) +
    facet_wrap(~ var) +
    ggsci::scale_fill_npg() +
    labs(x = "Group", y = "Value") +
    theme(legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_blank()
  if (pvalue) {
    fn <- fivenum(data$value)
    levels <- 2 * seq(length(unique(data$group)) - 1)
    hs <- fn[4] + abs(fn[5] - fn[1]) / levels
    p <- ggpval::add_pval(p, heights = hs, pval_text_adj = (fn[5] - fn[1]) / 15)
  }
  p
}

.map_boxplot2 <- function(data, pvalue, x = "group") {
  p <- ggplot(data, aes(x = !!rlang::sym(x), y = value, color = group)) +
    geom_boxplot(outlier.shape = NA, fill = "transparent") +
    geom_jitter(aes(x = !!rlang::sym(x), y = value, fill = group),
      stroke = 0, shape = 21, width = .1, color = "transparent") +
    facet_wrap(~ var) +
    scale_fill_manual(values = color_set()) +
    scale_color_manual(values = color_set()) +
    labs(x = "Group", y = "Value") +
    theme_light() +
    theme(legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_blank()
  if (pvalue) {
    fn <- fivenum(data$value)
    levels <- 2 * seq(length(unique(data$group)) - 1)
    hs <- fn[4] + abs(fn[5] - fn[1]) / levels
    p <- ggpval::add_pval(p, heights = hs, pval_text_adj = (fn[5] - fn[1]) / 15)
  }
  p
}

