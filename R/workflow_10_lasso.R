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
    info = c("...")
    ))

setGeneric("asjob_lasso", 
  function(x, ...) standardGeneric("asjob_lasso"))

setMethod("asjob_lasso", signature = c(x = "job_limma"),
  function(x, use.filter, use = "gene_name"){
    step_message("The default, 'job_limma' from 'job_tcga' were adapted to
      convertion.
      "
    )
    if (!is(object(x), 'EList'))
      stop("is(object(x), 'EList') == F")
    pos <- object(x)$genes[[use]] %in% use.filter
    object <- e(limma::`[.EList`(object(x), pos, ))
    mtx <- t(object$E)
    colnames(mtx) <- object$genes[[ use ]]
    x <- .job_lasso(object = mtx)
    x@params$metadata <- as_tibble(object$targets)
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

plot_colStack.efs <- function(raw, top = 30, lab_x = "Variable", lab_y = "Importance", lab_fill = "Method")
{
  data <- as_tibble(raw)
  data <- rename(data, group = rownames)
  data <- tidyr::gather(data, var, value, -group)
  tops <- sort(vapply(split(data$value, data$var), sum, double(1)), decreasing = T)
  tops <- head(names(tops), top)
  data <- filter(data, var %in% dplyr::all_of(tops))
  p <- ggplot(data, aes(x = reorder(var, value), y = value)) +
    geom_col(aes(fill = group), position = "stack", width = .6) +
    coord_flip() +
    scale_fill_manual(values = ggsci::pal_npg()(8)) +
    labs(x = lab_x, y = lab_y, fill = lab_fill)+
    theme_minimal()
  p <- wrap(p, 6, 3 + length(tops) * .2)
  attr(p, "tops") <- tops
  p
}
