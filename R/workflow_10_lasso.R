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
    analysis = "COX 回归"
    ))

setGeneric("asjob_lasso", 
  function(x, ...) standardGeneric("asjob_lasso"))

setMethod("asjob_lasso", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = "gene_name", from_normed = T,
    fun_scale = function(x) log2(x + 1), dup_method = c("max", "min", "mean"))
  {
    step_message("The default, 'job_limma' from 'job_tcga' were adapted to convertion.")
    if (from_normed && x@step >= 1L) {
      if (!is(x@params$normed_data, 'EList'))
        stop("is(x@params$normed_data, 'EList') == F")
      if (!is.null(use.filter)) {
        pos <- x@params$normed_data$genes[[use]] %in% use.filter
        object <- e(limma::`[.EList`(x@params$normed_data, pos, ))
        if (any(duplicated(object$genes[[ use ]]))) {
          cli::cli_alert_warning("Duplicated names (genes) founds.")
          lst <- .deduplicated_genes(object$E, object$genes, use, dup_method)
          object$E <- lst$counts
          object$genes <- lst$genes
        }
      } else {
        object <- x@params$normed_data
      }
      message("Ignore parameter of `fun_scale`")
      mtx <- t(object$E)
      colnames(mtx) <- object$genes[[ use ]]
      x <- .job_lasso(object = mtx)
      x@params$metadata <- as_tibble(object$targets)
    } else {
      if (x@step == 0L) {
        message(glue::glue("x@step == 0L, use raw counts."))
      }
      genes <- object(x)$genes
      counts <- object(x)$counts
      if (!is.null(use.filter)) {
        pos <- genes[[ use ]] %in% use.filter
        genes <- genes[pos, ]
        counts <- counts[pos, ]
        if (any(duplicated(genes[[ use ]]))) {
          cli::cli_alert_warning("Duplicated names (genes) founds.")
          lst <- .deduplicated_genes(counts, genes, use, dup_method)
          genes <- lst$genes
          counts <- lst$counts
        }
      }
      mtx <- t(counts)
      if (!is.null(fun_scale)) {
        message("Use `fun_scale` to scale the data.")
        mtx <- fun_scale(mtx)
        attr(mtx, "scaled:center") <- attr(mtx, "scaled:scale") <- NULL
      }
      colnames(mtx) <- genes[[ use ]]
      metadata <- as_tibble(object(x)$samples)
      x <- .job_lasso(object = mtx)
      x@params$metadata <- metadata
    }
    return(x)
  })

.deduplicated_genes <- function(counts, genes, use, dup_method) {
  counts <- data.frame(counts, check.names = F)
  counts <- dplyr::mutate(counts, gene = !!genes[[ use ]])
  counts <- dplyr::group_by(counts, gene)
  message(glue::glue("Use function `{dup_method}` for deduplicated."))
  dup_method <- match.fun(dup_method)
  counts <- dplyr::summarise(counts, dplyr::across(dplyr::everything(), ~ dup_method(.x)))
  counts <- dplyr::arrange(counts, gene)
  counts <- data.frame(counts[, -1], check.names = F)
  genes <- genes[!duplicated(genes[[ use ]]), ]
  genes <- dplyr::arrange(genes, !!rlang::sym(use))
  rownames(counts) <- rownames(genes)
  namel(counts, genes)
}

job_lasso <- function(data, metadata) {
  .job_lasso(object = data, params = list(metadata = metadata))
}

setMethod("step0", signature = c(x = "job_lasso"),
  function(x){
    step_message("Prepare your data with function `asjob_lasso`. ")
  })

setMethod("step1", signature = c(x = "job_lasso"),
  function(x, target = x@params$metadata$group, levels = c("Alive", "Dead"),
    time = x@params$metadata$days_to_last_follow_up, n.train = .8, seed = 555)
  {
    step_message("Data preparation.")
    len <- nrow(object(x))
    all <- 1:len
    set.seed(seed)
    wh.train <- sample(all, round(len * n.train))
    wh.valid <- all[!all %in% wh.train]
    x$train <- object(x)[wh.train, ]
    x$valid <- object(x)[wh.valid, ]
    x$train.target <- target[wh.train]
    x$valid.target <- target[wh.valid]
    x$train.time <- time[wh.train]
    x$valid.time <- time[wh.valid]
    x$time <- time
    x$target <- target
    x$levels <- levels
    x$wh.train <- wh.train
    p.all <- gtitle(new_pie(target)@data, "All")
    p.train <- gtitle(new_pie(x$train.target)@data, "Train")
    p.valid <- gtitle(new_pie(x$valid.target)@data, "Validate")
    p.consist <- wrap(xf(p.all = 1, p.train = 1, p.valid = 1), 6, 2)
    x@plots[[ 1 ]] <- namel(p.consist)
    return(x)
  })

setMethod("step2", signature = c(x = "job_lasso"),
  function(x, top = 30, efs = T, use_data = c("all", "train")){
    step_message("Evaluate variable (genes) importance.")
    if (efs) {
      use_data <- match.arg(use_data)
      if (use_data == "train") {
        target <- x@params$train.target
        target <- ifelse(target == x$levels[1], 0L, 1L)
        data <- data.frame(cbind(x$train, target))  
      } else {
        target <- x$target
        target <- ifelse(target == x$levels[1], 0L, 1L)
        data <- data.frame(cbind(object(x), target))  
      }
      efs <- e(EFS::ensemble_fs(data, ncol(data), runs = 10))
      colnames(efs) <- colnames(object(x))
      p.TopFeaturesSelectedByEFS <- plot_colStack.efs(efs, top = top)
      tops <- attr(p.TopFeaturesSelectedByEFS, "tops")
      x$efs <- efs
      x$efs_tops <- tops
      x@tables[[ 2 ]] <- namel(efs)
      x <- plotsAdd(x, p.TopFeaturesSelectedByEFS)
      x <- methodAdd(x, "以 R 包 `EFS` ({packageVersion('EFS')}) {cite_show('EfsAnEnsemblNeuman2017')} 筛选关键基因。")
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_lasso"),
  function(x, family = "cox", nfold = 10, use_tops = F, use_data = c("all", "train"),
    uni_cox = T, multi_cox = F, multi_cox.inherits = T,
    fun_multiCox = c("coxph", "glmnet", "cv.glmnet"),
    alpha = 1, ...)
  {
    if (use_tops && !is.null(x$efs_tops)) {
      fun <- function(mtx) mtx[, colnames(mtx) %in% x$efs_tops ]
      object(x) <- fun(object(x))
      x$train <- fun(x$train)
      x$valid <- fun(x$valid)
    }
    x$use_data <- use_data <- match.arg(use_data)
    if (use_data == "train") {
      data <- x$train
      target <- x$train.target
      time <- x$train.time
      valid <- x$valid
      valid.target <- x$valid.target
    } else if (use_data == "all") {
      valid <- data <- object(x)
      valid.target <- target <- x$target
      time <- x$time
    }
    if (family == "cox") {
      message(glue::glue("Use family of {family}, switch target to: 0 or 1."))
      target <- ifelse(target == x$levels[1], 0L, 1L)
      target <- survival::Surv(time = time, event = target)
    }
    if (family == "cox" && uni_cox) {
      # Univariate lasso
      cli::cli_alert_info("survival::coxph")
      uni_cox <- pbapply::pblapply(colnames(data),
        function(feature) {
          res <- survival::coxph(as.formula(paste0("target ~ `", feature, "`")),
              data.frame(data, check.names = F))
          summary(res)$coefficients[1, ]
        })
      uni_cox <- do.call(dplyr::bind_rows, uni_cox)
      x$uni_cox <- uni_cox <- dplyr::mutate(uni_cox, feature = !!colnames(data), .before = 1)
      sig.uni_cox <- dplyr::mutate(uni_cox, p.adjust = stats::p.adjust(`Pr(>|z|)`, "fdr"))
      sig.uni_cox <- dplyr::filter(sig.uni_cox, `Pr(>|z|)` < .05)
      x$sig.uni_cox <- dplyr::rename(sig.uni_cox, pvalue = `Pr(>|z|)`)
      x <- tablesAdd(x, t.sigUnivariateCoxCoefficients = x$sig.uni_cox)
      if (!is.null(x$efs_tops)) {
        x$sig.cox_efs <- dplyr::filter(sig.uni_cox, feature %in% x$efs_tops)
        x$sig.cox_efs <- .set_lab(x$sig.cox_efs, sig(x), "Significant Univariate COX results of EFS tops")
      }
      meth(x)$step3 <- glue::glue("以 R 包 `survival` ({packageVersion('survival')}) 进行单因素 COX 回归 (`survival::coxph`)。筛选 `Pr(>|z|)` < .05` 的基因。")
    }
    if (multi_cox) {
      if (multi_cox.inherits) {
        cli::cli_alert_info("Inherits results (significant feature) from Univariate COX.")
        data <- data[, colnames(data) %in% sig.uni_cox$feature ]
        valid <- valid[, colnames(valid) %in% sig.uni_cox$feature ]
        text_inherits <- "在单因素回归得到的基因的基础上，"
      } else {
        text_inherits <- ""
      }
      multi_cox <- list()
      # Multivariate
      fun_multiCox <- match.arg(fun_multiCox)
      if (fun_multiCox == "coxph") {
        res <- survival::coxph(target ~ ., data.frame(data, check.names = F))
        multi_cox$model <- res
        multi_cox$coef <- dplyr::rename(as_tibble(summary(res)$coefficients), feature = 1)
        x <- tablesAdd(x, t.sigMultivariateCoxCoefficients = multi_cox$coef)
        x <- methodAdd(x, "{text_inherits} 以 R 包 `survival` ({packageVersion('survival')}) 做多因素 COX 回归 (`survival::coxph`)。", T)
      } else if (fun_multiCox == "cv.glmnet") {
        methodName <- if (alpha) "lasso" else "ridge"
        multi_cox$model <- e(glmnet::cv.glmnet(data, target,
            alpha = alpha, family = family, nfold = nfold, ...
            ))  
        model <- multi_cox$model
        p.lassoCOX_model <- wrap(as_grob(expression(plot(model)), environment()), 6, 6)
        levels <- x$levels
        types.coef <- c("lambda.min", "lambda.1se")
        res <- sapply(types.coef, simplify = F,
          function(lam) {
            preds <- e(stats::predict(model, newx = valid, s = model[[ lam ]]))
            roc <- e(pROC::roc(valid.target, as.numeric(preds), plot = F, levels = levels))
            namel(preds, roc)
          })
        multi_cox$preds <- lapply(res, function(lst) lst$preds)
        multi_cox$roc <- lapply(res, function(lst) lst$roc)
        if (T) {
          message("Got coefficients.")
          x$lambda <- c(min = model$lambda.min, `1se` = model$lambda.1se)
          multi_cox$coef <- coefficients(model, s = c(lambda.min = model$lambda.min,
              lambda.1se = model$lambda.1se))
          if (all(!multi_cox$coef[, 2])) {
            message(crayon::red("lambda.1se has non coeffients, use Lambda.min."))
            sig.mul_cox <- dplyr::select(
              as_tibble(as.matrix(multi_cox$coef)),
              feature = 1, coef = 2
            )
            sig.mul_cox <- dplyr::filter(sig.mul_cox, coef != 0)
          } else {
            sig.mul_cox <- dplyr::rename(
              as_tibble(as.matrix(multi_cox$coef)),
              feature = 1, coef.lambda.min = 2, coef.lambda.1se = 3
            )
          }
          x <- tablesAdd(x, t.sigMultivariateCoxCoefficients = sig.mul_cox)
        }
        if (T) {
          n <- 0L
          p.lassoCOX_ROC <- lapply(multi_cox$roc,
            function(roc) {
              n <<- n + 1L
              plot_roc(roc)
            })
          names(p.lassoCOX_ROC) <- types.coef
        }
        if (T) {
          coef <- as_tibble(Matrix::as.matrix(multi_cox$coef))
          colnames(coef)[-1] <- types.coef
          coef <- arrange(coef, dplyr::desc(abs(lambda.min)))
          coef <- rename(coef, variable = rownames)
          message("Plot Significant coefficients.")
          p.lassoCOX_coeffients <- plot_sig(filter(coef, !grepl("Inter", variable)))
        }
        x <- plotsAdd(x, p.lassoCOX_model, p.lassoCOX_ROC, p.lassoCOX_coeffients)
        x <- methodAdd(x, "{text_inherits} 以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归，以 `{fun_multiCox}` 函数作 {nfold} 交叉验证获得模型。", T)
      } else if (fun_multiCox == "glmnet") {
        methodName <- if (alpha) "lasso" else "ridge"
        multi_cox$model <- e(glmnet::glmnet(data, target,
            alpha = alpha, family = family, ...))  
        multi_cox$coef <- coefficients(multi_cox$model)
        x <- methodAdd(x, "{text_inherits} 以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归。", T)
      }
      x$multi_cox <- multi_cox
    }
    if (!is.null(x$uni_cox) && !is.null(x$multi_cox$coef) && fun_multiCox == "coxph") {
      message("Try merge results.")
      if (nrow(x$uni_cox) == nrow(x$multi_cox$coef)) {
        message("The same row number of Univariate COX and Multivariate COX, try merge.")
        coefMerge <- mapply(list(x$uni_cox, x$multi_cox$coef), c("Uni_", "Multi_"), SIMPLIFY = F,
          FUN = function(data, name) {
            data <- dplyr::select(data, coefficients = coef, p = `Pr(>|z|)`)
            colnames(data) <- paste0(name, colnames(data))
            data
          })
        t.CoefficientsOfCOX <- dplyr::bind_cols(dplyr::select(x$uni_cox, feature),
          coefMerge[[1]], coefMerge[[2]])
        x <- tablesAdd(x, t.CoefficientsOfCOX)
      }
    }
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

.map_boxplot2 <- function(data, pvalue, x = "group", y = "value",
  xlab = "Group", ylab = "Value", ids = "var", test = "wilcox.test", ...)
{
  p <- ggplot(data, aes(x = !!rlang::sym(x), y = !!rlang::sym(y), color = !!rlang::sym(x))) +
    geom_boxplot(outlier.shape = NA, fill = "transparent") +
    geom_jitter(aes(x = !!rlang::sym(x), y = !!rlang::sym(y), fill = !!rlang::sym(x)),
      stroke = 0, shape = 21, width = .1, color = "transparent") +
    facet_wrap(ggplot2::vars(!!rlang::sym(ids)), ...) +
    scale_fill_manual(values = color_set()) +
    scale_color_manual(values = color_set()) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_blank()
  if (pvalue) {
    fn <- fivenum(data[[ y ]])
    levels <- 2 * seq(length(unique(data[[ x ]])) - 1)
    hs <- fn[4] + abs(fn[5] - fn[1]) / levels
    p <- ggpval::add_pval(p, heights = hs, pval_text_adj = (fn[5] - fn[1]) / 15, test = test)
  }
  p
}

