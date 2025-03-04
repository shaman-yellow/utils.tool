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

setGeneric("asjob_lasso", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_lasso"))

setMethod("asjob_lasso", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = .guess_symbol(x), from_normed = TRUE,
    fun_scale = function(x) scale(x, TRUE, TRUE),
    dup_method = c("max", "min", "mean"), 
    use.format = TRUE, ...)
  {
    step_message("The default, 'job_limma' from 'job_tcga' were adapted to convertion.")
    x <- dedup_by_rank.job_limma(x, ref.use = use, ...)
    dup_method <- match.arg(dup_method)
    project <- x$project
    if (x$normed) {
      x$normed_data <- new("EList", list(E = object(x), targets = x$metadata, genes = x$genes))
      from_normed <- TRUE
      x@step <- 1L
    }
    if (from_normed && x@step >= 1L) {
      snap <- ""
      if (!is(x@params$normed_data, 'EList'))
        stop("is(x@params$normed_data, 'EList') == FALSE")
      if (!is.null(use.filter)) {
        if (is(use.filter, "feature")) {
          message("Detected feature input.")
          snap <- glue::glue("将{snap(use.filter)}用于模型建立。")
          use.filter <- unlist(use.filter@.Data, use.names = FALSE)
        }
        if (use.format) {
          pos <- gname(x@params$normed_data$genes[[use]]) %in% use.filter
        } else {
          pos <- x@params$normed_data$genes[[use]] %in% use.filter
        }
        object <- e(limma::`[.EList`(x@params$normed_data, pos, ))
        if (any(duplicated(object$genes[[ use ]]))) {
          cli::cli_alert_warning("Duplicated names (genes) founds.")
          lst <- .deduplicated_genes(object$E, object$genes, use, dup_method)
          object$E <- lst$counts
          object$genes <- lst$genes
        }
        message(glue::glue("After filtered, dim: {bind(dim(object$E))}"))
      } else {
        object <- x@params$normed_data
      }
      mtx <- t(object$E)
      gnames <- object$genes[[ use ]]
      if (use.format) {
        colnames(mtx) <- gname(gnames)
      } else {
        colnames(mtx) <- gnames
      }
      mtx <- mtx[, !duplicated(colnames(mtx)) & !is.na(colnames(mtx))]
      x <- .job_lasso(object = fun_scale(mtx))
      x <- snapAdd(x, snap)
      x <- snapAdd(x, "共 {ncol(mtx)} 个基因在数据集 {project} 中找到 (根据基因名匹配)。")
      x@params$metadata <- as_tibble(object$targets)
    } else {
      stop("x@step == 0L, case deprecated.")
      # if (x@step == 0L) {
      #   message(glue::glue("x@step == 0L, use raw counts."))
      # }
      # genes <- object(x)$genes
      # counts <- object(x)$counts
      # if (!is.null(use.filter)) {
      #   pos <- genes[[ use ]] %in% use.filter
      #   genes <- genes[pos, ]
      #   counts <- counts[pos, ]
      #   if (any(duplicated(genes[[ use ]]))) {
      #     cli::cli_alert_warning("Duplicated names (genes) founds.")
      #     lst <- .deduplicated_genes(counts, genes, use, dup_method)
      #     genes <- lst$genes
      #     counts <- lst$counts
      #   }
      # }
      # mtx <- t(counts)
      # if (!is.null(fun_scale)) {
      #   message("Use `fun_scale` to scale the data.")
      #   mtx <- fun_scale(mtx)
      #   attr(mtx, "scaled:center") <- attr(mtx, "scaled:scale") <- NULL
      # }
      # colnames(mtx) <- genes[[ use ]]
      # metadata <- as_tibble(object(x)$samples)
      # x <- .job_lasso(object = mtx)
      # x@params$metadata <- metadata
    }
    return(x)
  })

.deduplicated_genes <- function(counts, genes, use, dup_method) {
  counts <- data.frame(counts, check.names = FALSE)
  counts <- dplyr::mutate(counts, gene = !!genes[[ use ]])
  counts <- dplyr::group_by(counts, gene)
  message(glue::glue("Use function `{dup_method}` for deduplicated."))
  dup_method <- match.fun(dup_method)
  counts <- dplyr::summarise(counts, dplyr::across(dplyr::everything(), ~ dup_method(.x)))
  counts <- dplyr::arrange(counts, gene)
  counts <- data.frame(counts[, -1], check.names = FALSE)
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
  function(x, target = "vital_status", levels = c("Alive", "Dead"),
    time = "days_to_last_follow_up", n.train = .8, seed = 555)
  {
    step_message("Data preparation.")
    set.seed(seed)
    if (!is.null(time)) {
      x$metadata[[ time ]] <- .mutate_last_follow_up_time(x$metadata[[ time ]])
    }
    if (TRUE) {
      keep <- (x$metadata[[ target ]] %in% levels)
      if (!is.null(time)) {
        keep <- keep & (x$metadata[[ time ]] > 0)
      }
      x$metadata <- x$metadata[keep, ]
      object(x) <- object(x)[keep, ]
    }
    num <- nrow(object(x))
    x$n.train
    x$wh.train <- sample(1:num, round(num * n.train))
    x$wh.valid <- setdiff(1:num, x$wh.train)
    x$target <- target
    x$time <- time
    x$levels <- levels
    x$get <- function(wh = c("train", "valid", "all"), x = get("x", envir = parent.frame(1)))
    {
      wh <- match.arg(wh)
      wh <- switch(wh, train = x$wh.train, valid = x$wh.valid, all = c(x$wh.train, x$wh.valid))
      data <- object(x)[wh, ]
      types <- x$metadata[[ x$target ]][ wh ]
      event <- ifelse(types == x$levels[1], 0L, 1L)
      if (!is.null(time)) {
        time <- x$metadata[[ time ]][ wh ]
        surv <- survival::Surv(time = time, event = event)
      } else {
        time <- surv <- NULL
      }
      namel(data, event, types, time, surv)
    }
    p.all <- gtitle(new_pie(x$get("a")$types)@data, "All")
    p.train <- gtitle(new_pie(x$get("t")$types)@data, "Train")
    p.valid <- gtitle(new_pie(x$get("v")$types)@data, "Validate")
    p.consist <- wrap(xf(p.all = 1, p.train = 1, p.valid = 1), 6, 2)
    x <- plotsAdd(x, p.consist)
    target <- x$get("a")$types
    x <- snapAdd(x, "所有数据生存状态 (去除生存状态未知的数据)，({try_snap(target)})。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_lasso"),
  function(x, use_data = c("all", "train"), top = 30, efs = FALSE)
  {
    step_message("Evaluate variable (genes) importance.")
    use_data <- match.arg(use_data)
    if (use_data == "train") {
      x <- snapAdd(x, "使用随机提取的 train 数据集 (train:valid = {x$n.train / (1-x$n.train)}:1)。")
    }
    x$use_data <- use_data
    if (efs) {
      data <- data.frame(cbind(x$get(use_data)$data, x$get(use_data)$event))  
      efs <- e(EFS::ensemble_fs(data, ncol(data), runs = 10))
      colnames(efs) <- colnames(object(x))
      p.TopFeaturesSelectedByEFS <- plot_colStack.efs(efs, top = top)
      tops <- attr(p.TopFeaturesSelectedByEFS, "tops")
      x$efs <- efs
      x$efs_tops <- tops
      x <- tablesAdd(x, t.efs = efs)
      x <- plotsAdd(x, p.TopFeaturesSelectedByEFS)
      x <- methodAdd(x, "以 R 包 `EFS` ({packageVersion('EFS')}) {cite_show('EfsAnEnsemblNeuman2017')} 筛选关键基因。")
      x <- snapAdd(x, "使用 EFS 对基因集筛选 TOP {top} 的基因。")
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_lasso"),
  function(x, use_data = x$use_data, family = "cox", use_tops = FALSE, cutoff = .05)
  {
    step_message("Univariate COX.")
    lst <- x$get(use_data)
    data <- lst$data
    if (family == "cox") {
      message(glue::glue("Use family of {family}"))
      target <- lst$surv
    } else {
      target <- lst$event
    }
    if (use_tops && !is.null(x$efs_tops)) {
      fun <- function(mtx) mtx[, colnames(mtx) %in% x$efs_tops ]
      data <- fun(data)
      x <- snapAdd(x, "使用 EFS 筛选的 TOP 基因，")
    }
    # Univariate
    cli::cli_alert_info("survival::coxph")
    uni_cox <- pbapply::pblapply(colnames(data),
      function(feature) {
        res <- survival::coxph(as.formula(paste0("target ~ `", feature, "`")),
          data.frame(data, check.names = FALSE))
        summary(res)$coefficients[1, ]
      })
    uni_cox <- do.call(dplyr::bind_rows, uni_cox)
    x$uni_cox <- uni_cox <- dplyr::mutate(uni_cox, feature = !!colnames(data), .before = 1)
    sig.uni_cox <- dplyr::mutate(uni_cox, p.adjust = stats::p.adjust(`Pr(>|z|)`, "fdr"))
    sig.uni_cox <- dplyr::filter(sig.uni_cox, `Pr(>|z|)` < cutoff)
    x$sig.uni_cox <- dplyr::rename(sig.uni_cox, pvalue = `Pr(>|z|)`)
    feature(x) <- x$sig.uni_cox$feature
    x <- tablesAdd(x, t.sigUnivariateCoxCoefficients = x$sig.uni_cox)
    x <- snapAdd(x, "执行单因素 COX 回归，筛选 P 值 &lt; {cutoff}，共筛选到 {nrow(x$sig.uni_cox)} 个基因。")
    if (!is.null(x$efs_tops)) {
      t.significantUnivariateCoxIntersectWithEFS_tops <- dplyr::filter(sig.uni_cox, feature %in% x$efs_tops)
      x <- tablesAdd(x, t.significantUnivariateCoxIntersectWithEFS_tops)
      x <- snapAdd(x, "EFS 的 TOP 基因与单因素 COX 回归的显著基因交集共 {nrow(x$sig.cox_efs)} 个。")
    }
    x <- methodAdd(x, "以 R 包 `survival` ({packageVersion('survival')}) 进行单因素 COX 回归 (`survival::coxph`)。筛选 `Pr(>|z|)` < .05` 的基因。")
    return(x)
  })

setMethod("step4", signature = c(x = "job_lasso"),
  function(x, use_data = x$use_data, use_valid = use_data, inherit_unicox = TRUE,
    inherit_unicox.cut.p = NULL,
    fun = c("coxph", "glmnet", "cv.glmnet"), nfold = 10,
    alpha = 1, family = "cox", type.measure = c(
      "deviance", "C", "default"
    ), ...)
  {
    step_message("Multivariate COX.")
    fun_multiCox <- match.arg(fun)
    type.measure <- match.arg(type.measure)
    data_lst <- x$get(use_data)
    data <- data_lst$data
    valid_lst <- x$get(use_data)
    valid <- valid_lst$data
    x$nfold <- nfold
    if (family == "cox") {
      target <- data_lst$surv
    } else {
      target <- data_lst$event
    }
    if (inherit_unicox) {
      cli::cli_alert_info("Inherits results (significant feature) from Univariate COX.")
      sig.uni_cox <- x$sig.uni_cox
      if (!is.null(inherit_unicox.cut.p)) {
        sig.uni_cox <- dplyr::filter(sig.uni_cox, pvalue < inherit_unicox.cut.p)
        x <- snapAdd(x, "在单因素回归得到的基因 (P &lt; {inherit_unicox.cut.p}) 的基础上，")
      } else {
        x <- snapAdd(x, "在单因素回归得到的基因 (P &lt; 0.05) 的基础上，")
      }
      data <- data[, colnames(data) %in% sig.uni_cox$feature ]
      valid <- valid[, colnames(valid) %in% sig.uni_cox$feature ]
    }
    multi_cox <- list()
    # Multivariate
    if (fun_multiCox == "coxph") {
      res <- survival::coxph(target ~ ., data.frame(data, check.names = FALSE))
      multi_cox$model <- res
      multi_cox$coef <- dplyr::rename(as_tibble(summary(res)$coefficients), feature = 1)
      x <- tablesAdd(x, t.sigMultivariateCoxCoefficients = multi_cox$coef)
      x <- methodAdd(x, "以 R 包 `survival` ({packageVersion('survival')}) 做多因素 COX 回归 (`survival::coxph`)。", TRUE)
      x <- snapAdd(x, "执行多因素 COX 回归 (`survival::coxph`)。")
    } else if (fun_multiCox == "cv.glmnet") {
      methodName <- if (alpha == 1) "lasso"
        else if (!alpha) "ridge"
        else if (alpha > 0 && alpha < 1) "Elastic Net"
      set.seed(x$seed)
      multi_cox$model <- model <- e(glmnet::cv.glmnet(data, target,
          alpha = alpha, family = family, nfold = nfold, type.measure = type.measure, ...))
      p.lassoCOX_model <- wrap(as_grob(expression(plot(model)), environment()), 6, 6, showtext = TRUE)
      lambdas <- c("lambda.min", "lambda.1se")
      res <- sapply(lambdas, simplify = FALSE,
        function(lam) {
          preds <- e(stats::predict(model, newx = valid, s = model[[ lam ]]))
          roc <- e(pROC::roc(valid_lst$types, as.numeric(preds), plot = FALSE, levels = x$levels))
          namel(preds, roc)
        })
      multi_cox$preds <- lapply(res, function(lst) lst$preds)
      multi_cox$roc <- lapply(res, function(lst) lst$roc)
      if (TRUE) {
        message("Got coefficients.")
        x$lambda <- c(min = model$lambda.min, `1se` = model$lambda.1se)
        x <- snapAdd(x, "使用 `glmnet::cv.glmnet` 作 {nfold} 倍交叉验证 (评估方式为 {model$name})，筛选 lambda 值。lambda.min, lambda.1se 值分别为 {bind(round(x$lambda, 3))} (R 随机种子为 {x$seed})。")
        multi_cox$coef <- coefficients(
          model, s = c(lambda.min = model$lambda.min, lambda.1se = model$lambda.1se)
        )
        x$sig.mul_cox <- sig.mul_cox <- dplyr::rename(
          as_tibble(as.matrix(multi_cox$coef[-1, ])),
          feature = 1, coef.lambda.min = 2, coef.lambda.1se = 3
        )
        if (all(!multi_cox$coef[-1, 2])) {
          message(crayon::red("lambda.1se has non coeffients, use Lambda.min (-> coef)."))
        }
        s.com <- vapply(1:2,
          function(x) {
            n <- multi_cox$coef[-1, x]
            length(n[ n != 0 ])
          }, integer(1))
        x$nfeature_lambdas <- s.com
        x <- snapAdd(x, "对应的特征数 (基因数) 分别为 {bind(s.com)}。")
        x <- tablesAdd(x, t.sigMultivariateCoxCoefficients = sig.mul_cox)
      }
      if (TRUE) {
        p.lassoCOX_ROC <- lapply(multi_cox$roc,
          function(roc) {
            plot_roc(roc)
          })
        names(p.lassoCOX_ROC) <- lambdas
      }
      if (TRUE) {
        coef <- as_tibble(Matrix::as.matrix(multi_cox$coef))
        colnames(coef)[-1] <- lambdas
        coef <- dplyr::rename(coef, variable = rownames)
        message("Plot Significant coefficients.")
        p.lassoCOX_coeffients <- sapply(
          lambdas, simplify = FALSE,
          function(lam) {
            coef <- dplyr::arrange(coef, dplyr::desc(abs(!!rlang::sym(lam))))
            plot_sig(dplyr::filter(coef, !grepl("Inter", variable)), y = lam)
          }
        )
      }
      x <- plotsAdd(x, p.lassoCOX_model, p.lassoCOX_ROC, p.lassoCOX_coeffients)
      x <- methodAdd(x, "以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归，以 `{fun_multiCox}` 函数作 {nfold} 交叉验证获得模型。", TRUE)
    } else if (fun_multiCox == "glmnet") {
      methodName <- if (alpha) "lasso" else "ridge"
      multi_cox$model <- e(glmnet::glmnet(data, target,
          alpha = alpha, family = family, ...))  
      multi_cox$coef <- coefficients(multi_cox$model)
      x <- methodAdd(x, "以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归。", TRUE)
    }
    x$multi_cox <- multi_cox
    x$fun_multiCox <- fun_multiCox
    return(x)
  })

setMethod("step5", signature = c(x = "job_lasso"),
  function(x){
    step_message("Try Merge Univariate cox and Multivariate cox results (coxph).")
    if (!is.null(x$uni_cox) && !is.null(x$multi_cox$coef) && x$fun_multiCox == "coxph") {
      message("Try merge results.")
      if (nrow(x$uni_cox) == nrow(x$multi_cox$coef)) {
        message("The same row number of Univariate COX and Multivariate COX, try merge.")
        coefMerge <- mapply(list(x$uni_cox, x$multi_cox$coef), c("Uni_", "Multi_"), SIMPLIFY = FALSE,
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

plot_roc <- function(roc, cols = ggsci::pal_npg()(9), groups = NULL) {
  if (is(roc, "list")) {
    pos <- seq(.05, .3, length.out = length(roc))
    if (!is.null(groups)) {
      if (!all(names(roc) %in% names(groups))) {
        stop('!all(names(roc) %in% names(groups)).')
      }
      names(roc) <- paste0(
        vapply(names(roc), function(name) groups[[ name ]], character(1)),
        "_", names(roc)
      )
      roc <- sort_by(roc, names(roc))
    }
    expr <- expression({
      lapply(seq_along(roc), 
        function(n) {
          plot(roc[[n]], col = cols[n], add = n > 1)
          graphics::text(
            .7, pos[n], paste0(names(roc)[n], " AUC: ", signif(roc[[n]]$auc[[1]], 2)),
            cex = 1.2, col = cols[n], adj = c(0, 0)
          )
        })
    })
  } else {
    expr <- expression({
      plot(roc)
      text(.1, .1, paste0("AUC: ", round(roc$auc[[1]], 2)), cex = 1.2)
    })
  }
  wrap(as_grob(expr, environment()), 6, 6)
}

plot_sig <- function(data, x = colnames(data)[1], y = colnames(data)[2],
  lab_x = "Variables", lab_y = "Coefficients") 
{
  data <- dplyr::filter(data, abs(!!rlang::sym(y)) != 0)
  p <- ggplot(data, aes(x = stats::reorder(!!rlang::sym(x),
        !!rlang::sym(y)), y = !!rlang::sym(y))) +
    geom_point(aes(color = !!rlang::sym(y))) +
    labs(x = lab_x, y = lab_y) +
    theme_minimal() +
    coord_flip() +
    guides(color = "none") +
    theme()
  wrap(p, 6, 1 + nrow(data) * .15)
}

plot_colStack.efs <- function(raw, top = 30,
  lab_x = "Variable", lab_y = "Importance", lab_fill = "Method")
{
  data <- as_tibble(raw)
  data <- rename(data, group = rownames)
  data <- tidyr::gather(data, var, value, -group)
  tops <- sort(vapply(split(data$value, data$var), sum, double(1)), decreasing = TRUE)
  tops <- head(names(tops), top)
  data <- filter(data, var %in% dplyr::all_of(tops))
  if (top > 50) {
    data_tops <- data.frame(var = tops, .rank = seq_along(tops))
    groups <- grouping_vec2list(data_tops$.rank, 50, TRUE)
    groups <- rep(paste0("Rank ", seq_along(groups)), lengths(groups))
    data_tops$.groups <- groups
    data <- tbmerge(data, data_tops, by = "var", all.x = TRUE)
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
  function(x, ref, pvalue = TRUE){
    data <- select(object(x), dplyr::all_of(ref))
    data$group <- x$metadata$group
    data <- tidyr::gather(data, var, value, -group)
    p <- .map_boxplot(data, pvalue)
    p
  })

setMethod("merge", signature = c(x = "job_lasso", y = "job_lasso"),
  function(x, y, scale = TRUE, ...)
  {
    common <- intersect(colnames(x$metadata), colnames(y$metadata))
    message(glue::glue("merge metadata of columns: {bind(common)}"))
    lst_metas <- lapply(list(x$metadata, y$metadata),
      function(meta) {
        dplyr::select(meta, dplyr::all_of(common))
      })
    sigs <- names(lst_metas) <- c(x@sig, y@sig)
    metadata <- frbind(lst_metas, idcol = "Dataset")
    gcommon <- intersect(colnames(object(x)), colnames(object(y)))
    object <- rbind(object(x)[, gcommon ], object(y)[, gcommon ])
    object <- scale(object)
    x <- .job_lasso(object = object, params = list(metadata = metadata))
    x <- snapAdd(x, "将数据集合并 ({bind(sigs)}) 。")
    return(x)
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

.mutate_last_follow_up_time <- function(x) {
  if (any(is.na(x))) {
    freqs <- table(is.na(x), useNA = "ifany")
    message(
      crayon::red(glue::glue("NA: {message_table(freqs)} in 'time' found, as max"))
    )
    max <- max(x, na.rm = TRUE)
    ifelse(is.na(x), max, x)
  } else {
    x
  }
}

multi_test_lasso <- function(x, ..., cut = 15,
  n = 50, seeds = sample(seq_len(50000), n))
{
  jobs <- pbapply::pblapply(seeds, 
    function(seed) {
      x@step <- 3L
      x$seed <- seed
      step4(x, ...)
    })
  ns <- vapply(jobs, 
    function(x) {
      any(x$nfeature_lambdas <= cut) && all(x$nfeature_lambdas > 0)
    }, logical(1))
  list(jobs = jobs[ns], seeds = seeds[ns])
}

message_table <- function(freqs) {
  paste0(
    glue::glue("{names(freqs)} (n = {unname(freqs)})"), 
    collapse = ", "
  )
}
