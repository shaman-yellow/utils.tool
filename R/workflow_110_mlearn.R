# ==========================================================================
# workflow of mlearn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mlearn <- setClass("job_mlearn", 
  contains = c("job"),
  prototype = prototype(
    pg = "mlearn",
    info = c(""),
    cite = "",
    method = "",
    tag = "mlearn",
    analysis = "Machine-Learning 机器学习"
    ))

setGeneric("asjob_mlearn",
  function(x, ...) standardGeneric("asjob_mlearn"))

setMethod("asjob_mlearn", signature = c(x = "job_deseq2"),
  function(x, ref, group = "group", levels = rev(.guess_compare_deseq2(x, 1L)), seed = 987456L)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    object <- x$vst
    if (is.null(object)) {
      stop('is.null(object).')
    }
    data <- SummarizedExperiment::assay(object)
    if (is(ref, "feature")) {
      snap <- snap(ref)
      ref <- resolve_feature(ref)
      if (length(ref) <= 2) {
        stop('length(ref) <= 2, too few genes.')
      }
      snapAdd_onExit("x", "将{snap}用于机器学习筛选关键基因。")
    }
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    data <- t(data[ rownames(data) %in% ref, ])
    metadata <- data.frame(object@colData)
    levels <- eval(levels)
    x <- .job_mlearn(object = data)
    x$metadata <- metadata
    x$levels <- levels
    x$target <- factor(metadata[[ group ]], levels = levels)
    x$seed <- seed
    return(x)
  })

setMethod("asjob_mlearn", signature = c(x = "job_limma"),
  function(x, ref, group = "group", levels = rev(.guess_compare_limma(x, 1L)), seed = 987456L)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    object <- x$normed_data
    if (is.null(object)) {
      stop('is.null(object).')
    }
    data <- object$E
    if (is(ref, "feature")) {
      snap <- snap(ref)
      ref <- resolve_feature(ref)
      if (length(ref) <= 2) {
        stop('length(ref) <= 2, too few genes.')
      }
      snapAdd_onExit("x", "将{snap}用于机器学习筛选关键基因。")
    }
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    data <- t(data[ rownames(data) %in% ref, ])
    metadata <- data.frame(object$targets)
    levels <- eval(levels)
    x <- .job_mlearn(object = data)
    x$metadata <- metadata
    x$levels <- levels
    x$target <- factor(metadata[[ group ]], levels = levels)
    x$seed <- seed
    return(x)
  })

setMethod("step0", signature = c(x = "job_mlearn"),
  function(x){
    step_message("Prepare your data with function `job_mlearn`.")
  })

setMethod("step1", signature = c(x = "job_mlearn"),
  function(x, subset_sizes = 15:50, n = 10, method = "cv", kernel = "linear", seed = x$seed,
    workers = NULL, ..., rerun = FALSE, skip = FALSE)
  {
    step_message("SVM-RFE.")
    if (skip) {
      return(x)
    }
    data <- object(x)
    target <- x$target
    args <- as.list(environment())
    args$rerun <- args$x <- NULL
    dir.create("tmp", FALSE)
    svm_rfe <- expect_local_data(
      "tmp", "svm_rfe", .run_svm_rfe, args, rerun = rerun
    )
    x$svm_rfe <- svm_rfe
    p.svm <- caret:::ggplot.rfe(svm_rfe) +
      lims(x = range(subset_sizes)) +
      theme_classic()
    p.svm <- set_lab_legend(
      wrap(p.svm, 5.5, 4),
      glue::glue("{x@sig} SVM-RFE candidate subset sizes evaluation"),
      glue::glue("SVM-RFE候选子集正确率曲线|||{n}折交叉验证准确率（{n}x CV Accuracy）随特征数量变化的趋势。")
    )
    svm_rfe_res <- list(best_size = svm_rfe$bestSubset, features = caret::predictors(svm_rfe))
    x$svm_rfe_res <- svm_rfe_res
    x$t.svm_rfe_accuracy <- dplyr::arrange(
      as_tibble(svm_rfe$results), dplyr::desc(Accuracy)
    )
    x <- plotsAdd(x, p.svm)
    x <- methodAdd(x, "以 R 包 `e1071` ⟦pkgInfo('e1071')⟧ 构建支持向量机递归特征消除模型 (SVM-RFE)。核函数设定为线性核 (kernel = {kernel}) ；以 R 包 `caret` ⟦pkgInfo('caret')⟧ 实现递归特征消除流程，采用 {n} 折交叉验证评估模型分类性能，依据分类准确率 (Accuracy) 从高到低排序筛选最优特征基因子集。")
    x <- snapAdd(
      x, "SVM-RFE 最佳子集数为 {svm_rfe_res$best_size}{aref(p.svm)}，准确率 (Accuracy) 为 {round(x$t.svm_rfe_accuracy$Accuracy[1], 3)}，误差值 (AccuracySD) 为 {round(x$t.svm_rfe_accuracy$AccuracySD[1], 3)}，对应 feature 为：{bind(svm_rfe_res$features)}。\n\n\n\n"
    )
    return(x)
  })

setMethod("step2", signature = c(x = "job_mlearn"),
  function(x, n = 10, lambda.type = c("1se", "min"), seed = x$seed, ...)
  {
    step_message("Lasso")
    lambda.type <- match.arg(lambda.type)
    lambda.type <- paste0("lambda.", lambda.type)
    data <- object(x)
    target <- x$target
    set.seed(seed)
    cv_lasso <- e(glmnet::cv.glmnet(
        x = data, y = target,
        family = "binomial", alpha = 1, nfolds = n,
        # type.measure = "deviance",
        standardize = TRUE, parallel = FALSE
        ))
    lambda <- cv_lasso[[ lambda.type ]]
    coefs <- coef(cv_lasso, s = lambda)
    coefs_matrix <- as.matrix(coefs)
    whichCoefs <- which(coefs_matrix[, 1] != 0 & rownames(coefs_matrix) != "(Intercept)")
    selected <- rownames(coefs_matrix)[ whichCoefs ]
    # coef_values <- coefs_matrix[selected, 1]
    x$lasso_res <- list(
      cv_lasso = cv_lasso, coefs = coefs_matrix, features = selected, type = lambda.type
    )
    expr <- expression({
      fun <- function() {
        cv <- cv_lasso
        requireNamespace("glmnet")
        suffix <- c("1se", "min")
        types <- paste0("lambda.", suffix)
        y <- max(cv$cvm)
        lambdas <- vapply(types, function(x) cv[[x]], double(1))
        x <- log(lambdas)
        labels <- glue::glue("log(λ) ({suffix})\n = {signif(log(lambdas), 2)}")
        plot(cv, sign.lambda = 1)
        text(x, y, labels, adj = 1)
      }
      fun()
    })
    p.lasso_cv <- as_grob(expr, environment())
    p.lasso_cv <- set_lab_legend(
      wrap(p.lasso_cv, 5.5, 4, showtext = TRUE),
      glue::glue("{x@sig} LASSO Cross Validation"),
      glue::glue("LASSO 交叉验证误差|||Lasso 回归模型的交叉验证图，用于选择正则化参数 λ。图中展示了不同 λ 值下的二项式偏差（Binomial Deviance）。横坐标是log(λ)，即正则化参数 λ 的对数值。随着 λ 值的增加，模型的复杂度降低，正则化强度增加。纵坐标是二项式偏差。")
    )
    expr <- expression({
      fun <- function() {
        cv <- cv_lasso
        requireNamespace("glmnet")
        suffix <- c("1se", "min")
        types <- paste0("lambda.", suffix)
        lambdas <- vapply(types, function(x) cv[[x]], double(1))
        x <- log(lambdas)
        if (any(formalArgs(glmnet:::plot.glmnet) == "sign.lambda")) {
          plot(cv$glmnet.fit, sign.lambda = 1, label = FALSE, xvar = "lambda")
        } else {
          plot(cv$glmnet.fit, label = FALSE, xvar = "lambda")
        }
        abline(v = x, lty = 2)
        labels <- glue::glue("log(λ) ({suffix})\n = {signif(log(lambdas), 2)}")
        text(x, par("usr")[4] * 0.7, labels, adj = 1)
      }
      fun()
    })
    p.coefs_path <- as_grob(expr, environment())
    p.coefs_path <- set_lab_legend(
      wrap(p.coefs_path, 5.5, 4, showtext = TRUE),
      glue::glue("{x@sig} Lasso Coefficient path"),
      glue::glue("LASSO 系数路径|||Lasso 回归系数路径图，展示了不同特征的系数随正则化参数 log(λ) 变化的情况。横坐标是 log(λ)，纵坐标是模型中各个特征的系数值。随着 λ 值的增加（从右到左），更多的特征系数被压缩至零，这是Lasso回归的特征选择过程。")
    )
    x <- plotsAdd(x, p.lasso_cv, p.coefs_path)
    prin <- if (lambda.type == "lambda.1se") "1-SE" else "最小误差"
    x <- methodAdd(x, "以 R 包 glmnet ⟦pkgInfo('glmnet')⟧ 开展 LASSO 逻辑回归分析。设置 α = 1 实现 L1 正则化，通过 {n} 折交叉验证结合 {prin} 准则确定最优 λ 值。")
    x <- snapAdd(x, "LASSO 筛选的核心 feature（非零系数）数量为 {length(selected)}{aref(p.lasso_cv)}，对应为：{bind(selected)}。\n\n\n\n")
    return(x)
  })

setMethod("step3", signature = c(x = "job_mlearn"),
  function(x, ntree = 1000, top = 10, seed = x$seed, ...)
  {
    step_message("Random Forest.")
    data <- object(x)
    target <- x$target
    mtry = floor(sqrt(ncol(data)))
    set.seed(seed)
    rf_model <- e(randomForest::randomForest(
      x = data, y = target, ntree = ntree, mtry = mtry,
      importance = TRUE, proximity = FALSE,
      oob.prox = FALSE, keep.forest = TRUE
    ))
    error_data <- as_tibble(rf_model$err.rate)
    error_data <- dplyr::mutate(error_data, trees = seq_len(nrow(error_data)))
    error_data <- tidyr::pivot_longer(error_data, -trees, names_to = "Error_Type", values_to = "Error_Rate")
    p.error <- ggplot(error_data, aes(x = trees, y = Error_Rate, color = Error_Type)) +
      geom_line() +
      labs(x = "Number of trees", y = "Error Rate") +
      theme_minimal() + theme(legend.title = element_blank())
    p.error <- set_lab_legend(
      wrap(p.error, 5, 3),
      glue::glue("{x@sig} Trend of random forest error rate"),
      glue::glue("随机森林误差率随树数量变化趋势图|||在训练过程中模型对不同组别识别的错误概。OOB 为总体袋外误差（OOB error），即所有类别的平均误差率。随着树的数量增加，总体误差率逐渐趋于稳定。")
    )
    importance_df <- as_tibble(
      randomForest::importance(rf_model), idcol = "feature"
    )
    importance_df <- dplyr::arrange(importance_df, dplyr::desc(MeanDecreaseGini))
    t.tops <- head(importance_df, top)
    t.tops <- set_lab_legend(
      t.tops,
      glue::glue("{x@sig} top importance feature"),
      glue::glue("按 MeanDecreaseGini 降低排序的Top Feature。")
    )
    x$rf_res <- list(rf_model = rf_model, features = t.tops$feature)
    x <- tablesAdd(x, t.tops)
    x <- plotsAdd(x, p.error)
    x <- methodAdd(x, "以 R 包 `randomForest` ⟦pkgInfo('randomForest')⟧ 构建随机森林分类模型，设定决策树数量（ntree）为 {ntree}，特征选择数 (mtry) 为基因总数的平方根，通过袋外数据 (OOB) 评估模型误差率，计算 Feature 重要性评分，筛选相对重要性 top {top}；同时分析分类树数量与误差率的关联趋势，确定模型最优复杂度。")
    x <- snapAdd(x, "随机森林特征重要性 Top {top} 基因：{bind(t.tops$feature)}{aref(p.error)}。\n\n\n\n")
    return(x)
  })

setMethod("step4", signature = c(x = "job_mlearn"),
  function(x, n = 10, seed = x$seed, rerun = FALSE){
    step_message("XGBoost")
    data <- object(x)
    target <- x$target
    fun_xgb <- function(...) {
      .mlearn_alter_xgboost(
        data, target, nfold = n, seed = seed
      )
    }
    res <- expect_local_data(
      "tmp", "xgboost", fun_xgb,
      list(n, seed, target), rerun = rerun
    )
    eval_log <- res$cv$evaluation_log
    p.importance <- ggplot(res$importance, aes(x = reorder(Feature, Gain), y = Gain)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      xlab("Gene") +
      ylab("Importance (Gain)") +
      ggtitle("Feature Importance")
    p.importance <- set_lab_legend(
      wrap_scale(p.importance, 20, nrow(res$importance), size = .1),
      glue::glue("{x@sig} XGBoost Feature Importance"),
      glue::glue("XGBoost 模型特征重要性排序|||展示对分类任务贡献度最高的基因（按 Gain 值衡量），横轴表示特征在模型中的相对重要性，纵轴为基因名称，重要性越高表示该基因在模型决策中贡献越大。")
    )
    data <- dplyr::select(
      eval_log, iter, Validate = test_auc_mean, Train = train_auc_mean
    )
    data <- tidyr::pivot_longer(data, -iter, names_to = "type", values_to = "value")
    p.auc <- ggplot(data, aes(x = iter, y = value, color = type)) +
      geom_line() +
      theme_bw() +
      labs(x = "Iteration", y = "AUC", color = "Type") +
      ggtitle("Cross-validation AUC")
    p.auc <- set_lab_legend(
      wrap(p.auc, 5, 3.5),
      glue::glue("{x@sig} XGBoost Cross-validation AUC"),
      glue::glue("交叉验证模型性能迭代曲线图|||横轴为迭代轮数（Number of Trees），纵轴为模型在验证集上的 AUC 值，用于评估模型随训练过程的收敛趋势及最优迭代轮数的选择。")
    )
    x <- methodAdd(
      x, "以 R 包 `xgboost` ⟦pkgInfo('xgboost')⟧ 构建梯度提升树二分类模型，设定最大迭代轮数（nrounds）为 100，学习率（eta）为 0.05，最大树深（max_depth）为 4，并结合 {n} 折交叉验证与早停策略（early stopping）确定最优迭代轮数；基于模型计算特征重要性（Gain），筛选重要基因；同时分析迭代轮数与模型性能（AUC）变化趋势，以评估模型收敛过程与复杂度。"
    )
    features <- res$selected_genes
    x <- snapAdd(x, "XGBoost 所有重要基因 (n = {length(features)})：{bind(features)}{aref(p.importance)}。\n\n\n\n")
    x$xgb_res <- list(
      rf_model = res$model, rf_cv = res$cv, features = features
    )
    x <- plotsAdd(x, p.auc, p.importance)
    return(x)
  })

.mlearn_alter_xgboost <- function(
  data, target,
  nrounds = 100, nfold = 10,
  early_stopping_rounds = 5,
  seed = 123)
{
  set.seed(seed)
  # -----------------------------
  # Prepare data
  # -----------------------------
  X <- as.matrix(data)
  if (is.factor(target)) {
    y <- as.integer(target) - 1
  } else {
    y <- target
  }

  dtrain <- e(xgboost::xgb.DMatrix(data = X, label = y))

  # -----------------------------
  # Parameters (robust defaults)
  # -----------------------------
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 4,
    eta = 0.05,
    subsample = 0.7,
    colsample_bytree = 0.4,
    lambda = 1,
    alpha = 0
  )

  # -----------------------------
  # Cross-validation
  # -----------------------------
  cv <- e(xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    nfold = nfold,
    early_stopping_rounds = early_stopping_rounds,
    verbose = 1
  ))

  eval_log <- cv$evaluation_log
  best_nrounds <- which.max(eval_log$test_auc_mean)

  # -----------------------------
  # Train final model
  # -----------------------------
  model <- e(xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 1
  ))

  # -----------------------------
  # Feature importance
  # -----------------------------
  importance <- e(xgboost::xgb.importance(
    model = model,
    feature_names = colnames(X)
  ))
  # -----------------------------
  # Output
  # -----------------------------
  list(
    model = model, cv = cv, best_nrounds = best_nrounds,
    importance = importance,
    selected_genes = importance$Feature
  )
}

run_mlean_with_seeds <- function(x, ..., ntry = 100, expect = 2)
{
  x@sig <- "test"
  x <- copy_job(x)
  lapply(seq_len(ntry),
    function(n) {
      seeds <- sample(1:100000, 3)
      capture.output({
        x <- suppressMessages(step1(x, seed = seeds[1], ...))
        x <- suppressMessages(step2(x, seed = seeds[2], ...))
        x <- suppressMessages(step3(x, seed = seeds[3], ...))
      })
      alls <- list(
        SVM_RFE = x$svm_rfe_res$features,
        LASSO = x$lasso_res$features,
        Random_Forest = x$rf_res$features
      )
      res <- ins(lst = alls)
      message(glue::glue("N = {n}, use seeds: {bind(seeds)}, got: {bind(res)}"))
      if (length(res) >= expect) {
        return(seeds)
      } else {
        NULL
      }
    })
}

setMethod("asjob_venn", signature = c(x = "job_mlearn"),
  function(x){
    job_venn(lst = feature(x))
  })

setMethod("feature", signature = c(x = "job_mlearn"),
  function(x){
    lst <- list(SVM_RFE = x$svm_rfe_res$features,
      LASSO = x$lasso_res$features,
      Random_Forest = x$rf_res$features,
      XGBoost = x$xgb_res$features
    )
    lst <- lst[ !vapply(lst, is.null, logical(1)) ]
    as_feature(lst, "Machine Learning")
  })

.run_svm_rfe <- function(data, target, n,
  method = "cv", kernel = "linear", cost = 10,
  subset_sizes, workers = NULL, seed = 123, ...)
{
  # -----------------------------
  # Basic checks & preprocessing
  # -----------------------------
  set.seed(seed)

  if (!is.factor(target)) {
    stop('!is.factor(target).')
  }

  # Filter valid subset sizes
  subset_sizes <- subset_sizes[subset_sizes <= ncol(data)]
  if (length(subset_sizes) == 0) {
    stop("No valid subset_sizes after filtering.")
  }

  # -----------------------------
  # Seeds for reproducibility
  # -----------------------------
  seeds <- vector("list", length = n + 1)
  size <- length(subset_sizes) + 1
  for (i in seq_len(n)) {
    seeds[[i]] <- sample.int(100000L, size)
  }
  seeds[[n + 1]] <- sample.int(100000L, 1)

  # -----------------------------
  # Custom SVM functions for RFE
  # -----------------------------
  svm_funcs <- caret::caretFuncs

  svm_funcs$fit <- function(x, y, first, last, ...) {
    e1071::svm(x = x, y = y,
      kernel = kernel, scale = TRUE, cost = cost,
      probability = FALSE, ...
    )
  }

  svm_funcs$pred <- function(object, x) {
    predict(object, x)
  }

  svm_funcs$rank <- function(object, x, y) {
    # Feature importance only valid for linear kernel
    if (kernel != "linear" || is.null(object$coefs)) {
      return(data.frame(
        var = colnames(x),
        Overall = rep(NA_real_, ncol(x))
      ))
    }

    # Compute weight vector: w = t(coefs) %*% SV
    w <- t(object$coefs) %*% object$SV
    importance <- abs(as.vector(w))

    # Safety check
    if (length(importance) != ncol(x)) {
      importance <- rep(NA_real_, ncol(x))
    }

    data.frame(var = colnames(x), Overall = importance)
  }

  # -----------------------------
  # RFE control
  # -----------------------------
  ctrl <- caret::rfeControl(
    functions = svm_funcs,
    method = method,
    number = n,
    seeds = seeds,
    verbose = TRUE,
    allowParallel = TRUE
  )

  # -----------------------------
  # Parallel setup
  # -----------------------------
  if (!is.null(workers)) {
    workers <- min(workers, parallel::detectCores() - 1)

    cl <- e(parallel::makeCluster(workers))
    lib_paths <- .libPaths()
    parallel::clusterExport(cl, "lib_paths", envir = environment())

    # Ensure workers have required packages
    e(parallel::clusterEvalQ(cl, {
      .libPaths(lib_paths)
      require(caret)
      require(e1071)
    }))

    e(doParallel::registerDoParallel(cl))
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  }

  # -----------------------------
  # Run RFE
  # -----------------------------
  res <- e(caret::rfe(x = data, y = target,
    sizes = subset_sizes, rfeControl = ctrl,
    metric = "Accuracy", maximize = TRUE
  ))

  return(res)
}



# .run_svm_rfe <- function(data, target, n, method, kernel, 
#   subset_sizes, workers, seed)
# {
#   set.seed(seed)
#   seeds <- vector(mode = "list", length = n + 1)
#   size <- length(subset_sizes) + 1
#   for (i in seq_len(n)) {
#     seeds[[i]] <- sample.int(100000L, size)
#   }
#   seeds[[ n + 1 ]] <- sample.int(100000L, 1)
#   ctrl <- e(caret::rfeControl(functions = caret::caretFuncs,
#       method = method, number = n, seeds = seeds, verbose = TRUE, allowParallel = TRUE))
#   svm_funcs <- caret::caretFuncs
#   svm_funcs$fit <- function(x, y, first, last, ...) {
#     e1071::svm(x, y, kernel = kernel, scale = TRUE, probability = TRUE, ...)
#   }
#   svm_funcs$pred <- function(object, x) {
#     predict(object, x)
#   }
#   svm_funcs$rank <- function(object, x, y) {
#     # Calculate feature weights (coefficients of linear kernels)
#     if (is.null(object$coefs)) {
#       # If there are no coefficients, return a random ranking
#       data.frame(var = colnames(x), Overall = runif(ncol(x)))
#     } else {
#       # Calculate feature weights: w = t(x) %*% coefs
#       w <- t(object$coefs) %*% object$SV
#       if (ncol(w) == ncol(x)) {
#         importance <- abs(as.numeric(w))
#       } else {
#         importance <- rep(0, ncol(x))
#       }
#       data.frame(var = colnames(x), Overall = importance)
#     }
#   }
#   subset_sizes <- subset_sizes[subset_sizes <= (ncol(data) - 1)]
#   # run
#   # if (missing(workers) || is.null(workers)) {
#   #   n_cores <- e(parallel::detectCores()) - 1
#   # } else {
#   #   n_cores <- workers
#   # }
#   # if (n_cores > 10) {
#   #   stop('n_cores > 10, too many cores set to run.')
#   # }
#   if (!is.null(workers)) {
#     cl <- e(parallel::makeCluster(workers))
#     e(doParallel::registerDoParallel(cl))
#   }
#   res <- e(caret::rfe(x = data, y = target,
#       sizes = subset_sizes, rfeControl = ctrl,
#       metric = "Accuracy", maximize = TRUE, funcs = svm_funcs
#       ))
#   if (!is.null(workers)) {
#     e(parallel::stopCluster(cl))
#     e(foreach::registerDoSEQ())
#   }
#   res
# }

setMethod("set_remote", signature = c(x = "job_mlearn"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
