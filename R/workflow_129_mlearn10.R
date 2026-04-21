# ==========================================================================
# workflow of mlearn10
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mlearn10 <- setClass("job_mlearn10", 
  contains = c("job"),
  prototype = prototype(
    pg = "mlearn10",
    info = c(""),
    cite = "",
    method = "",
    tag = "mlearn10",
    analysis = ""
    ))

setMethod("step0", signature = c(x = "job_mlearn10"),
  function(x){
    step_message("Prepare your data with function `job_mlearn10`.")
  })

setGeneric("asjob_mlearn10",
  function(x, ...) standardGeneric("asjob_mlearn10"))

setMethod("asjob_mlearn10", signature = c(x = "job_mlearn"),
  function(x, ...)
  {
    .job_mlearn10(x)
  })

setMethod("step1", signature = c(x = "job_mlearn10"),
  function(x, n = 10L, nthread = 1L)
  {
    step_message("Run caret")
    x$res_ml_train10 <- ml_train10(
      object(x), x$target, x$levels[2], cv_folds = n,
      nthread = nthread
    )
    return(x)
  })


# =============================================================================
# Main Function
# =============================================================================
ml_train10 <- function(
  data, target,
  positive_class = NULL,
  cv_folds = 10L,
  tune_length = 5L,
  nthread = 1L,
  seed = 123L
)
{
  # ---------------------------------------------------------------------------
  # package check
  # ---------------------------------------------------------------------------
  pkgs <- c(
    "caret", "glmnet", "randomForest", "kernlab",
    "xgboost", "gbm", "nnet", "rpart"
  )

  miss_pkg <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
    ]

  if (length(miss_pkg) > 0) {
    stop("Missing package(s): ", paste(miss_pkg, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # data check
  # ---------------------------------------------------------------------------
  data <- as.data.frame(data)
  target <- as.factor(target)

  if (nrow(data) != length(target)) {
    stop("nrow(data) != length(target)")
  }

  if (length(levels(target)) != 2) {
    stop("target must be binary")
  }

  lv <- levels(target)

  if (is.null(positive_class)) {
    positive_class <- lv[2]
  }

  if (!positive_class %in% lv) {
    stop("positive_class not found in target")
  }

  target <- factor(
    target,
    levels = c(
      positive_class,
      setdiff(lv, positive_class)
    )
  )

  dat <- data.frame(Class = target, data)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  set.seed(seed)

  # ---------------------------------------------------------------------------
  # caret control
  # ---------------------------------------------------------------------------
  ctrl <- caret::trainControl(
    method = "cv",
    number = cv_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = TRUE
  )

  # ---------------------------------------------------------------------------
  # model config
  # ---------------------------------------------------------------------------
  cfg <- list(

    SVM = list(
      method = "svmLinear",
      tuneLength = tune_length
      ),

    RF = list(
      method = "rf",
      ntree = 1000L,
      tuneLength = 3L
      ),

    EN = list(
      method = "glmnet",
      tuneLength = tune_length
      ),

    XGBoost = list(
      method = "xgbTree",
      tuneLength = tune_length,
      verbose = FALSE,
      nthread = nthread
      ),

    Lasso = list(
      method = "glmnet",
      tuneGrid = expand.grid(
        alpha = 1,
        lambda = seq(0.0001, 1, length.out = 20)
      )
      ),

    ANN = list(
      method = "nnet",
      trace = FALSE,
      MaxNWts = 1000000,
      tuneLength = tune_length
      ),

    DT = list(
      method = "rpart",
      tuneLength = tune_length
      ),

    RR = list(
      method = "glmnet",
      tuneGrid = expand.grid(
        alpha = 0,
        lambda = seq(0.0001, 1, length.out = 20)
      )
      ),

    GBM = list(
      method = "gbm",
      verbose = FALSE,
      tuneLength = tune_length
    )
  )

  # ---------------------------------------------------------------------------
  # train caret models
  # ---------------------------------------------------------------------------
  model_list <- list()

  # for cache
  identity_args <- list(
    rownames(data), colnames(data), target,
    cv_folds, lv, positive_class
  )

  base_args <- list(
    form = Class ~ .,
    data = dat,
    metric = "ROC",
    trControl = ctrl
  )

  fun_train <- function(...) {
    tryCatch(
      do.call(caret::train, c(base_args, args)),
      error = function(e) {
        message(glue::glue("Error in training: {nm}"))
        NULL
      }
    )
  }

  for (nm in names(cfg)) {
    args <- cfg[[nm]]
    message(glue::glue("Training model: {nm}"))

    if (nm != "XGBoost") {
      fit <- expect_local_data(
        "tmp", glue::glue("mlearn10_{nm}"),
        fun_train, identity_args
      )
    } else {
      fit <- NULL
    }

    model_list[[nm]] <- fit
  }

  # ---------------------------------------------------------------------------
  # LightGBM
  # ---------------------------------------------------------------------------
  if (requireNamespace("lightgbm", quietly = TRUE)) {

    x_mat <- as.matrix(dat[, -1, drop = FALSE])

    y_vec <- ifelse(
      dat$Class == levels(dat$Class)[1],
      1, 0
    )

    dtrain <- lightgbm::lgb.Dataset(
      data = x_mat,
      label = y_vec
    )

    params <- list(
      objective = "binary",
      metric = "auc",
      learning_rate = 0.05,
      num_leaves = 31,
      feature_fraction = 0.9,
      bagging_fraction = 0.8,
      bagging_freq = 5,
      num_threads = nthread,
      verbose = -1
    )

    cv_res <- tryCatch(
      lightgbm::lgb.cv(
        params = params,
        data = dtrain,
        nrounds = 100,
        nfold = cv_folds,
        stratified = TRUE,
        verbose = -1
        ),
      error = function(e) NULL
    )

    if (!is.null(cv_res)) {

      auc_vec <- unlist(
        cv_res$record_evals$valid$auc$eval
      )

      best_auc <- max(auc_vec, na.rm = TRUE)

      best_iter <- cv_res$best_iter

      if (is.null(best_iter)) {
        best_iter <- which.max(auc_vec)
      }

      lgb_model <- lightgbm::lgb.train(
        params = params,
        data = dtrain,
        nrounds = best_iter,
        verbose = -1
      )

      lgb_wrapper <- list(
        model = lgb_model,
        bestTune = data.frame(
          nrounds = best_iter
          ),
        modelInfo = list(
          label = "LightGBM"
          ),
        levels = levels(dat$Class),
        feature_names = colnames(x_mat),
        cv_auc = best_auc
      )

      class(lgb_wrapper) <- c(
        "lgb.Booster.train",
        "caret.train.like"
      )

      model_list$LightGBM <- lgb_wrapper

    } else {

      model_list$LightGBM <- NULL
    }

  } else {

    model_list$LightGBM <- NULL
  }

  # ---------------------------------------------------------------------------
  # summary
  # ---------------------------------------------------------------------------
  res <- data.frame(
    Model = names(model_list),
    CV_AUC = NA_real_,
    stringsAsFactors = FALSE
  )

  features <- list()

  for (i in seq_along(model_list)) {

    fit <- model_list[[i]]
    if (is.null(fit)) next
    nm <- names(model_list)[i]
    if (any(nm == c("Lasso", "EN", "RR"))){
      co <- coef(fit$finalModel, s = fit$bestTune$lambda)
      co <- as.matrix(co)
      gene <- rownames(co)[co[, 1] != 0]
      gene <- setdiff(gene, "(Intercept)")
      features[[nm]] <- gene
    } else if(any(nm == c("RF", "GBM", "DT"))){
      vi <- caret::varImp(fit)$importance
      score <- vi[, 1]
      names(score) <- rownames(vi)
      score <- sort(score, decreasing = TRUE)
      features[[nm]] <- names(score)[1:min(10, length(score))]
    } else if (any(nm == "LightGBM")) {
      imp <- lightgbm::lgb.importance(fit$model)
      features$LightGBM <- imp$Feature
    }
    if (inherits(fit, "train")) {
      res$CV_AUC[i] <- max(
        fit$results$ROC,
        na.rm = TRUE
      )
    } else if (inherits(fit, "lgb.Booster.train")) {
      res$CV_AUC[i] <- fit$cv_auc
    }
  }
  res <- res[order(-res$CV_AUC), ]
  # ---------------------------------------------------------------------------
  # return
  # ---------------------------------------------------------------------------
  out <- list(
    features = features,
    models = model_list,
    summary = res,
    data = dat,
    positive_class = positive_class,
    levels = levels(target)
  )
  class(out) <- "ml_train10"
  return(out)
}

# =============================================================================
# Custom S3 predict method for LightGBM wrapper
# =============================================================================
predict.lgb.Booster.train <- function(object, newdata, type = "prob", ...) {

  if (!inherits(newdata, "matrix")) {
    newdata <- as.matrix(newdata)
  }

  miss_col <- setdiff(object$feature_names, colnames(newdata))
  if (length(miss_col) > 0) {
    stop("Missing columns in newdata: ",
      paste(miss_col, collapse = ", "))
  }

  newdata <- newdata[, object$feature_names, drop = FALSE]

  raw_pred <- predict(object$model, data = newdata)

  if (type == "prob") {

    out <- data.frame(
      V1 = 1 - raw_pred,
      V2 = raw_pred
    )

    colnames(out) <- object$levels
    return(out)

  } else if (type == "raw") {

    return(raw_pred)

  } else {

    stop("type must be 'prob' or 'raw'")
  }
}

