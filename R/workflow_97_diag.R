# ==========================================================================
# workflow of lst_diag
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_diag <- setClass("job_diag", 
  contains = c("job_lasso"),
  prototype = prototype(
    pg = "lst_diag",
    info = c("..."),
    cite = "",
    method = "",
    tag = "filter:diag",
    analysis = "Lasso 诊断模型建立"
    ))

job_diag <- function()
{
  .job_diag()
}

setGeneric("asjob_diag", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_diag"))

setMethod("asjob_diag", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = .guess_symbol(x), from_normed = TRUE,
    fun_scale = function(x) scale(x, TRUE, TRUE),
    dup_method = c("max", "min", "mean"), 
    use.format = TRUE, ...)
  {
    args <- as.list(environment())
    x <- do.call(asjob_lasso, args)
    ref <- .job_diag()
    x <- .job_diag(
      x, tag = ref@tag, analysis = ref@analysis, cite = ref@cite, 
      info = ref@info, method = ref@method
    )
    return(x)
  })

setMethod("step0", signature = c(x = "job_diag"),
  function(x){
    step_message("Prepare your data with function `job_diag`.")
  })

setMethod("step1", signature = c(x = "job_diag"),
  function(x, target = "group", levels = "guess", pattern_control = "control|ctrl|normal",
    n.train = .8, seed = 555)
  {
    if (identical(levels, "guess")) {
      allLevels <- unique(x$metadata$group)
      if (length(allLevels) == 2) {
        if (any(isThat <- grpl(allLevels, pattern_control, TRUE))) {
          if (length(which(isThat)) == 1) {
            levels <- c(allLevels[isThat], allLevels[!isThat])
          } else {
            stop(glue::glue("Don't know ... {bind(allLevels)}"))
          }
        } else {
          stop("No any 'control' or 'normal' matched.")
        }
      } else {
        stop('length(allLevels) != 2')
      }
    }
    x <- callNextMethod(x, target = target, levels = levels,
      time = NULL, n.train = .8, seed = 555)
    x@snap$step1 <- ""
    return(x)
  })

setMethod("step2", signature = c(x = "job_diag"),
  function(x, use_data = c("all", "train"), top = 30, efs = FALSE){
    x <- callNextMethod()
    return(x)
  })

setMethod("step3", signature = c(x = "job_diag"),
  function(x){
    step_message("Do nothing.")
    return(x)
  })

setMethod("step4", signature = c(x = "job_diag"),
  function(x, use_data = x$use_data, use_valid = use_data,
    fun = c("cv.glmnet", "glmnet"), nfold = 10,
    alpha = 1, family = "binomial", type.measure = c("default"), ...)
  {
    step_message("Lasso ...")
    fun_model <- match.arg(fun)
    type.measure <- match.arg(type.measure)
    data_lst <- x$get(use_data)
    data <- data_lst$data
    valid_lst <- x$get(use_data)
    valid <- valid_lst$data
    x$nfold <- nfold
    target <- data_lst$event
    lst_diag <- list()
    if (fun_model == "cv.glmnet") {
      methodName <- if (alpha == 1) {
        "lasso"
      } else if (!alpha) {
        "ridge"
      } else if (alpha > 0 && alpha < 1) {
        "Elastic Net"
      }
      set.seed(x$seed)
      lst_diag$model <- model <- e(glmnet::cv.glmnet(data, target,
          alpha = alpha, family = family, nfold = nfold, type.measure = type.measure, ...))
      p.lasso_model <- wrap(as_grob(expression(plot(model)), environment()), 6, 6, showtext = TRUE)
      lambdas <- c("lambda.min", "lambda.1se")
      res <- sapply(lambdas, simplify = FALSE,
        function(lam) {
          preds <- e(stats::predict(model, newx = valid, s = model[[ lam ]]))
          roc <- e(pROC::roc(valid_lst$types, as.numeric(preds), plot = FALSE, levels = x$levels))
          namel(preds, roc)
        })
      lst_diag$preds <- lapply(res, function(lst) lst$preds)
      lst_diag$roc <- lapply(res, function(lst) lst$roc)
      if (TRUE) {
        message("Got coefficients.")
        x$lambda <- c(min = model$lambda.min, `1se` = model$lambda.1se)
        x <- snapAdd(x, "使用 `glmnet::cv.glmnet` 作 {nfold} 倍交叉验证 (评估方式为 {model$name})，筛选 lambda 值。lambda.min, lambda.1se 值分别为 {bind(round(x$lambda, 3))} (R 随机种子为 {x$seed})。")
        lst_diag$coef <- coefficients(
          model, s = c(lambda.min = model$lambda.min, lambda.1se = model$lambda.1se)
        )
        x$sig.diag <- sig.diag <- dplyr::rename(
          as_tibble(as.matrix(lst_diag$coef[-1, ])),
          feature = 1, coef.lambda.min = 2, coef.lambda.1se = 3
        )
        if (all(!lst_diag$coef[-1, 2])) {
          message(crayon::red("lambda.1se has non coeffients, use Lambda.min (-> coef)."))
        }
        s.com <- vapply(1:2,
          function(x) {
            n <- lst_diag$coef[-1, x]
            length(n[ n != 0 ])
          }, integer(1))
        x$nfeature_lambdas <- s.com
        x <- snapAdd(x, "对应的特征数 (基因数) 分别为 {bind(s.com)}。")
        x <- tablesAdd(x, t.sigCoefficients = sig.diag)
      }
      if (TRUE) {
        p.lasso_ROC <- lapply(lst_diag$roc,
          function(roc) {
            plot_roc(roc)
          })
        names(p.lasso_ROC) <- lambdas
      }
      if (TRUE) {
        coef <- as_tibble(Matrix::as.matrix(lst_diag$coef))
        colnames(coef)[-1] <- lambdas
        coef <- dplyr::rename(coef, variable = rownames)
        message("Plot Significant coefficients.")
        p.lasso_coeffients <- sapply(
          lambdas, simplify = FALSE,
          function(lam) {
            coef <- dplyr::arrange(coef, dplyr::desc(abs(!!rlang::sym(lam))))
            plot_sig(dplyr::filter(coef, !grepl("Inter", variable)), y = lam)
          }
        )
      }
      x <- plotsAdd(x, p.lasso_model, p.lasso_ROC, p.lasso_coeffients)
      x <- methodAdd(x, "以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归，以 `{fun_model}` 函数作 {nfold} 交叉验证获得模型。", TRUE)
    } else if (fun_model == "glmnet") {
      methodName <- if (alpha) "lasso" else "ridge"
      lst_diag$model <- e(glmnet::glmnet(data, target,
          alpha = alpha, family = family, ...))  
      lst_diag$coef <- coefficients(lst_diag$model)
      x <- methodAdd(x, "以 R 包 `glmnet` ({packageVersion('glmnet')}) 作 {methodName} 处罚的 {family} 回归。", TRUE)
    }
    x$lst_diag <- lst_diag
    x$fun_model <- fun_model
    return(x)
  })

setMethod("map", signature = c(x = "job_diag", ref = "job_diag"),
  function(x, ref, lambda = c("min", "1se"))
  {
    message("Validate the model using external dataset.")
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    if (ref@step < 4L) {
      stop('ref@step < 4L.')
    }
    model <- ref@params$lst_diag$model
    used_vars <- rownames(coefficients(model))[-1]
    valid_lst <- x$get("all")
    if (!identical(used_vars, colnames(valid_lst$data))) {
      message(glue::glue('!identical(used_vars, colnames(valid_lst$data)), try match...'))
      matched <- match(used_vars, colnames(valid_lst$data))
      if (any(is.na(matched))) {
        notGot <- used_vars[ !used_vars %in% colnames(valid_lst$data) ]
        sigCoeff <- dplyr::filter(
          ref@params$sig.diag, dplyr::if_any(dplyr::where(is.double), ~ . != 0)
        )
        if (all(!notGot %in% sigCoeff$feature)) {
          fun_color <- crayon::yellow
          supp <- matrix(
            0, nrow = dim(valid_lst$data)[1], dimnames = list(rownames(valid_lst$data), notGot)
          )
          valid_lst$data <- cbind(valid_lst$data, supp)
        } else {
          fun_color <- crayon::red
        }
        message(fun_color(glue::glue("Some var can not match ({bind(notGot)})")))
      } else {
        valid_lst$data <- valid_lst$data[, matched]
      }
    }
    valid <- valid_lst$data
    lambdas <- paste0("lambda.", lambda)
    x$valid_results <- sapply(lambdas, simplify = FALSE,
      function(lam) {
        preds <- e(stats::predict(model, newx = valid, s = model[[ lam ]]))
        roc <- e(pROC::roc(valid_lst$types, as.numeric(preds), plot = FALSE, levels = x$levels))
        namel(preds, roc)
      })
    x$.map_heading <- "使用外部数据验证"
    return(x)
  })

