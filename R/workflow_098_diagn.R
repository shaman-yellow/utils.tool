# ==========================================================================
# workflow of diagn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_diagn <- setClass("job_diagn", 
  contains = c("job"),
  prototype = prototype(
    pg = "diagn",
    info = c(""),
    cite = "",
    method = "",
    tag = "",
    analysis = "RNA 数据集验证"
    ))

setGeneric("asjob_diagn", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_diagn"))

setMethod("asjob_diagn", signature = c(x = "job_vennDEGs"),
  function(x, genes = feature(x), names = x$degs_versus@object_names,
    pattern_control = "guess", env = .GlobalEnv)
  {
    objects <- lapply(names, get, envir = env)
    y <- job_diagn(objects, genes = genes)
    y@analysis <- "Lasso 诊断模型建立"
    if (identical(pattern_control, "guess")) {
      versus <- strsplit(unlist(x$degs_versus), "\\s?-\\s?")
      versus <- vapply(versus, function(x) x[[ 2 ]], character(1))
      names(versus) <- names(object(y))
      y$pattern_control <- versus
    }
    return(y)
  })

job_diagn <- function(x, genes = NULL, ...)
{
  valid_job_list(x, "job_limma", 1L)
  projects <- vapply(x, function(x) x$project, character(1))
  jobs <- lapply(x,
    function(x) {
      asjob_diag(x, genes, ...)
    })
  names(jobs) <- projects
  x <- .job_diagn(object = jobs)
  x$projects <- projects
  x$genes <- genes
  return(x)
}

setMethod("step0", signature = c(x = "job_diagn"),
  function(x){
    step_message("Prepare your data with function `job_diagn`.")
  })

setMethod("step1", signature = c(x = "job_diagn"),
  function(x, pattern_control = x$pattern_control %||% "control|ctrl|normal",
    target = "group", levels = "guess", ...)
  {
    step_message("Prepare data.")
    if (length(pattern_control) > 1) {
      if (is.null(names(pattern_control))) {
        stop('is.null(names(pattern_control)).')
      }
      if (!all(names(object(x)) %in% names(pattern_control))) {
        stop('!all(names(object(x)) %in% names(pattern_control)).')
      }
      object(x) <- sapply(
        names(object(x)), simplify = FALSE,
        function(name) {
          message(glue::glue("Dataset {name} use pattern_control: {pattern_control[[ name ]]}"))
          step1(
            object(x)[[ name ]], 
            pattern_control = pattern_control[[ name ]], 
            target = target, levels = levels, ...
          )
        }
      )
    } else {
      object(x) <- lapply(object(x), 
        function(x) {
          step1(
            x, pattern_control = pattern_control, target = target, 
            levels = levels, ...
          )
        })
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_diagn"),
  function(x, use_data = c("all", "train"), top = 30, 
    efs = FALSE, ...)
  {
    step_message("Evaluate variable (genes) importance.")
    object(x) <- lapply(object(x), 
      function(x) {
        step2(x, use_data = use_data, top = top, efs = efs, ...)
      })
    return(x)
  })

setMethod("step3", signature = c(x = "job_diagn"),
  function(x, ...){
    step_message("Do nothing.")
    object(x) <- lapply(object(x), 
      function(x) {
        step3(x, ...)
      })
    return(x)
  })

setMethod("step4", signature = c(x = "job_diagn"),
  function(x, fun = c("cv.glmnet", "glmnet"), nfold = 10,
    alpha = 1, family = "binomial", type.measure = c("default"), ...)
  {
    step_message("Lasso ...")
    object(x) <- lapply(object(x), 
      function(x) {
        step4(
          x, fun = fun, nfold = nfold, alpha = alpha, 
          family = family, type.measure = type.measure, ...
        )
      })
    return(x)
  })

setMethod("step5", signature = c(x = "job_diagn"),
  function(x){
    step_message("Integration")
    t.res <- lapply(object(x), 
      function(x) {
        sets <- x$get(x$use_data)
        events <- sets$event
        rocs <- x$lst_diag$roc
        coefs <- x@tables$step4$t.sigCoefficients
        tibble::tibble(
          n = length(events), pos = length(which(events == 1)),
          lambda.min.auc = as.double(rocs$lambda.min$auc),
          lambda.1se.auc = as.double(rocs$lambda.1se$auc),
          lambda.min.coef = length(which(coefs$coef.lambda.min != 0)),
          lambda.1se.coef = length(which(coefs$coef.lambda.1se != 0))
        )
      })
    t.res <- rbind_list(t.res, .id = "project")
    t.res <- setLegend(t.res, "为不同数据集建立的诊断模型的参数。")
    x <- tablesAdd(x, t.allDatasetModels = t.res)
    return(x)
  })

setMethod("step6", signature = c(x = "job_diagn"),
  function(x, mode = c("maxSamples")){
    step_message("")
    mode <- match.arg(mode)
    t.res <- x@tables$step5$t.allDatasetModels
    map_res <- list()
    if (mode == "maxSamples") {
      projModel <- t.res$project[ which.max(t.res$n) ]
      y <- set_independence(x)
      ref <- object(x)[[ projModel ]]
      object(y) <- object(y)[ names(object(y)) != projModel ]
      y <- map(y, ref)
      map_res <- y@params[ names(y@params) %in% c("p.hps", "p.rocs", "rocs") ]
      x$projModel <- projModel
      x <- methodAdd(x, "{ref@meth$step4}")
      x <- snapAdd(x, "以 {projModel} 创建诊断模型。{ref@snap$step4}{y@snap$stepm}")
      if (TRUE) {
        z <- set_independence(ref)
        coefs <- ref@tables$step4$t.sigCoefficients
        features <- sapply(
          c("coef.lambda.min", "coef.lambda.1se"), simplify = FALSE,
          function(cname) {
            coefs$feature[ coefs[[ cname ]] != 0 ]
          }
        )
        feature(x) <- features
        z <- map(z, ref)
        x$self_res <- c(list(p.hp = z@params$p.hp_valid), ref@plots$step4)
      }
    }
    x$map_res <- map_res
    return(x)
  })

setMethod("map", signature = c(x = "job_diagn", ref = "job_diag"),
  function(x, ref, ...){
    object(x) <- lapply(object(x), 
      function(x) {
        map(x, ref, ...)
      })
    x$valid_results <- lapply(object(x), function(x) x$valid_results)
    lambdas <- c("lambda.min", "lambda.1se")
    rocs <- sapply(lambdas, simplify = FALSE, 
      function(lam) {
        roc <- lapply(x$valid_results, 
          function(object) {
            object[[ lam ]]$roc
          })
      })
    p.rocs <- lapply(rocs, plot_roc)
    p.rocs <- .set_lab(p.rocs, sig(x), names(p.rocs), "ROC")
    p.rocs <- setLegend(p.rocs, glue::glue("为 {names(p.rocs)} ROC 曲线。"))
    x <- snapAdd(
      x, "使用数据集 {bind(names(object(x)))} 验证。", add = FALSE, step = "m"
    )
    x$p.hps <- lapply(object(x), function(x) x$p.hp_valid)
    x$p.hps <- .set_lab(x$p.hps, sig(x), names(x$p.hps), "feature heatmap in validation dataset")
    x$p.hps <- setLegend(x$p.hps, glue::glue("为特征基因在验证数据集 ({names(x$p.hps)}) 的表达热图。"))
    x$p.rocs <- p.rocs
    x$rocs <- rocs
    return(x)
  })



