# ==========================================================================
# workflow of mime
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mime <- setClass("job_mime", 
  contains = c("job"),
  prototype = prototype(
    pg = "mime",
    info = c("https://github.com/l-magnificence/Mime"),
    cite = "[@Mime_A_flexibl_Liu_H_2024]",
    method = "",
    tag = "mime",
    analysis = "Mime 机器学习构建模型"
    ))

setGeneric("asjob_mime", 
  function(x, ...) standardGeneric("asjob_mime"))

setMethod("asjob_mime", signature = c(x = "job_lasso"),
  function(x){
    train_data <- as_tibble(x$train)
    train_data <- dplyr::rename(train_data, ID = rownames)
    train_data <- dplyr::mutate(train_data,
      OS.time = !!x$train.time,
      OS = !!x$train.target,
      .after = 1
    )
    valid_data <- as_tibble(x$valid)
    valid_data <- dplyr::rename(valid_data, ID = rownames)
    valid_data <- dplyr::mutate(valid_data,
      OS.time = !!x$valid.time,
      OS = !!x$valid.target,
      .after = 1
    )
    object <- namel(train_data, valid_data)
    .job_mime(object = object)
  })

setMethod("step0", signature = c(x = "job_mime"),
  function(x){
    step_message("Prepare your data with function `job_mime`.")
  })

setMethod("step1", signature = c(x = "job_mime"),
  function(x, unicox = T, mode = c("all", "single", "double"), seed = 555, ...){
    step_message("Predict results.")
    if (is.null(x$res)) {
      x$mode <- match.arg(mode)
      x$res <- e(Mime1::ML.Dev.Prog.Sig(
          object(x)$train_data,
          object(x),
          candidate_genes = colnames(object(x)$train_data)[-(1:3)],
          unicox.filter.for.candi = unicox,
          unicox_p_cutoff = 0.05,
          mode = x$mode,
          nodesize = 5,
          seed = seed, ...))
    }
    if (unicox) {
      text <- "候选基因以单因素 COX 回归经过初步筛选，设置 p.value cut-off 为 0.05。"
    } else {
      text <- ""
    }
    nt <- nrow(object(x)$train_data)
    nv <- nrow(object(x)$valid_data)
    x <- methodAdd(x, "以 R 包 `Mime` ({packageVersion('Mime1')}) 使用多种机器学习算法 (`Mime1::ML.Dev.Prog.Sig` 默认包含的所有算法)，构建预测模型。训练集和验证集的比例为 {nt/nv}:1。{text}")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_mime"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
