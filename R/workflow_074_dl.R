# ==========================================================================
# workflow of dl
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_dl <- setClass("job_dl", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "dl",
    info = c("https://github.com/JinYSun/D-GCAN"),
    cite = "[@PredictionOfDSunJ2022]",
    method = "The deep learning method `D-GCAN` (python) was used for prediction of drug-likeness",
    tag = "dl",
    analysis = "D-GCAN Drug-Likeness 预测"
    ))

job_dl <- function(smiles, force_cpu = TRUE)
{
  x <- .job_dl()
  if (is.data.frame(smiles)) {
    message("The input is data.frame. Use first column as names, next column as smiles.")
    smiles <- nl(smiles[[1]], smiles[[2]], FALSE)
  }
  smiles <- smiles[ !duplicated(smiles) ]
  ## check smiles character
  if (any(isDot <- grpl(smiles, "\\."))) {
    message("SMILES with character '.' will be excluded.\n", "Excluded: ", length(smiles[isDot]))
    x$excluded <- smiles[ isDot ]
    smiles <- smiles[ !isDot ]
  }
  object(x) <- smiles
  activate_base("dl")
  sys <- reticulate::import("sys")
  sys$path <- c(sys$path, pg("dl"))
  if (force_cpu) {
    message("Use CPU for computing.")
    x$force_cpu <- TRUE
    fun <- function(name) {
      file_py <- file.path(pg("dl"), paste0(name, ".py"))
      raw <- readLines(file_py)
      revise <- gs(raw, "torch.device\\('cuda'\\)", "torch.device('cpu')")
      revise <- gs(revise, "print\\('The code uses a GPU!'\\)", "print('The code uses a CPU!')")
      message("Revised: \n\t", file_py)
      writeLines(revise, file_py)
      namel(raw, revise, file_py)
    }
    x$scripts <- sapply(c("train", "predict", "preprocess", "DGCAN"), fun, simplify = FALSE)
  } else {
    x$force_cpu <- FALSE
  }
  predict <- reticulate::import("predict")
  train <- reticulate::import("train")
  if (force_cpu) {
    fun <- function(x) {
      writeLines(x$raw, x$file_py)
      message("The codes has been recovered:\n\t", x$file_py)
    }
    lapply(x$scripts, fun)
  }
  x$sys <- sys
  x$predict <- predict
  x$train <- train
  return(x)
}

setMethod("step0", signature = c(x = "job_dl"),
  function(x){
    step_message("Prepare your data with function `job_dl`.")
  })

setMethod("step1", signature = c(x = "job_dl"),
  function(x, batch_size = 3L) {
    step_message("Prepare the trained model.")
    file_model <- file.path(pg("dl"), "model", "model.pth")
    if (!file.exists(file_model)) {
      message("Model file not exists, training herein.")
      dir.create(file.path(pg("dl"), "model"))
      owd <- getwd()
      setwd(pg("dl"))
      x$res_train <- tryCatch(x@params$train$train(
          file.path(pg("dl_dataset"), "bRo5.txt"),
          radius = 1L,
          dim = 52L,
          layer_hidden = 4L,
          layer_output = 10L,
          dropout = 0.45,
          batch_train = batch_size,
          batch_test = batch_size,
          lr = 3e-4,
          lr_decay = 0.85,
          decay_interval = 25L,
          iteration = 140L,
          N = 5000L,
          dataset_train = file.path(pg("dl_dataset"), "data_train.txt")
          ), finally = setwd(owd))
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_dl"),
  function(x){
    step_message("Use Model to predict.")
    file_smiles <- tempfile("smiles", fileext = ".txt")
    writeLines(object(x), file_smiles)
    owd <- getwd()
    setwd(pg("dl"))
    res <- tryCatch(x@params$predict$predict(
      normalizePath(file_smiles),
      # the following parameter is useless and can be ignore.
      radius = 1L, property = FALSE, dim = 52L,
      layer_hidden = 4L, layer_output = 10L, dropout = 0.45,
      batch_train = 8L, batch_test = 8L, lr = 3e-4,
      lr_decay = 0.85, decay_interval = 25L, iteration = 140L, N = 5000L
    ), finally = setwd(owd))
    if (length(res) != length(object(x))) {
      stop("length(res) != length(object(x))")
    }
    message("Get length of results: ", length(res))
    t.res <- tibble::tibble(name = names(object(x)), smiles = object(x), drugLike = ifelse(res, TRUE, FALSE))
    t.res <- dplyr::mutate(t.res, isOK = drugLike)
    x <- methodAdd(x, "以 Python 工具 `D-GCAN` {cite_show('PredictionOfDSunJ2022')} 预测 Drug-likeness。
      `D-GCAN` 的训练和使用参数参考 <https://github.com/JinYSun/D-GCAN>。")
    x <- snapAdd(x, "以化合物结构式 (SMILES) 通过 `D-GCAN` 程序预测是否与药物相似 (Drug-likeness)。")
    x <- snapAdd(x, "`D-GCAN` 预测结果，所有用于预测的化合物 (含有结构式信息的) {nrow(t.res)} 个，与药物相似的 {length(which(t.res$isOK))} 个 (注：根据唯一结构式统计)。")
    t.res <- .set_lab(t.res, sig(x), "prediction of Drug-Likeness data")
    p.res <- new_pie(c(as.character(t.res$drugLike), rep("Excluded", length(x$excluded))))
    p.res <- .set_lab(p.res, sig(x), "prediction of Drug-Likeness")
    x@tables[[ 2 ]] <- namel(t.res)
    x@plots[[ 2 ]] <- namel(p.res)
    return(x)
  })

setMethod("res", signature = c(x = "job_dl"),
  function(x){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    x@tables$step2$t.res
  })


