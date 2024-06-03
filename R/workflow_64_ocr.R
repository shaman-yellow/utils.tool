# ==========================================================================
# workflow of ocr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ocr <- setClass("job_ocr", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/mindee/doctr"),
    cite = "",
    method = "Python tool `doctr` <https://github.com/mindee/doctr> used for interpreted character from images",
    tag = "ocr"
    ))

job_ocr <- function()
{
  .job_ocr()
}

setMethod("step0", signature = c(x = "job_ocr"),
  function(x){
    step_message("Prepare your data with function `job_ocr`.")
  })

setMethod("step1", signature = c(x = "job_ocr"),
  function(x, use = "images")
  {
    step_message("Load the Python library")
    use <- match.arg(use)
    # e(reticulate::import_builtins())
    np <- e(reticulate::import("numpy", convert = F))
    doctr.io <- e(reticulate::import("doctr.io"))
    ocr <- e(reticulate::import("doctr.models"))
    cli::cli_alert_info("ocr$ocr_predictor")
    model <- ocr$ocr_predictor(det_arch = "db_resnet50", reco_arch = "crnn_vgg16_bn", pretrained = T)
    ## get the function
    x$model <- model
    x$np <- np
    ## only support unit8
    x$format <- function(x) {
      list(reticulate::np_array(x[[1]], np$uint8))
    }
    if (use == "images") {
      x$parse <- doctr.io$DocumentFile$from_images
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_ocr", ref = "character"),
  function(x, ref)
  {
    if (!all(file.exists(ref))) {
      stop("all(file.exists(`ref`)) == F")
    }
    ref <- normalizePath(ref)
    res <- pbapply::pblapply(ref,
      function(img) {
        obj <- x$parse(img)
        that <- x$model(x$format(obj))
        that$export()
      })
    res
  })

