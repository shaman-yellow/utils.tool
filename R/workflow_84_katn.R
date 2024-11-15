# ==========================================================================
# workflow of katn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_katn <- setClass("job_katn", 
  contains = c("job_kat"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype())

setGeneric("asjob_katn", 
  function(x, ...) standardGeneric("asjob_katn"))

setMethod("asjob_katn", signature = c(x = "job_seurat"),
  function(x, use = names(x@object@assays)[[1]], layer = "counts"){
    metadata <- object(x)@meta.data
    cellsGroup <- split(rownames(metadata), metadata$orig.ident)
    data <- do.call(`$`, list(object(x)@assays[[ use ]], layer))
    data <- pbapply::pblapply(cellsGroup,
      function(cells) {
        object <- data[, cells]
        rownames(object) <- gs(rownames(object), "\\.[0-9]*$", "")
        .job_kat(object = object)
      })
    x <- .job_katn(object = data)
    x$metadata <- metadata
    return(x)
  })

setMethod("step0", signature = c(x = "job_katn"),
  function(x){
    step_message("Prepare your data with function `job_katn`.")
  })

setMethod("step1", signature = c(x = "job_katn"),
  function(x, workers = 5, paths = "copykat_batch")
  {
    step_message("Quality control (QC).")
    n <- 0L
    paths <- paste0(paths, "_", names(object(x)))
    lapply(paths,
      function(path) {
        if (dir.exists(path)) {
          stop(glue::glue("{path} exists."))
        } else {
          dir.create(path)
        }
      })
    meth(x)$step1 <- glue::glue("R 包 `CopyKAT` 用于鉴定恶性细胞 {cite_show('DelineatingCopGaoR2021')}。`CopyKAT` 可以区分整倍体与非整倍体，其中非整倍体被认为是肿瘤细胞，而整倍体是正常细胞 {cite_show('CausesAndConsGordon2012')}。由于 `CopyKAT` 不适用于多样本数据 (无法应对批次效应) ，因此，对各个样本独立执行 `CopyKAT`。")
    object(x) <- pbapply::pblapply(object(x),
      function(obj) {
        n <<- n + 1L
        cli::cli_alert_info(glue::glue("Running {n}: {getwd()}"))
        obj <- step1(obj, workers, paths[n])
        obj <- step2(obj)
        obj <- step3(obj)
        gc()
        return(obj)
      })
    return(x)
  })
