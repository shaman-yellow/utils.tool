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

setGeneric("asjob_katn", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_katn"))

setMethod("asjob_katn", signature = c(x = "job_seurat"),
  function(x, refs = NULL, group.by = x$group.by,
    use = names(x@object@assays)[[1]], split.by = "orig.ident", layer = "counts")
  {
    metadata <- object(x)@meta.data
    cellsGroup <- split(rownames(metadata), metadata[[ split.by ]])
    data <- do.call(`$`, list(object(x)@assays[[ use ]], layer))
    data <- pbapply::pbmapply(
      cellsGroup, names(cellsGroup), FUN = function(cells, name) {
        object <- data[, cells]
        rownames(object) <- gname(rownames(object))
        refs <- dplyr::filter(
          metadata, !!rlang::sym(group.by) %in% refs, !!rlang::sym(split.by) == name
        )
        if (!nrow(refs)) {
          stop('!nrow(refs).')
        }
        refs <- rownames(refs)
        job_kat(object, refs)
      }
    )
    x <- .job_katn(object = data)
    x$metadata <- metadata
    return(x)
  })

setValidity("job_katn",
  function(object){
    .valid_pblapply_res(object)
  })

.valid_pblapply_res <- function(object) {
  if (is(object(object), "list")) {
    valid <- !any(isError <- vapply(object(object),
        function(obj) {
          inherits(obj, "try-error")
        }, logical(1)))
    if (!valid) {
      message(crayon::red(glue::glue("valid: {try_snap(!isError)}")))
      print(object(object)[isError])
    }
    valid
  } else {
    TRUE
  }
}

setMethod("step0", signature = c(x = "job_katn"),
  function(x){
    step_message("Prepare your data with function `job_katn`.")
  })

setMethod("step1", signature = c(x = "job_katn"),
  function(x, workers = 5, space = glue::glue("copykat_batch_{x@sig}"), cl = 5, ...)
  {
    step_message("Running...")
    x$space <- space
    dir.create(space, FALSE)
    if (is.remote(x)) {
      x <- set_remote_for_sub_jobs(x, space)
      object(x) <- pbapply::pblapply(object(x), cl = cl, 
        function(obj) {
          obj <- step1(obj, workers = workers, ...)
        })
    } else {
      paths <- file.path(space, names(object(x)))
      lapply(paths,
        function(path) {
          if (dir.exists(path)) {
            stop(glue::glue("{path} exists."))
          }
        })
      n <- 0L
      object(x) <- pbapply::pblapply(object(x),
        function(obj) {
          n <<- n + 1L
          wd <- getwd()
          writeLines("")
          cli::cli_alert_info(glue::glue("Running {n}: {wd}"))
          obj <- try(step1(obj, workers, paths[n]))
          if (inherits(obj, "try-error")) {
            setwd(wd)
            return(obj)
          }
          obj <- step2(obj)
          saveRDS(obj, file.path(path, "res.rds"))
          return(obj)
        })
    }
    x <- snapAdd(x, "以 `CopyKAT` 鉴定恶质细胞。")
    x <- methodAdd(x, "R 包 `CopyKAT` ({packageVersion('copykat')}) 用于鉴定恶性细胞 {cite_show('DelineatingCopGaoR2021')}。`CopyKAT` 可以区分整倍体与非整倍体，其中非整倍体被认为是肿瘤细胞，而整倍体是正常细胞 {cite_show('CausesAndConsGordon2012')}。由于 `CopyKAT` 不适用于多样本数据 (批次效应的存在) ，因此，对各个样本独立鉴定。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_katn"),
  function(x, workers = 20, cl = 5, ignore = FALSE, ...){
    step_message("Plot heatmap.")
    if (is.remote(x)) {
      object(x) <- pbapply::pblapply(object(x), cl = cl, 
        function(obj) {
          if (ignore) {
            obj@step <- 1L
          }
          obj <- step2(
            obj, workers = workers, inherits = TRUE, 
            ignore = ignore, ...
          )
        })
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_katn"),
  function(x){
    step_message("Collate results (remote).")
    if (is.remote(x)) {
      file_backup <- file.path(x$space, "objects.rds")
      if (!file.exists(file_backup)) {
        saveRDS(object(x), file_backup)
      }
      object(x) <- lapply(object(x), 
        function(obj) {
          obj$res_copykat <- NULL
          return(obj)
        })
      plots <- lapply(object(x), function(x) x@plots$step2$p.copykat)
      lab(plots) <- glue::glue("{sig(x)} all malignant cells heatmap")
      tables <- lapply(object(x), function(x) x@tables$step2$res_copykat)
      res <- rbind_list(tables, "orig.ident")
      res <- set_lab_legend(
        res,
        glue::glue("{x@sig} all copyKAT prediction data"),
        glue::glue("为 copyKAT 所有样本注释结果附表。")
      )
      p.props <- plot_cells_proportion(
        res, "copykat.pred", "orig.ident", FALSE
      )
      p.props <- set_lab_legend(
        p.props,
        glue::glue("{x@sig} proportions of aneuploid and diploid"),
        glue::glue("为 copyKAT 注释的所有样本中 aneuploid (Malignant cell) 与
          diploid (Benign cell) 的细胞比例。")
      )
      x <- plotsAdd(x, all_heatmap = plots, p.props = p.props)
      x <- tablesAdd(x, t.res_copykat = res)
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_katn"),
  function(x, ref, from = "scsa_cell", to = "copykat_cell")
  {
    if (ref@step < 3L) {
      stop('ref@step < 3L.')
    }
    if (is.null(res <- ref@tables$step3$t.res_copykat)) {
      stop('is.null(res <- ref@tables$step3$t.res_copykat).')
    }
    fun_method <- selectMethod("map", c("job_seurat", "job_kat"))
    ref <- .job_kat(tables = list(step2 = list(res_copykat = res)))
    fun_method(x, ref)
  })

setMethod("set_remote", signature = c(x = "job_katn"),
  function(x, wd = glue::glue("~/katn_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

