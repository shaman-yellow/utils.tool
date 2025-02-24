# ==========================================================================
# workflow of nscfea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_nscfea <- setClass("job_nscfea", 
  contains = c("job"),
  prototype = prototype(
    pg = "nscfea",
    info = c("https://github.com/changwn/scFEA/tree/master"),
    cite = "[@AGraphNeuralAlgham2021]",
    method = "",
    tag = "scrna:flux",
    analysis = "scFEA 单细胞数据的代谢通量预测"
    ))

setGeneric("asjob_nscfea", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_nscfea"))

setMethod("asjob_nscfea", signature = c(x = "job_seurat5n"),
  function(x, org = c("mouse", "human"), sample = "orig.ident", 
    dir = "nscfea", testWhich = NULL, ...)
  {
    org <- match.arg(org)
    meta <- x@object@meta.data
    group_cells <- split(seq_len(nrow(meta)), meta[[ sample ]])
    if (!is.null(testWhich)) {
      group_cells <- group_cells[ testWhich ]
    }
    if (file.exists(dir)) {
      if (sureThat("Remove exists `dir`?")) {
        unlink(dir, TRUE)
        dir.create(dir, FALSE)
      }
    } else {
      dir.create(dir, FALSE)
    }
    dirs <- file.path(
      dir, paste0(seq_along(group_cells), "_scfea")
    )
    objects <- mapply(group_cells, dirs, SIMPLIFY = FALSE,
      FUN = function(cells, dir) {
        asjob_scfea(x, cells, dir = dir, org = org, ...)
      })
    x <- .job_nscfea(object = objects)
    x <- methodAdd(x, objects[[1]])
    x <- snapAdd(x, "将 `Seurat` 对象中的各个样本分别以 `scFEA` 预测代谢通量 (避免批次效应)。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_nscfea"),
  function(x){
    step_message("Prepare your data with function `job_nscfea`.")
  })

setMethod("step1", signature = c(x = "job_nscfea"),
  function(x){
    step_message("Run remote in batch.")
    if (!is.remote(x)) {
      stop('!is.remote(x).')
    }
    object(x) <- lapply(object(x), 
      function(object) {
        step1(object)
      })
    return(x)
  })

setMethod("step2", signature = c(x = "job_nscfea"),
  function(x){
    step_message("Collate results")
    object(x) <- lapply(object(x), 
      function(object) {
        step2(object)
      })
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_nscfea"),
  function(x, wd = "nscfea")
  {
    x$wd <- "."
    if (rem_file.exists(wd)) {
      isThat <- sureThat("Dir exists, remove all files ?")
      if (isThat) {
        cdRun("ssh ", x$remote, " 'rm -r ", wd, "'")
      }
    }
    rem_dir.create(wd)
    x$wd <- wd
    wds <- file.path(wd, paste0(seq_along(object(x)), "_scfea"))
    object(x) <- mapply(object(x), wds, SIMPLIFY = FALSE,
      FUN = function(object, wd) {
        set_remote(object, wd)
      })
    return(x)
  })
