# ==========================================================================
# workflow of iobr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_iobr <- setClass("job_iobr",
  contains = c("job"),
  prototype = prototype(
    pg = "iobr",
    info = c("https://iobr.github.io/book/index.html"),
    cite = "[@Enhancing_immun_Zeng_2024]",
    method = "",
    tag = "iobr",
    analysis = "IOBR 肿瘤免疫微环境分析"
    ))

job_iobr <- function(object, metadata)
{
  x <- .job_iobr(object = object)
  x$metadata <- metadata
  return(x)
}

setGeneric("asjob_iobr", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_iobr"))

setMethod("asjob_iobr", signature = c(x = "job_limma"),
  function(x, use = .guess_symbol(x))
  {
    snap <- snap(x)
    object <- extract_unique_genes.job_limma(x)
    mtx <- object$E
    rownames(mtx) <- gname(object$genes[[ use ]])
    x <- snapAdd(
      job_iobr(mtx, metadata = object$targets), 
      "以数据集 ({x$project}, dataset: {x@sig}) 进行 IOBR 免疫微环境评估。"
    )
    x <- snapAdd(x, snap)
    return(x)
  })

setMethod("step0", signature = c(x = "job_iobr"),
  function(x){
    step_message("Prepare your data with function `job_iobr`.")
  })

setMethod("step1", signature = c(x = "job_iobr"),
  function(x, method = c("cibersort"))
  {
    step_message("Calculating ...")
    method <- match.arg(method)
    require(IOBR)
    if (is.null(x$res)) {
      x$res <- e(IOBR::deconvo_tme(
        eset = object(x), method = method, arrays = FALSE
      ))
    }
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_iobr"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
