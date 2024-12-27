# ==========================================================================
# workflow of mn2kat
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mn2kat <- setClass("job_mn2kat", 
  contains = c("job_monocle"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: ...")
    ))

setMethod("do_monocle", signature = c(x = "job_seurat", ref = "job_kat"),
  function(x, ref, dims = 1:15, resolution = 1.2){
    x <- getsub(x, cells = grp(x@object@meta.data$scsa_copykat, "Cancer"))
    x@step <- 2L
    sr_cancer <- step3(x, dims, resolution)
    mn <- asjob_monocle(sr_cancer, "seurat_clusters")
    x <- .job_mn2kat()
    x@params <- mn@params
    x@object <- mn@object
    x$sr_cancer <- sr_cancer
    x$kat <- ref
    return(x)
  })

setMethod("step0", signature = c(x = "job_mn2kat"),
  function(x){
    step_message("Prepare your data with function `job_mn2kat`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_mn2kat"),
  function(x, ...){
    x <- callNextMethod(x, ...)
    p.cancer_position <- suppressMessages(map(x$sr_cancer, x$kat, cutree = TRUE))
    x@plots[[ 1 ]] <- c(x@plots[[ 1 ]], namel(p.cancer_position))
    return(x)
  })

setMethod("map", signature = c(x = "job_mn2kat"),
  function(x, cutree){
    map(x$sr_cancer, x$kat, cutree = cutree)
  })

setMethod("vis", signature = c(x = "job_mn2kat"),
  function(x, group.by = "seurat_clusters"){
    vis(x$sr_cancer, "seurat_clusters")
  })

setMethod("asjob_seurat", signature = c(x = "job_mn2kat"),
  function(x, k, rename = NULL){
    x <- regroup(x$sr_cancer, x$cellClass_tree.gene_module, k, rename = rename)
    show(vis(x, "regroup.hclust", 1.5))
    return(x)
  })
