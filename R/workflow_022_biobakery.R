# ==========================================================================
# workflow of biobakery
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_biobakery <- setClass("job_biobakery", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "biobakery_workflows",
    info = paste0("http://huttenhower.sph.harvard.edu/biobakery_workflows",
      "\nhttps://github.com/biobakery/biobakery/wiki/biobakery_workflows#2-metagenome-profiling"),
    cite = "[@BiobakeryAMeMciver2017]",
    method = "`Biobakery` used for metagenome preprocessing",
    tag = "metagenome",
    analysis = "Biobakery 宏基因组数据前处理"
    ))

setGeneric("asjob_biobakery", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_biobakery"))

setMethod("asjob_biobakery", signature = c(x = "job_fastp"),
  function(x){
    y <- .job_biobakery()
    y$wd <- object(x)
    y$metadata <- x$metadata
    y$from_job_fastp <- TRUE
    return(y)
  })

setMethod("step0", signature = c(x = "job_biobakery"),
  function(x){
    step_message("Prepare your data with function `job_biobakery`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_biobakery"),
  function(x, workers = 9){
    step_message("One-stop analysis.")
    if (!is.null(x$from_job_fastp)) {
      options <- " --bypass-quality-control "
      pair_id <- "_1.QC."
    } else {
      stop("This workflow only accept job converted from job_fastp.")
    }
    rem_run(pg(x), " wmgx ",
      " --input ",
      " --remove-intermediate-output ",
      " --output res_biobakery ",
      " --pair-identifier ", pair_id,
      " --threads ", workers,
      " ", options
    )
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_biobakery"),
  function(x, wd, postfix = NULL, db = "/data/hlc/biobakery_workflows_databases"){
    x$postfix <- postfix
    x$run_after_cd <- paste0("export KNEADDATA_DB_HUMAN_GENOME=", db)
    x$wd <- wd
    return(x)
  })
