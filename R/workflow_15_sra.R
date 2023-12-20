# ==========================================================================
# workflow of sra
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_sra <- setClass("job_sra", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("...")
    ))

job_sra <- function(project, wd = "sra_data")
{
  dir.create(wd)
  x <- .job_sra(object = project)
  x@params$wd <- wd
  x
}

setMethod("step0", signature = c(x = "job_sra"),
  function(x){
    step_message("Prepare your data with function `job_sra`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_sra"),
  function(x){
    step_message("Prepare for dowloading SRA data (16s rRNA data).")
    cli::cli_alert_info("esearch -db sra -query")
    cdRun("esearch -db sra -query ", x@object,
        " | efetch -format runinfo > ", x@params$wd, "/info.csv")
    x@params$info <- ftibble(paste0(x@params$wd, "/info.csv"))
    return(x)
  })

setMethod("step2", signature = c(x = "job_sra"),
  function(x){
    step_message("Download load SRA data.")
    cli::cli_alert_info("bash: prefetch ...")
    pbapply::pblapply(x@params$info$Run,
        function(id) {
          path <- x@params$wd
          exists <- file.exists(paste0(path, "/", id))
          while (!exists) {
            cdRun("prefetch ", id, " ", "--output-directory ", path)
            exists <- file.exists(paste0(path, "/", id))
          }
        })
    x@params$all_sra <- list.files(x@params$wd, "\\.sra$", full.names = T, recursive = T)
    return(x)
  })

setMethod("step3", signature = c(x = "job_sra"),
  function(x, pattern = "fastq\\.gz$"){
    step_message("Format as fastq file")
    x@params$pattern <- pattern
    cli::cli_alert_info("fastq-dump --gzip --split-3")
    pbapply::pblapply(x@params$all_sra,
        function(file) {
          path <- get_path(file)
          if (length(list.files(path, pattern)) == 0)
            cdRun("fastq-dump --gzip --split-3 ", file, " -O ", path)
        })
    return(x)
  })

setMethod("step4", signature = c(x = "job_sra"),
  function(x, filter = T){
    step_message("Try to format `x@params$info` as metadata.")
    info <- dplyr::mutate(x@params$info, SampleName = gs(SampleName, "_", "."))
    filepath <- lapply(info$Run,
      function(id) {
        path <- list.files(paste0(x@params$wd, "/", id), x@params$pattern, full.names = T)
        normalizePath(path)
      })
    metadata <- dplyr::select(info, "sample-id" = SampleName, Run)
    metadata <- try_fqs_meta(metadata, filepath, filter = filter)
    x@params$metadata <- metadata
    print(metadata)
    return(x)
  })

setGeneric("asjob_qiime", 
  function(x, ...) standardGeneric("asjob_qiime"))

setMethod("asjob_qiime", signature = c(x = "job_sra"),
  function(x, wd = "qiime_data"){
    dir.create(wd)
    job_qiime(x@params$metadata, wd)
  })

setGeneric("asjob_sra", 
  function(x, ...) standardGeneric("asjob_sra"))

setMethod("asjob_sra", signature = c(x = "df"),
  function(x, path){
    metadata <- x
    x <- .job_sra()
    x@step <- 3L
    x@params$info <- metadata
    x@params$wd <- path
    x <- step4(x)
    return(x)
  })
