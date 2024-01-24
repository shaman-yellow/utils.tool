# ==========================================================================
# workflow of bowtie2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_bowtie2 <- setClass("job_bowtie2", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/BenLangmead/bowtie2"),
    cite = "[@ScalingReadAlLangme2019; @FastGappedReaLangme2012; @UltrafastAndMLangme2009]",
    method = "The tool of `Bowtie2` was used for alignment of DNA sequences to the human genome"
    ))

job_bowtie2 <- function(metadata, wd = ".", workers = 10, command = "bowtie2")
{
  .job_bowtie2(object = metadata, params = list(wd = wd, workers = workers), pg = command)
}

setMethod("step0", signature = c(x = "job_bowtie2"),
  function(x){
    step_message("Prepare your data with function `job_bowtie2`.")
  })

setMethod("step1", signature = c(x = "job_bowtie2"),
  function(x, path_ref = "../ref", file_ref = "hg38.fa", file_index = gs(file_ref, "\\.fa$", ""))
  {
    step_message("Prepare index file")
    x$ref <- paste0(path_ref, "/", file_ref)
    if (!rem_file.exists(x$ref)) {
      rem_run("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz -P ", path_ref)
      rem_run("gunzip ", path_ref, "hg38.fa.gz")
    }
    x$index <- paste0(path_ref, "/", file_index)
    if (!rem_file.exists(paste0(x$index, ".1.bt2l"))) {
      rem_run(pg(x), "-build --large-index ",
        " --threads ", x$workers,
        " ", x$ref,
        " ", x$index)
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_bowtie2"),
  function(x, mode = c("notHost", "host")){
    step_message("Alignment to reference genome.")
    mode <- match.arg(mode)
    res <- pbapply::pbapply(object(x), 1,
      function(vec) {
        name <- get_realname(vec[[ "forward-absolute-filepath" ]])
        notHost <- paste0(name, "_rmHost")
        output <- paste0(name, ".bam")
        notHost_file <- paste0(notHost, ".", 1:2, ".fq.gz")
        if (!any(rem_file.exists(c(output, notHost_file)))) {
          rem_run(pg(x),
            " -x ", x$index,
            " -1 ", vec[[ "forward-absolute-filepath" ]],
            " -2 ", vec[[ "reverse-absolute-filepath" ]],
            if (mode == "notHost") {
              paste0(" --un-conc-gz ", notHost)
            } else NULL,
            " --threads ", x$workers,
            " | ", pg("samtools", is.remote(x)), " view -Sb - > ", output
            # " -S ", output
          )
        }
        if (mode == "notHost") {
          vec[[ "f_notHost" ]] <- notHost_file[1]
          vec[[ "r_notHost" ]] <- notHost_file[2]
        } else {
          vec[[ "host" ]] <- output
        }
      })
    x$res <- res
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_bowtie2"),
  function(x, wd)
  {
    x$wd <- wd
    x$set_remote <- T
    return(x)
  })

setGeneric("asjob_bowtie2", 
  function(x, ...) standardGeneric("asjob_bowtie2"))

setMethod("asjob_bowtie2", signature = c(x = "job_fastp"),
  function(x, workers = 10){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    wd <- object(x)
    x <- .job_bowtie2(object = x$metadata, params = x@params)
    x@pg <- "bowtie2"
    x$wd <- wd
    x$workers <- workers
    return(x)
  })
