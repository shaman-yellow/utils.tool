# ==========================================================================
# workflow of fastp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_fastp <- setClass("job_fastp", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/OpenGene/fastp")
    ))

job_fastp <- function(path)
{
  .job_fastp(object = path)
}

setMethod("step0", signature = c(x = "job_fastp"),
  function(x){
    step_message("Prepare your data with function `job_fastp`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_fastp"),
  function(x, f1_r2 = c("f1", "r2"), suffix = ".fastq.gz", workers = 7){
    step_message("Quality control (QC).
      "
    )
    if (!is.null(x@params$set_remote)) {
      fqs <- list.remote(object(x), paste0(suffix, "$"), x@params$remote)
    } else {
      fqs <- list.files(object(x), paste0(suffix, "$"), full.names = T, recursive = T)
    }
    fqs_path <- unique(get_path(fqs))
    fqs_path <- fqs_path[!grepl("fastp_qc$", fqs_path)]
    if (!is.null(x@params$set_remote)) {
      pbapply::pblapply(fqs_path, fastp_pair.remote,
        f1_r2 = f1_r2, suffix = suffix, workers = workers, x = x)
    } else {
      pbapply::pblapply(fqs_path, fastp_pair, f1_r2 = f1_r2,
        suffix = suffix, workers = workers)
    }
    x@params$fqs_path <- fqs_path
    return(x)
  })

setMethod("step2", signature = c(x = "job_fastp"),
  function(x, pattern_report = "html", pattern_fq = "QC.fastq.gz$"){
    step_message("Collate results which succeed with report as metadata.")
    if (is.remote(x)) {
      reports <- list.remote(object(x), pattern_report)
    } else {
      reports <- list.files(object(x), pattern_report, full.names = T, recursive = T)
    }
    dirs <- get_path(get_path(reports))
    metadata <- data.frame(SampleName = get_realname(gs(dirs, "/$", "")),
      dirs = dirs, reports = reports)
    metadata <- mutate(metadata, Run = SampleName)
    if (is.remote(x)) {
      filepath <- list.remote(dirs, pattern_fq)
    } else {
      filepath <- lapply(dirs,
        function(dir) {
          list.files(dir, pattern_fq, full.names = T, recursive = T)
        })
    }
    metadata <- try_fqs_meta(metadata, filepath)
    metadata <- as_tibble(metadata)
    failed <- x@params$fqs_path[!x@params$fqs_path %in% metadata$dirs]
    if (length(failed) > 0) {
      x@params$failed <- failed
      message("Failed:")
      print(failed)
    }
    x@params$metadata <- metadata
    if (is.remote(x)) {
      if (!file.exists(to <- x@params$map_local)) {
        dir.create(to)
      }
      cp.remote(metadata$reports, to)
    }
    return(x)
  })

set_remote.default <- function(x, tmpdir, map_local, remote) {
  x@params$remote <- remote
  x@params$set_remote <- T
  x@params$map_local <- map_local
  x@params$tmpdir <- tmpdir
  return(x)
}

setMethod("set_remote", signature = c(x = "job_fastp"),
  function(x, tmpdir = "/data/hlc/tmp", map_local = "fastp_local", remote = "remote")
  {
    set_remote.default(x, tmpdir, map_local, remote)
    x@params$postfix <- function(x) {
      x[1] <- gs(x[1], "^fastp", "~/miniconda3/bin/conda run -n base fastp")
      x
    }
    return(x)
  })

cp.remote <- function(from, to, remote = "remote") {
  if (missing(remote)) {
    x <- get("x", envir = parent.frame(1))
    remote <- x@params$remote
  }
  pbapply::pblapply(from,
    function(path) {
      system(paste0("scp ", remote, ":", path, " ", to))
    })
}

list.remote <- function(path, pattern, remote = "remote") {
  if (missing(remote)) {
    x <- get("x", envir = parent.frame(1))
    remote <- x@params$remote
  }
  if (length(path) == 1) {
    files <- system(paste0("ssh ", remote, " 'find ", path, "'"), intern = T)
    files[ grepl(pattern, files) ]
  } else if (length(path) > 1) {
    files <- system(paste0("ssh ", remote, " ",
        "'for i in ", paste0(path, collapse = " "),
        "; do find $i; echo -----; done'"),
      intern = T)
    files <- sep_list(files, "^-----$")
    files <- lapply(files,
      function(files) {
        files <- files[ -length(files) ]
        files[ grepl(pattern, files) ]
      })
    names(files) <- path
    files
  } else {
    stop("The path may be character(0).")
  }
}

fastp_pair <- function(path, suffix = ".fastq.gz", pattern = paste0("[1-2]", suffix, "$"),
  f1_r2 = c("R1", "R2"),
  names = unique(gs(list.files(path, pattern), pattern, "")), workers = 4,
  copy_report = T, cdRun = get_fun("cdRun"))
{
  dir.create(paste0(path, "/fastp_qc"), F)
  dir.create(paste0(path, "/fastp_report"), F)
  pbapply::pblapply(names,
    function(name) {
      i1 <- paste0(name, f1_r2[1], suffix)
      i2 <- paste0(name, f2_r2[2], suffix)
      o1 <- paste0("fastp_qc/", name, "1.QC.fastq.gz")
      o2 <- paste0("fastp_qc/", name, "2.QC.fastq.gz")
      report <- paste0("fastp_report/", name, ".html")
      if (!file.exists(paste0(path, "/", report))) {
        cdRun("conda run -n base fastp -i ", i1,
          " -I ", i2,
          " -o ", o1,
          " -O ", o2,
          " -h ", report,
          " --detect_adapter_for_pe ",
          " -w ", workers,
          path = path
        )
      }
    })
  if (copy_report) {
    file.copy(paste0(path, "/fastp_report"), ".", recursive = T)
  }
  message("Job finished.")
}

fastp_pair.remote <- function(path, suffix = ".fastq.gz", pattern = paste0(".[1-2]", suffix, "$"),
  f1_r2 = c("R1", "R2"), workers = 4, copy_report = T, cdRun = get_fun("cdRun"), x)
{
  names <- unique(gs(list.remote(path, pattern), pattern, ""))
  expr <- paste0("mkdir fastp_qc fastp_report;")
  lapply(names,
    function(name) {
      name <- get_realname(name)
      i1 <- paste0(name, f1_r2[1], suffix)
      i2 <- paste0(name, f1_r2[2], suffix)
      o1 <- paste0("fastp_qc/", name, "1.QC.fastq.gz")
      o2 <- paste0("fastp_qc/", name, "2.QC.fastq.gz")
      report <- paste0("fastp_report/", name, ".html")
      if (!remote_file.exists(paste0(path, "/", report))) {
        remoteRun("fastp -i ", i1,
          " -I ", i2,
          " -o ", o1,
          " -O ", o2,
          " -h ", report,
          " --detect_adapter_for_pe ",
          " -w ", workers,
          run_after_cd = expr,
          path = path
        )
      }
    })
  message("Job finished.")
}

setMethod("asjob_qiime", signature = c(x = "job_fastp"),
  function(x, metadata = x@params$metadata, wd = "qiime_data"){
    .check_columns(metadata, c("SampleName", "group"))
    metadata <- rename(metadata, `sample-id` = SampleName)
    if (!file.exists(wd)) {
      dir.create(wd)
    }
    x <- job_qiime(metadata, wd)
    x@params$pattern_fq <- "QC.fastq.gz$"
    x
  })
