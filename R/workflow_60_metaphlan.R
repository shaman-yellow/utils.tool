# ==========================================================================
# workflow of mpa
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mpa <- setClass("job_mpa", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/biobakery/MetaPhlAn/wiki/MetaPhlAn-4"),
    cite = "[@ExtendingAndIBlanco2023]",
    method = "`MetaPhlAn` used for profiling the composition of microbial communities from metagenomic data"
    ))

setMethod("step0", signature = c(x = "job_mpa"),
  function(x){
    step_message("Prepare your data with function `asjob_mpa`.")
  })

setMethod("step1", signature = c(x = "job_mpa"),
  function(x, path_db = "../bowtie2db")
  {
    step_message("Dowload (prepare) the database files.")
    if (!rem_file.exists(path_db)) {
      rem_dir.create(path_db)
    }
    x$path_db <- path_db
    latest <- paste0(path_db, "/mpa_latest")
    url <- "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/"
    if (!rem_file.exists(latest)) {
      rem_run("wget ", url, "mpa_latest", " -P ", path_db)
    }
    version <- readLines(get_file_from_remote(latest, x$wd, x$remote))
    files <- rem_list.files(path_db, version)
    if (!length(files)) {
      thats <- paste0(url, "/", version, c(".md5", ".tar", "_marker_info.txt.bz2", "_species.txt.bz2"))
      indexs <- paste0(url, "/bowtie2_indexes/", version, c("_bt2.md5", "_bt2.tar"))
      alls <- c(thats, indexs)
      pbapply::pblapply(alls,
        function(url) {
          rem_run("wget -c ", url, " -P ", path_db)
        })
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_mpa"),
  function(x) {
    step_message("Mapping... Time consumed.")
    res <- pbapply::pbapply(x$metadata, 1, simplify = T,
      function(vec) {
        name <- get_realname(vec[[ "forward-absolute-filepath" ]])
        output <- paste0(name, "_metagenome.txt")
        if (!rem_file.exists(output)) {
          rem_run(pg(x), " ",
            vec[[ "forward-absolute-filepath" ]], ",",
            vec[[ "reverse-absolute-filepath" ]],
            " --bowtie2out ", name, ".bowtie2.bz2",
            " --nproc ", x$workers,
            " --input_type fastq",
            " --bowtie2db ", x$path_db,
            " -o ", output
          )
        }
        vec[[ "mpa_output" ]] <- output
        vec
      })
    x$res <- res
    return(x)
  })

setMethod("step3", signature = c(x = "job_mpa"),
  function(x){
    rem_run(
      "merge_metaphlan_tables.py",
      " ", paste0(files),
      "> output/merged_abundance_table.txt")
  })

setMethod("set_remote", signature = c(x = "job_mpa"),
  function(x, wd, postfix = NULL, run_after_cd = NULL, tmpdir = NULL){
    x$wd <- wd
    x$set_remote <- T
    return(x)
  })


setGeneric("asjob_mpa", 
  function(x, ...) standardGeneric("asjob_mpa"))

setMethod("asjob_mpa", signature = c(x = "job_fastp"),
  function(x, workers = 10){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    wd <- object(x)
    x <- .job_mpa(object = x$metadata, params = x@params)
    x@pg <- "metaphlan"
    x$wd <- wd
    x$workers <- workers
    return(x)
  })

