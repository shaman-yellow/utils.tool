# ==========================================================================
# workflow of kall
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_kall <- setClass("job_kall", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("Tutorial: http://pachterlab.github.io/kallisto/manual.html",
      "\nhttps://ftp.ensembl.org/pub/release-110/fasta/"),
    cite = "[@NearOptimalPrBray2016]",
    method = "`Kallisto` used for RNA-seq mapping and quantification"
    ))

setGeneric("asjob_kall", 
  function(x, ...) standardGeneric("asjob_kall"))

setMethod("asjob_kall", signature = c(x = "job_fastp"),
  function(x){
    metadata <- x$metadata
    x <- .job_kall()
    x$metadata <- metadata
    return(x)
  })

setMethod("step0", signature = c(x = "job_kall"),
  function(x){
    step_message("Prepare your data with function `asjob_kall`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_kall"),
  function(x, idx = .prefix("hg38_mrna.idx", "db"), ref_file = .prefix("Homo_sapiens.GRCh38.cdna.all.fa.gz", "db"))
    {
    step_message("Prepare kallisto gene reference index file.")
    if (!file.exists(idx)) {
      kall_index(ref_file, idx)
    }
    x$idx <- normalizePath(idx)
    return(x)
  })

setMethod("step2", signature = c(x = "job_kall"),
  function(x, workers = 8, output = "kallisto_quantification")
  {
    step_message("Mapping to reference genome.")
    pbapply::pbapply(x$metadata, 1,
      function(info) {
        path <- get_path(info[[ "forward-absolute-filepath" ]])
        dir.create(normed_output <- paste0(path, "/", output), F)
        i1 <- get_filename(info[[ "forward-absolute-filepath" ]])
        i2 <- get_filename(info[[ "reverse-absolute-filepath" ]])
        dir <- paste0(output, "/", get_realname(info[[ "forward-absolute-filepath" ]]))
        if (!dir.exists(normed_dir <- paste0(path, "/", dir))) {
          cdRun("kallisto quant -i ", x$idx,
            " -o ", dir,
            " ", i1, " ", i2,
            " -t ", workers,
            path = path
          )
        }
        file.copy(normed_output, ".", recursive = T)
      })
    x$output <- output
    return(x)
  })

setMethod("step3", signature = c(x = "job_kall"),
  function(x, path = x$output){
    step_message("Collate all quantification results.")
    res <- read_kall_quant(path)
    x@tables[[ 3 ]] <- res
    return(x)
  })
