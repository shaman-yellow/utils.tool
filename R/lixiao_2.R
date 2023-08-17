# ==========================================================================
# extend of lixiao.R
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

kall_index <- function(file, index = "../hg38_mrna.idx") {
  cdRun("kallisto index",
    " -i ", index,
    " ", file)
}

kall_quant <- function(path, pattern = "[1-2]\\.QC\\.fastq\\.gz$",
  names = unique(gs(list.files(path, pattern), pattern, "")), workers = 7,
  index = "~/outline/lixiao/hg38_mrna.idx", output_dir = paste0("quant_", get_realname(index)),
  copy_report = T)
{
  if (!file.exists(index)) {
    stop("Please download the cDNA reference file, and use `kall_index` to build the index file.")
  }
  if (!length(names)) {
    stop("length(names) == 0")
  }
  dir.create(paste0(path, "/", output_dir), F)
  pbapply::pblapply(names,
    function(name) {
      i1 <- paste0(name, "1.QC.fastq.gz")
      i2 <- paste0(name, "2.QC.fastq.gz")
      output <- paste0(output_dir, "/", name)
      cdRun("kallisto quant -i ", index,
        " -o ", output,
        " ", i1, " ", i2,
        " -t ", workers,
        path = path
      )
    })
  if (copy_report) {
    file.copy(paste0(path, "/", output_dir), ".", recursive = T)
  }
  message("Job finished.")
}

read_kall_quant <- function(path) {
  files <- list.files(path, "abundance.tsv$", full.names = T, recursive = T)
  metadata <- data.frame(file = files)
  metadata <- mutate(metadata, directory = get_path(file),
    sample = get_realname(directory)
  )
  n <- 0
  counts <- lapply(metadata$file,
    function(file) {
      n <<- n + 1
      data <- ftibble(file)
      data <- mutate(data, count = est_counts)
      if (n == 1)
        select(data, target_id, count)
      else
        select(data, count)
    })
  counts <- suppressMessages(do.call(dplyr::bind_cols, counts))
  counts <- mutate(counts, target_id = gs(target_id, "\\.[0-9]{1,}$", ""))
  counts <- distinct(counts, target_id, .keep_all = T)
  colnames(counts)[-1] <- metadata$sample
  namel(counts, metadata)
}

fastp_pair <- function(path, pattern = "[1-2]\\.fastq\\.gz$",
  names = unique(gs(list.files(path, pattern), pattern, "")), workers = 4,
  copy_report = T)
{
  dir.create(paste0(path, "/fastp_qc"), F)
  gir.create(paste0(path, "/fastp_report"), F)
  pbapply::pblapply(names,
    function(name) {
      i1 <- paste0(name, "1.fastq.gz")
      i2 <- paste0(name, "2.fastq.gz")
      o1 <- paste0("fastp_qc/", name, "1.QC.fastq.gz")
      o2 <- paste0("fastp_qc/", name, "2.QC.fastq.gz")
      report <- paste0("fastp_report/", name, ".html")
      cdRun("fastp -i ", i1,
        " -I ", i2,
        " -o ", o1,
        " -O ", o2,
        " -h ", report,
        " --detect_adapter_for_pe ",
        " -w ", workers,
        path = path
      )
    })
  if (copy_report) {
    file.copy(paste0(path, "/fastp_report"), ".", recursive = T)
  }
  message("Job finished.")
}

.fasta <- setClass("fasta", 
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

setMethod("show", 
  signature = c(object = "fasta"),
  function(object){
    print(stringr::str_trunc(head(object, n = 10), 40))
    cat("...\n")
  })

setGeneric("group", 
  function(x, ...) standardGeneric("group"))

setMethod("group", signature = c(x = "fasta"),
  function(x, each = 500){
    n <- -1
    group <- vapply(x, FUN.VALUE = numeric(1),
      function(x) {
        n <<- n + 1
        n %/% 2
      })
    x <- unname(split(x, group))
    x <- grouping_vec2list(x, each, T)
    x <- lapply(x,
      function(x) {
        .fasta(unlist(x))
      })
    x
  })

setGeneric("write", 
  function(x, ...) standardGeneric("write"))

setMethod("write", signature = c(x = "fasta"),
  function(x, max = 500){
    dir.create("fasta", F)
    name <- paste0("fasta/", as.character(substitute(x, parent.frame(1))))
    n <- 0
    lapply(group(x),
      function(x) {
        n <<- n + 1
        file <- paste0(name, "_", n, ".fasta")
        writeLines(gs(x, "\\*$", ""), file)
      })
  })

setMethod("write", signature = c(x = "grob.obj"),
  function(x, name, width, height, dev = "png", factor = 100){
    file <- paste0(get_savedir("figs"), "/", name, ".", dev)
    match.fun(dev)(file, width * factor, height * factor)
    grid.draw(x)
    dev.off()
  })

get_seq.pro <- function(ids, mart, unique = T, fasta = T, from = "hgnc_symbol", to = "peptide") {
  ids <- unique(ids)
  data <- e(biomaRt::getSequence(id = ids, type = from, seqType = to, mart = mart))
  data <- filter(data, !grepl("unavailable", !!rlang::sym(to)))
  if (unique) {
    data <- mutate(data, long = nchar(!!rlang::sym(to)))
    data <- arrange(data, dplyr::desc(long))
    data <- distinct(data, hgnc_symbol, .keep_all = T)
  }
  if (fasta) {
    fasta <- apply(data, 1,
      function(vec) {
        c(paste0(">", vec[[2]]), vec[[1]])
      })
    fasta <- .fasta(unlist(fasta))
    namel(data, fasta)
  } else {
    data
  }
}

get_seq.rna <- function(ids, mart, unique = T, fasta = T, from = "hgnc_symbol", to = "cdna") {
  do.call(get_seq.pro, as.list(environment()))
}
