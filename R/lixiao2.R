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

gtitle <- function(grob, title, fill = "#E18727FF") {
  into(grecti3(title, tfill = fill), grob)
}

# ==========================================================================
# combine pdf picture
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.fig_frame <- setClass("fig_frame", 
  contains = c(),
  representation = representation(file = "character",
    width = "numeric", height = "numeric",
    panel.width = "numeric"),
  prototype = prototype(
    panel.width = .95
    ))

.figs_frame <- setClass("figs_frame", 
  contains = c("fig_frame"),
  representation = representation(
    figs_width = "numeric"),
  prototype = prototype(
    panel.width = .88
    ))

.column <- setClass("column", 
  contains = c("list"),
  representation = representation(
    rel.width = "numeric", rel.height = "numeric"),
  prototype = NULL)

setValidity("column", 
  function(object){
    allif <- vapply(object, function(x) is(x, "fig_frame"), logical(1))
    if (all(allif)) T else stop("Object must be 'fig_frame'.")
  })

.columns <- setClass("columns", 
  contains = c("list"),
  representation = representation(
    rel.width = "numeric", rel.height = "numeric",
    rel.cl.widths = "numeric"),
  prototype = NULL)

setValidity("columns", 
  function(object){
    allif <- vapply(object, function(x) is(x, "column"), logical(1))
    if (all(allif)) T else stop("Object must be 'column'.")
  })

setMethod("show", signature = c(object = "column"),
  function(object){
    message("Object: ", length(object), "\n",
      "Relative width: ", object@rel.width, "\n",
      "Relative height: ", round(object@rel.height, 2)
    )
  })

setMethod("show", signature = c(object = "columns"),
  function(object){
    message("Object: ", length(object), "\n",
      "Relative width: ", object@rel.width, "\n",
      "Relative height: ", round(object@rel.height, 2), "\n",
      "Relative column width: ", 
      paste0(round(object@rel.cl.widths, 2), collapse = ", ")
    )
  })

setMethod("show", signature = c(object = "fig_frame"),
  function(object){
    message(object@file, "\n",
      "Width: ", object@width, "\nHeight: ", round(object@height, 2)
    )
  })

setMethod("show", signature = c(object = "figs_frame"),
  function(object){
    message(paste0(object@file, collapse = ", "), "\n",
      "Width: ", object@width, "\nHeight: ", round(object@height, 2),
      "\nFigures width: ",
      paste0(round(object@figs_width, 2), collapse = ",")
    )
  })

new_fig_frame <- function(file) {
  if (!file.exists(file))
    stop("file.exists(file) == F")
  info <- pdftools::pdf_pagesize(file)
  if (nrow(info) > 1)
    stop("Only one page pdf file supported.")
  .fig_frame(file = file, width = 1, height = info$height / info$width)
}

new_figs_frame <- function(...) {
  figs <- lapply(list(...), new_fig_frame)
  files <- vapply(figs, function(x) x@file, character(1))
  heights <- vapply(figs, function(x) x@height, numeric(1))
  scale.widths <- 1 / heights
  rel.widths <- scale.widths / sum(scale.widths) - .01
  height <- rel.widths[1] * heights[1]
  .figs_frame(file = files, width = 1, height = height,
    figs_width = rel.widths)
}

rw <- new_figs_frame

cl <- function(...) {
  lst <- lapply(list(...),
    function(x) {
      if (is(x, "character"))
        new_fig_frame(x)
      else if (is(x, "figs_frame"))
        x
      else stop("Target can not be added in.")
    })
  sum.height <- sum(vapply(lst, function(x) x@height, numeric(1)))
  .column(lst, rel.width = 1, rel.height = sum.height)
}

cls <- function(...) {
  lst <- list(...)
  cl.heights <- vapply(lst, function(x) x@rel.height, numeric(1))
  cl.widths <- 1 / cl.heights
  rel.cl.widths <- cl.widths / sum(cl.widths)
  rel.height <- rel.cl.widths[1] * cl.heights[1]
  .columns(lst, rel.width = 1, rel.height = rel.height, 
    rel.cl.widths = rel.cl.widths
  )
}

setGeneric("render", 
  function(x, ...) standardGeneric("render"))

setMethod("render", signature = c(x = "columns"),
  function(x, name, engine = "xelatex"){
    lines <- astex(x)
    if (missing(name))
      name <- paste0(substitute(x, parent.frame(1)), ".tex")
    writeLines(lines, name)
    system(paste0(engine, " ", name))
  })

setMethod("render", signature = c(x = "column"),
  function(x, name, engine = "xelatex"){
    if (missing(name))
      name <- paste0(substitute(x, parent.frame(1)), ".tex")
    x <- cls(x)
    render(x, name, engine)
  })

setGeneric("astex", 
  function(x, ...) standardGeneric("astex"))

setMethod("astex", signature = c(x = "columns"),
  function(x, width = 20, fun_head = pictureMergeHead){
    head <- fun_head(width, width * x@rel.height)
    contents <- unlist(mapply(astex, x, x@rel.cl.widths, SIMPLIFY = F))
    c(head, "", "\\begin{document}", "",
      "\\begin{figure}[htb]\\centering", "",
      contents, "",
      "\\label{fig:main}",
      "\\end{figure}",
      "", "\\end{document}")
  })

setMethod("astex", signature = c(x = "column"),
  function(x, width){
    contents <- unlist(lapply(x, astex))
    c(paste0("\\begin{col}{", width, "\\textwidth}"), contents,
    "\\end{col}")
  })

setMethod("astex", signature = c(x = "fig_frame"),
  function(x){
    c("\\begin{minipage}[t]{.95\\linewidth}",
      paste0("\\sidesubfloat[]{\\includegraphics[width=", x@panel.width,
        "\\textwidth]{", x@file, "}}"),
      "\\end{minipage}"
    )
  })

setMethod("astex", signature = c(x = "figs_frame"),
  function(x){
    lines <- lapply(1:length(x@file), 
      function(n) {
        c("\\begin{minipage}[t]{", x@figs_width[n], "\\linewidth}",
          paste0("\\sidesubfloat[]{\\includegraphics[width=", x@panel.width,
            "\\textwidth]{", x@file[n], "}}"),
          "\\end{minipage}"
        )
      })
    c("\\begin{minipage}[t]{.98\\linewidth}",
      "", unlist(lines), "",
      "\\end{minipage}")
  })

pictureMergeHead <- function(width = 20, height = 20) {
  head <- readLines(paste0(.expath, "/pictureMergeHead.tex"))
  page <- paste0("\\geometry{paperheight = ", height, "cm, paperwidth = ", width, "cm, ",
    "left = 3mm, right = 3mm, top = 3mm, bottom = 3mm}")
  c(head, page)
}


