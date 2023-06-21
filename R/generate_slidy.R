# ==========================================================================
# slidy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

generate_slidy <- function(name, path = "/mnt/data/wizard/Documents/zcmu_reports")
{
  file <- paste0(path, "/", name, ".Rmd")
  meta <- generate_slidy_meta()
  cat(meta, file = file)
}

make_slidy <- function(name, path = "/mnt/data/wizard/Documents/zcmu_reports")
{
  file <- paste0(path, "/", name, ".Rmd")
  rmarkdown::render(file)
}

preprocess_bib <- function(file = paste0(.expath, "/library.bib")) {
  lst <- read_bib(file)
  which <- grepl("^@[0-9]*$", names(lst))
  lst[which] <- lapply(lst[which],
    function(l) {
      l[1] <- gsub("\\{", "{FIXN", l[1])
      return(l)
    })
  writeLines(unlist(lst), "library.bib")
}

reload_bib <- function(file = "library.bib",
  exclude = c("doi", "urldate", "issn", "address", "isbn"))
{
  bib <- bibtex::read.bib(file)
  bib <- fix_bib(bib, exclude)
  assign(".bib", bib, envir = topenv())
}

fix_bib <- function(bib, exclude = c("doi", "urldate", "issn", "address", "isbn"))
{
  bib <- lapply(bib,
    function(b) {
      expr <- paste0(paste0("b$", exclude, " <- NULL"), collapse = "\n")
      eval(parse(text = expr))
      return(b)
    })
  bib <- do.call(c, bib)
  return(bib)
}

citethis <- function(..., trunc = T, trunc.author = 1, trunc.title = F,
  prefix = "\\tiny ", sep = "\\vspace{0.5em} \\newline\n",
  pkgs = NULL, pkgs.fix = fix_bib, exbibentry = NULL, postFun = NULL)
{
  keys <- list(...)
  if (length(keys) > 0) {
    keys <- vapply(keys,
      function(ch) {
        if (is.numeric(ch) | grepl("^[0-9]*$", ch)) {
          return(paste0("FIXN", ch))
        } else {
          ch
        }
      }, character(1))
    keys <- paste0(keys, ".", keys)
    if (!exists(".bib", where = topenv())) {
      reload_bib()
    }
    all <- names(.bib)
    bib <- .bib[ match(keys, all) ]
  } else {
    bib <- bibentry()
  }
  if (!is.null(pkgs)) {
    bib.pkg <- lapply(pkgs,
      function(pkg) {
        c(citation(pkg))
      })
    bib.pkg <- do.call(c, bib.pkg)
    if (!is.null(pkgs.fix)) {
      bib.pkg <- pkgs.fix(bib.pkg)
    }
  } else {
    bib.pkg <- bibentry()
  }
  bib <- c(bib, bib.pkg)
  if (!is.null(exbibentry)) {
    bib <- c(bib, exbibentry)
  }
  if (trunc) {
    bib <- lapply(bib,
      function(b) {
        if (length(b$author) > trunc.author) {
          b$author <- c(b$author[1:trunc.author], person("et al."))
        }
        if (trunc.title) {
          b$title <- stringr::str_trunc(gsub("\\{|\\}", "", b$title), 20)
        } else {
          b$title <- "#";
        }
        b$number <- NULL
        b$volume <- NULL;
        b$pages <- NULL
        return(b)
      })
    bib <- do.call(c, bib)
  }
  if (!is.null(bib))
    writeLines(prefix)
  text <- format(bib, "text")
  if (!trunc.title) {
    text <- gsub("[\"“”#]\\.*", "", text)
  }
  if (!is.null(postFun)) {
    text <- vapply(text, postFun, character(1))
  }
  writeLines(text, sep = sep)
}

testSection <- function(file, pattern, level = 2, render = rmarkdown::render) {
  lines <- readLines(file)
  yml <- getyml(lines)
  getSecPos <- function(lines, level) {
    pattern <- paste0("^", paste0("#{", level, "}[^#]"))
    grep(pattern, lines)
  }
  sec.pos <- getSecPos(lines, level)
  tar.pos <- grep(pattern, lines)
  if (length(tar.pos) == 0) {
    stop("length(tar.pos) == 0")
  }
  if (length(tar.pos) > 1) {
    warning("length(tar.pos) > 1")
    tar.pos <- tar.pos[1]
  }
  if (tar.pos < (yend <- attr(yml, "pos")$end)) {
    stop("tar.pos < attr(yml, 'pos')$end")
  }
  if (any(sec.pos == tar.pos)) {
    start.pos <- tar.pos
    if (level > 1) sec.pos <- getSecPos(lines, paste0("1,", level))
    end.pos <- sec.pos[head(which(sec.pos > tar.pos), n = 1)]
  } else {
    start.pos <- sec.pos[tail(which(sec.pos < tar.pos), n = 1)]
    if (level > 1) sec.pos <- getSecPos(lines, paste0("1,", level))
    end.pos <- sec.pos[head(which(sec.pos > tar.pos), n = 1)]
  }
  end.pos <- end.pos - 1
  if (length(start.pos) == 0) {
    start.pos <- yend + 1
  }
  if (start.pos > length(lines)) {
    stop("start.pos > length(lines)")
  }
  if (length(end.pos) == 0) {
    end.pos <- length(lines)
  }
  if (end.pos < start.pos) {
    stop("end.pos < start.pos")
  }
  content <- lines[ start.pos:end.pos ]
  lines <- c(yml, "", content)
  writeLines(lines, nfile <- paste0("_temptest_", get_filename(file)))
  if (!is.null(render)) {
    output <- render(nfile)
    op(output)
  }
}

getyml <- function(lines) {
  pos <- grep("^---", lines)
  if (length(pos) != 2) {
    stop("length(pos) != 2")
  }
  yml <- lines[pos[1]:pos[2]]
  attr(yml, "pos") <- list(start = pos[1], end = pos[2])
  yml
}

arrange_figsPath <- function(file, pattern = "^!\\[.*\\]", to = NULL, overwrite = F)
{
  path <- get_path(file)
  filename <- get_filename(file)
  if (!is.null(to)) {
    if (!dir.exists(to)) {
      stop("dir.exists(to) == F")
    } else {
      dir <- to
    }
  } else {
    dir <- paste0(path, "/", gsub("\\.[a-zA-Z]*$", "", filename))
  }
  lines <- readLines(file)
  pos <- grep(pattern, lines)
  fun <- function(n) {
    line <- lines[n]
    file <- gsub("\"", "", stringr::str_extract(line, "(?<=\\]\\().*(?=\\))"))
    filename <- get_filename(file)
    nfile <- paste0(dir, "/", filename)
    line <- gsub("\\]\\(.*\\)", paste0("](", nfile, ")"), line)
    file.copy(file, dir, overwrite)
    line
  }
  for (i in pos) {
    lines[i] <- fun(i)
  }
  writeLines(lines, file)
}

latexfig <- function(file, caption = NULL, scale.width = 0.7, scale.height = 0.45)
{
  if (grepl("\\.pdf$", file)) {
    info <- pdftools::pdf_pagesize(file)
  } else {
    info <- bitmap_info(file)
  }
  ratio <- info$width / info$height
  width <- normSize(ratio, list(width = scale.width, height = scale.height))$width
  md <- c(paste0("::: {.col data-latex=\"{", round(width, 2), "\\textwidth}\"}"),
    paste0("![", caption, "](", file, ")"),
    ":::")
  writeLines(md)
}

hlName <- function(text, name = "Huang L") {
  gsub(paste0("(", name, ")"), "\\\\underline{\\\\textbf{\\1}}", text)
}

normSize <- function(ratio, scale)
{
  if ((scale$width / scale$height) >= ratio) {
    ## height as reference
    height <- scale$height
    width <- height * ratio
  } else {
    ## width as reference
    width <- scale$width
    height <- width / ratio
  }
  list(height = height, width = width)
}

bitmap_info <- function(file) {
  img <- magick::image_read(file)
  info <- magick::image_info(img)
  magick::image_destroy(img)
  info
}

get_title <- function(file) {
  lines <- readLines(file)
  pos <- grepl("^#", lines)
  lines[pos]
} 
