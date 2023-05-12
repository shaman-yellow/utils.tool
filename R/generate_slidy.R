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
  bib <- lapply(bib,
    function(b) {
      expr <- paste0(paste0("b$", exclude, " <- NULL"), collapse = "\n")
      eval(parse(text = expr))
      return(b)
    })
  bib <- do.call(c, bib)
  assign(".bib", bib, envir = topenv())
}

citethis <- function(..., trunc = T, prefix = "\\tiny ", sep = "\\vspace{0.5em} \\newline ") {
  keys <- list(...)
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
  if (trunc) {
    bib <- lapply(bib,
      function(b) {
        if (length(b$author) > 3) {
          b$author <- c(b$author[1:3], person("et al."))
        }
        return(b)
      })
    bib <- do.call(c, bib)
  }
  knitr::opts_current$set(echo = F, eval = T, results = "asis")
  writeLines(prefix)
  writeLines(c(format(bib, "text")), sep = sep)
}

