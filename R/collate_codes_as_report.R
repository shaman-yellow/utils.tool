# ==========================================================================
# for codes as report output
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rough_records <- function(dir, rdname, title, files= NULL,
  pattern = "\\.R$", sep = "^# =*\\s*$")
{
  codes <- get_codes.dir(dir, files = files, pattern = pattern)
  chunks <- as_chunk.codes(codes, sep = sep)
  lines <- c("---", "---", "", chunks)
  report <- as_report.rough(lines)
  # write_thesisDocx(report, "codes_of_mcnebula2Docx.Rmd", "R codes of MCnebula2")
  write_articlePdf(report, rdname, title)
}

get_codes.dir <- function(dir, files = NULL, pattern = "\\.R$") {
  if (is.null(files)) {
    files <- list.files(dir, full.names = T, recursive = T)
  }
  files <- files[grepl(pattern, files)]
  codes <- sapply(files, readLines, simplify = F)
  names(codes) <- vapply(names(codes), get_filename, character(1))
  codes
}

as_chunk.codes <- function(lst, as_lines = T, sep = "^# =*\\s*$") {
  name <- names(lst)
  lines <- lst
  fun <- function(lines) {
    lst <- sep_list(lines, sep, before = T)
    lst <- lapply(lst,
      function(ch) {
        c('```{r eval = F, echo = T}', ch, '```', "")
      })
    unlist(lst)
  }
  chunks <- lapply(1:length(name),
    function(n) {
      c(paste0("# File: ", name[n]), "", fun(lines[[ n ]]))
    })
  if (as_lines) {
    chunks <- unlist(chunks)
  }
  chunks
}


