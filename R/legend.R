# ==========================================================================
# for summary text files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sumText <- function(num, simplify = T, files = paste0("ch", num, ".md")) {
  sapply(files, simplify = simplify,
    function(file) {
      md <- paste0(readLines(file), collapse = "\n")
      md <- stringr::str_extract_all(md, "[\u4e00-\u9fa5]{1,}")[[ 1 ]]
      md <- paste0(md, collapse = "")
      nchar(md)
    })
}

