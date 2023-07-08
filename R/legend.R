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

format_chname <- function(dir, pattern = "ch[0-9]{1,}\\.md", num = 3){
  files <- list.files(dir, pattern, full.names = T)
  names <- get_filename(files)
  count <- nchar(stringr::str_extract(names, "[0-9]{1,}"))
  mapply(files, count, USE.NAMES = F,
    FUN = function(file, c){
      if (c <= num) {
        fix <- paste0(rep("0", num - c), collapse = "")
        new <- gsub("ch([0-9]{1,})\\.md", paste0("ch", fix, "\\1.md"), file)
        file.rename(file, new)
      }
    })
}

