# ==========================================================================
# for ZFvim ...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.date <- function() {
  res <- cat(as.character(Sys.Date()))
}

read_dic <- function(file) {
  db <- readLines(file)
  db.i <- stringr::str_extract(db, "^[a-z]*")
  tibble::tibble(n = seq_along(db), i = db.i, db = db)
}
