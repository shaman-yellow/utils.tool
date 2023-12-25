# ==========================================================================
# description for conception (colnames)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.get_des <- function(ref) {
  db <- .get_db_des()
  lst <- db[ names(db) %in% ref ]
  vapply(names(lst), function(x) lst[[ x ]], character(1))
}

.get_db_des <- function(fresh = F) {
  if (is.null(maybe <- getOption("des")) | fresh) {
    des <- .parse_md()
    options(des = des)
    return(des)
  }
  maybe
}

.parse_md <- function(dir = "~/outline/lixiao/", pattern = "^help.*md$|description.*md$",
  excludes = c("id", ".id"))
{
  files <- list.files(dir, pattern, full.names = T, recursive = T)
  lst <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (grpl(name <- get_realname(file), "description")) {
        lines <- lines[ grpl(lines, "^[0-9\\-].*:") ]
        lines <- gs(lines, "^-\\s*|`", "")
      } else if (grpl(name, "help")) {
        lines <- sep_list(lines, "^$")
        lines <- unlist(lapply(lines, function(x) paste0(x, collapse = " ")))
        lines <- gs(lines, "^\\s*|\\s*$", "")
        lines <- gs(lines, "\\s{2,}", " ")
      }
      lines <- strsplit(lines, ":")
      lines
    })
  lst <- unlist(lst, recursive = F)
  lst <- lapply(lst, function(x) if (length(x) > 1) x else NULL)
  lst <- lst_clear0(lst)
  names(lst) <- NULL
  lst <- nl(vapply(lst, function(x) x[1], character(1)),
    vapply(lst, function(x) x[2], character(1)))
  lst <- lst[ !duplicated(names(lst)) ]
  lst <- lst[ !names(lst) %in% excludes ]
  lapply(lst, function(x) gs(x, "([{}])", "\\\\\\1"))
}
