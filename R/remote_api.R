# ==========================================================================
# for remote files operation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rem_list.files <- function(path, pattern,
  all.files = F, full.names = F, recursive = F)
{
  if (!check_remote()) {
    if (length(path) > 1) {
      sapply(path, list.files, pattern = pattern, simplify = F)
    } else {
      list.files(path, pattern, all.files = all.files,
        full.names = full.names, recursive = recursive)
    }
  } else {
    x <- get("x", envir = parent.frame(1))
    list.remote(path, pattern, all.files = all.files,
      full.names = full.names, recursive = recursive,
      remote = x@params$remote, wd = x@params$wd
    )
  }
}

rem_readLines <- function(file, ...) {
  if (!check_remote()) {
    readLines(file, ...)
  } else {
    x <- get("x", envir = parent.frame(1))
    file <- get_file_from_remote(file, x$wd, x$remote)
    readLines(file, ...)
  }
}

get_file_from_remote <- function(file, wd, remote = "remote")
{
  tempfile <-  tempfile()
  cdRun("scp ", remote, ":", wd, "/", file, " ", tempfile)
  tempfile
}

rem_ftibble <- function(file, ...) {
  if (!check_remote()) {
    ftibble(file, ...)
  } else {
    x <- get("x", envir = parent.frame(1))
    file <- get_file_from_remote(file, x$wd, x$remote)
    ftibble(file, ...)
  }
}

rem_file.copy <- function(from, to, recursive = F, ...) {
  if (!check_remote()) {
    file.copy(from, to, recursive = recursive, ...)
  } else {
    x <- get("x", envir = parent.frame(1))
    if (recursive) {
      options <- "-r"
    } else {
      options <- ""
    }
    cdRun("ssh ", x$remote, " '",
      "cd ", x$wd, "; ",
      "if [ -d ", to, " ]; then arg=-t; fi;",
      "cp ", options, " ", from, " $arg ", to,
      "'")
  }
}

rem_file.rename <- function(from, to) {
  if (!check_remote()) {
    file.rename(from, to)
  } else {
    x <- get("x", envir = parent.frame(1))
    cdRun("ssh ", x$remote, " '",
      "cd ", x$wd, "; ",
      "mv ", from, " $arg ", to,
      "'")
  }
}

rem_file.remove <- function(...) {
  if (!check_remote()) {
    file.remove(...)
  } else {
    x <- get("x", envir = parent.frame(1))
    cdRun("ssh ", x$remote, " '",
      "cd ", x$wd, "; ",
      "rm ", paste0(unlist(list(...)), collapse = " "),
      "'")
  }
}

rem_file.exists <- function(file) {
  if (!check_remote()) {
    file.exists(file)
  } else {
    x <- get("x", envir = parent.frame(1))
    res <- system(paste0("ssh ", remote, " '",
        "cd ", x$wd, "; ",
        expr_sys.file.exists(file), "'"),
      intern = T)
    if (res == "T") T
    else F
  }
}

rem_normalizePath <- function(...) {
  if (!check_remote()) {
    normalizePath(...)
  } else {
    x <- get("x", envir = parent.frame(1))
    system(paste0("ssh ", x$remote, " '",
        "cd ", x$wd, "; ",
        "readlink -f ", paste0(unlist(list(...)), collapse = " "),
        "'"), intern = T)
  }
}

rem_unlink <- function(...) {
  if (!check_remote()) {
    unlink(...)
  } else {
    x <- get("x", envir = parent.frame(1))
    stop("...")
  }
}

rem_dir.create <- function(path, ...) {
  if (!check_remote()) {
    dir.create(path, ...)
  } else {
    x <- get("x", envir = parent.frame(1))
    cdRun("ssh ", x$remote, " '",
      "cd ", x$wd, "; ",
      "mkdir ", path,
      "'")
  }
}

check_remote <- function(n = 1) {
  remote <- F
  if (exists("x", envir = parent.frame(n))) {
    x <- get("x", envir = parent.frame(n))
    if (is(x, "job")) {
      if (is.remote(x)) {
        remote <- T
      }
    }
  }
  remote
}

set_if_null <- function(object, value) {
  if (is.null(object)) 
    value
  else
    object
}
