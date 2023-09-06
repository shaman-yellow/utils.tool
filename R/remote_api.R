# ==========================================================================
# for remote files operation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rem_list.files <- function(path, pattern, ...) {
  if (!check_remote()) {
    if (length(path) > 1) {
      sapply(path, list.files, pattern = pattern, simplify = F)
    } else {
      list.files(path, pattern, ...)
    }
  } else {
    list.remote(path, pattern,
      remote = x@params$remote, wd = x@params$wd
    )
  }
}

rem_readLines <- function(...) {
  if (!check_remote()) {
    readLines(...)
  } else {
    stop("...")
  }
}

rem_file.copy <- function(...) {
  if (!check_remote()) {
    file.copy(...)
  } else {
    stop("...")
  }
}

rem_file.rename <- function(...) {
  if (!check_remote()) {
    file.rename(...)
  } else {
    stop("...")
  }
}

rem_file.remove <- function(...) {
  if (!check_remote()) {
    file.remove(...)
  } else {
    stop("...")
  }
}

rem_file.exists <- function(...) {
  if (!check_remote()) {
    file.exists(...)
  } else {
    stop("...")
  }
}

rem_normalizePath <- function(...) {
  if (!check_remote()) {
    normalizePath(...)
  } else {
    stop("...")
  }
}

rem_unlink <- function(...) {
  if (!check_remote()) {
    unlink(...)
  } else {
    stop("...")
  }
}

rem_list.files <- function(...) {
  if (!check_remote()) {
    list.files(...)
  } else {
    stop("...")
  }
}

rem_dir.create <- function(...) {
  if (!check_remote()) {
    dir.create(...)
  } else {
    stop("...")
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
