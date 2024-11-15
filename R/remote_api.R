# ==========================================================================
# for remote files operation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rem_run <- function(...) {
  x <- get("x", envir = parent.frame(1))
  if (is.null(x$wd)) {
    message("Use `x$wd` as '.'")
    x$wd <- "."
  }
  if (!check_remote()) {
    cdRun(..., path = x$wd)
  } else {
    remoteRun(..., path = x$wd, tmpdir = x$tmpdir,
      remote = x$remote, postfix = x$postfix,
      run_after_cd = x$run_after_cd, x = x
    )
  }
}

remoteRun <- function(..., path, run_after_cd = NULL,
  postfix = NULL, remote = "remote", tmpdir = NULL, x)
{
  expr <- paste0(unlist(list(...)), collapse = "")
  if (missing(x)) {
    x <- get("x", parent.frame(1))
  }
  if (missing(remote)) {
    remote <- x@params$remote
  }
  if (missing(postfix)) {
    postfix <- x@params$postfix
  }
  if (!is.null(postfix)) {
    expr <- postfix(expr)
  }
  if (missing(tmpdir)) {
    tmpdir <- x@params$tmpdir
  }
  if (!is.null(tmpdir))
    expr <- c(paste0("export TMPDIR=", tmpdir), expr)
  if (missing(path)) {
    path <- x@params$wd
  }
  if (!is.null(run_after_cd)) {
    expr <- c(run_after_cd, expr)
  }
  if (!is.null(path)) {
    expr <- c(paste0("cd ", path), expr)
  }
  heading <- getOption("scriptHeading")
  if (!is.null(heading)) {
    expr <- c(heading, "", expr)
  }
  script <- tempfile("remote_Script_", fileext = ".sh")
  writeLines(expr, script)
  writeLines(crayon::yellow(paste0("The script file for remote is: ", script)))
  remote_script <- file.path(path, basename(script))
  system(paste0("scp ", script, " ", remote, ":", remote_script))
  cmd <- paste0("ssh ", remote, " '", getOption("remoteRun", "bash"), " ", remote_script, "'")
  cli::cli_alert_info(cmd)
  system(cmd)
}

set_remoteRun <- function(remoteRun = "bash", scriptHeading = NULL, scriptPrefix = NULL) {
  options(remoteRun = remoteRun,
    scriptHeading = scriptHeading,
    scriptPrefix = scriptPrefix
  )
}

scriptPrefix <- function(x) {
  if (!is.remote(x)) {
    return()
  }
  prefix <- getOption("scriptPrefix")
  if (is.null(prefix)) {
    return(NULL)
  } else {
    return(paste0(prefix, " "))
  }
}

rem_list.files <- function(path, pattern,
  all.files = F, full.names = F, recursive = F, x)
{
  if (!check_remote()) {
    if (length(path) > 1) {
      sapply(path, list.files, pattern = pattern, simplify = F)
    } else {
      list.files(path, pattern, all.files = all.files,
        full.names = full.names, recursive = recursive)
    }
  } else {
    if (missing(x)) {
      x <- get("x", envir = parent.frame(1))
    }
    list.remote(path, pattern, all.files = all.files,
      full.names = full.names, recursive = recursive,
      remote = x@params$remote, x = x
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

rem_get <- function(file) {
  if (!check_remote()) {
    file
  } else {
    x <- get("x", envir = parent.frame(1))
    dir.create(x$map_local, F)
    local <- paste0(x$map_local, "/", file)
    get_file_from_remote(file, x$wd, local)
    local
  }
}

get_file_from_remote <- function(file, wd, to = NULL, remote = "remote")
{
  if (is.null(to)) {
    to <- tempfile()
  }
  cdRun("scp ", remote, ":", wd, "/", file, " ", to)
  to
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
    res <- system(paste0("ssh ", x$remote, " '",
        "cd ", x$wd, "; ",
        paste0(expr_sys.file.exists(file), collapse = "; "), "'"),
      intern = T)
    ifelse(res == "T", T, F)
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

check_remote <- function(n = 2) {
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

list.remote <- function(path, pattern, remote = "remote",
  all.files = F, full.names = F, recursive = F, x)
{
  if (missing(x))
    x <- get("x", envir = parent.frame(1))
  if (missing(remote)) {
    remote <- x@params$remote
  }
  options <- " -mindepth 1 "
  if (!recursive) {
    options <- paste0(options, " -maxdepth 1 ")
  }
  if (!all.files) {
    options <- paste0(options, " -not -name \".*\" ")
  }
  before <- paste0("cd ", x@params$wd, ";")
  if (length(path) == 1) {
    if (!full.names) {
      before <- paste0(before, " cd ", path, ";")
      path <- "."
    }
    files <- system(paste0("ssh ", remote, " '", before,
        " find ", path, " ", options, "'"),
      intern = T)
    files[ grepl(pattern, basename(files)) ]
  } else if (length(path) > 1) {
    if (!full.names) {
      files <- system(paste0("ssh ", remote, " ",
          "'", before, " for i in ", paste0(path, collapse = " "),
          "; do cd $i; find . ", options, "; echo -----; done'"),
        intern = T)
    } else {
      files <- system(paste0("ssh ", remote, " ",
          "'", before, " for i in ", paste0(path, collapse = " "),
          "; do find $i ", options, "; echo -----; done'"),
        intern = T)
    }
    files <- sep_list(files, "^-----$")
    files <- lapply(files,
      function(files) {
        files <- files[ -length(files) ]
        files[ grepl(pattern, basename(files)) ]
      })
    names(files) <- path
    files
  } else {
    stop("The path may be character(0).")
  }
}
