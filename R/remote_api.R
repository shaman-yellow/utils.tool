# ==========================================================================
# for remote files operation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

rem_run <- function(..., .script = NULL, .append = FALSE, .test = FALSE) {
  x <- get("x", envir = parent.frame(1))
  if (is.null(x$wd)) {
    message("Use `x$wd` as '.'")
    x$wd <- "."
  }
  if (!check_remote(1L) && !.test) {
    cdRun(..., path = x$wd)
  } else {
    remoteRun(..., path = x$wd, tmpdir = x$tmpdir,
      remote = x$remote, postfix = x$postfix,
      run_after_cd = x$run_after_cd, x = x, .script = .script, 
      .append = .append, .test = .test
    )
  }
}

remoteRun <- function(..., path, run_after_cd = NULL,
  postfix = NULL, remote = "remote", tmpdir = NULL, 
  .script = NULL, .append = FALSE, .test = FALSE, x)
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
  if (is.null(.script)) {
    script <- tempfile("remote_Script_", fileext = ".sh")
  } else {
    script <- .script
  }
  if (.append) {
    cat(expr, "\n\n", sep = "\n", file = script, append = TRUE)
  } else {
    writeLines(expr, script)
  }
  writeLines(crayon::yellow(paste0("The script file for remote is: ", script)))
  if (!.append && !.test) {
    .run_script_in_remote(script, path, remote = remote)
  }
  script
}

.run_script_in_remote <- function(script, path,
  remote_script = file.path(path, basename(script)), remote = "remote")
{
  system(paste0("scp ", script, " ", remote, ":", remote_script))
  cmd <- paste0("ssh ", remote, " '", getOption("remoteRun", "bash"), " ", remote_script, "'")
  cli::cli_alert_info(cmd)
  system(cmd)
}

set_remoteRun.bosai <- function(core = 16) {
  set_remoteRun(glue::glue("sbatch -p v6_384 -N 1 -n 1 -c {core}"), "#!/bin/bash")
}

set_remoteRun <- function(remoteRun = "bash", scriptHeading = NULL, scriptPrefix = NULL) {
  options(remoteRun = remoteRun,
    scriptHeading = scriptHeading,
    scriptPrefix = scriptPrefix
  )
}

testRem_file.exists <- function(x, file, wait = 10,
  cancel = "testRem", env_cancel = .GlobalEnv, later = TRUE, show_queque = TRUE)
{
  getFun <- function() {
    res <- try(get(cancel, envir = env_cancel), TRUE)
    if (inherits(res, "try-error") | !is.logical(res)) {
      res <- TRUE
    }
    res
  }
  testFun <- function() {
    if ((notHasThat <- !rem_file.exists(file)) && getFun()) {
      if (show_queque) {
        system(glue::glue("ssh {x$remote} 'squeue'"))
      }
      if (later) {
        message(
          glue::glue(
            "'!rem_file.exists(file), try again in {wait} minutes,\nUse `{cancel} <- FALSE` to cancel."
          )
        )
        later::later(testFun, wait * 60)
      } else {
        message(glue::glue("'!rem_file.exists(file), try again in {wait} minutes."))
        Sys.sleep(wait * 60)
        testFun()
      }
    } else if (!notHasThat && Sys.which("notify-send") != "") {
      cdRun("notify-send 'Got File: ", x$wd, "/", file, "'")
    }
  }
  testFun()
}

scriptPrefix <- function(x) {
  if (!is.remote(x)) {
    return()
  }
  prefix <- getOption("scriptPrefix")
  if (is.null(prefix)) {
    return(NULL)
  } else {
    return(paste0(prefix, "\n\n"))
  }
}

rem_list.files <- function(path, pattern,
  all.files = FALSE, full.names = FALSE, recursive = FALSE, x)
{
  if (!check_remote(1L)) {
    if (length(path) > 1) {
      sapply(path, list.files, pattern = pattern, simplify = FALSE)
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
    file <- get_file_from_remote(file, x$wd, remote = x$remote)
    readLines(file, ...)
  }
}

rem_get <- function(file) {
  if (!check_remote()) {
    file
  } else {
    x <- get("x", envir = parent.frame(1))
    dir.create(x$map_local, FALSE)
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

rem_file.copy <- function(from, to, recursive = FALSE, ...) {
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

rem_file.exists <- function(file, wd) {
  if (!check_remote()) {
    file.exists(file)
  } else {
    x <- get("x", envir = parent.frame(1))
    if (missing(wd)) {
      wd <- x$wd
    }
    res <- system(paste0("ssh ", x$remote, " '",
        "cd ", wd, "; ",
        paste0(expr_sys.file.exists(file), collapse = "; "), "'"),
      intern = TRUE)
    ifelse(res == "TRUE", TRUE, FALSE)
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
        "'"), intern = TRUE)
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

rem_dir.create <- function(path, ..., wd) {
  if (!check_remote()) {
    dir.create(path, ...)
  } else {
    x <- get("x", envir = parent.frame(1))
    if (missing(wd)) {
      wd <- x$wd
    }
    cdRun("ssh ", x$remote, " '",
      "cd ", wd, "; ",
      "mkdir ", path,
      "'")
  }
}

check_remote <- function(n = 2) {
  remote <- FALSE
  if (exists("x", envir = parent.frame(n))) {
    x <- get("x", envir = parent.frame(n))
    if (is(x, "job")) {
      if (is.remote(x)) {
        remote <- TRUE
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
  all.files = FALSE, full.names = FALSE, recursive = FALSE, x)
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
      intern = TRUE)
    files[ grepl(pattern, basename(files)) ]
  } else if (length(path) > 1) {
    if (!full.names) {
      files <- system(paste0("ssh ", remote, " ",
          "'", before, " for i in ", paste0(path, collapse = " "),
          "; do cd $i; find . ", options, "; echo -----; done'"),
        intern = TRUE)
    } else {
      files <- system(paste0("ssh ", remote, " ",
          "'", before, " for i in ", paste0(path, collapse = " "),
          "; do find $i ", options, "; echo -----; done'"),
        intern = TRUE)
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

.args_transformer <- function(text, envir) {
  res <- glue::identity_transformer(text, envir)
  if (length(res) != 1) {
    stop('.args_transformer, length(res) != 1')
  }
  if (is.character(res)) {
    paste0("'", res, "'")
  } else {
    res
  }
}

run_job_remote <- function(x, expression, ..., envir = parent.frame(1), wd = x$wd,
  name = attr(x@sig, 'name'), rds = glue::glue("{wd}/{name}.rds"),
  ignore_remote_cache = FALSE, ignore_local_cache = FALSE,
  inherit_last_result = FALSE, step_just = TRUE, wait = 1L)
{
  # prepare script
  if (!is.remote(x)) {
    stop('!is.remote(x).')
  }
  if (is.null(name)) {
    stop('is.null(name). step < 1L?')
  }
  expression <- as.expression(substitute(expression))
  # glue with "{arg}", if character, -> 'test' (with quote), else if is numeric, -> 1 (without quote)
  fun_glue <- function(expression, envir = parent.frame(1)) {
    glue::glue(as.character(expression),
      .transformer = .args_transformer,
      .open = "\"{", .close = "}\"",
      .envir = envir)
  }
  # set a file for judging whether finished.
  file_end <- glue::glue("{name}_{x@step}.finish")
  try_continue <- FALSE
  if (rem_file.exists(file_end)) {
    if (!ignore_remote_cache) {
      try_continue <- TRUE
    }
    rem_file.remove(file_end)
  }
  finish <- fun_glue(expression({
    saveRDS(x, "{rds}")
    writeLines("END", "{file_end}")
  }))
  file_upload <- add_filename_suffix(rds, "UPLOAD")
  setup <- fun_glue(expression({
    devtools::load_all("~/utils.tool")
    if ("{inherit_last_result}" || "{try_continue}") {
      x <- readRDS("{rds}")
      if ("{try_continue}") {
        if (if_job_step_expected(x, "{x@step}")) {
          writeLines("END", "{file_end}")
          quit("no")
        } else {
          x <- readRDS("{file_upload}")
        }
      }
    } else {
      x <- readRDS("{file_upload}")
    }
  }))
  main <- fun_glue(expression, envir)
  script <- shQuote(paste0(c(setup, main, finish), collapse = "\n"))
  hash <- digest::digest(c(x@sig, script))
  # prepare data (job)
  if (is.null(x$map_local)) {
    stop('is.null(x$map_local).')
  }
  dir.create(x$map_local, FALSE)
  file_local <- file.path(x$map_local, basename(rds))
  file_res <- add_filename_suffix(file_local, hash)
  if (!file.exists(file_res) || ignore_local_cache) {
    if (!inherit_last_result) {
      if (step_just) {
        x@step <- x@step - 1L
      }
      # `not_remote`, Prevent recursion
      message(glue::glue("Saving file in local: {file_local}"))
      saveRDS(not_remote(x), file_local)
      if (step_just) {
        x@step <- x@step + 1L
      }
      cdRun(glue::glue("scp {file_local} {x$remote}:{x$wd}/{basename(file_upload)}"))
    } else {
      message(crayon::red(glue::glue("`inherit_last_result` == TRUE, use previous result file.")))
    }
    rem_run(glue::glue("{pg('Rscript', TRUE)} -e {script}"), ...)
    # wait for finish (later = FALSE)
    testRem_file.exists(x, file_end, wait, later = FALSE)
    get_file_from_remote(basename(rds), x$wd, file_res, x$remote)
  } else {
    message(glue::glue("File already exists: {file_res}"))
  }
  res <- as_remote(readRDS(file_res))
  if (res@step != x@step) {
    warning(crayon::red('res@step != x@step.'))
  } else {
    x <- res
  }
  return(x)
}

if_job_step_expected <- function(x, expect) {
  if (x@step == expect) {
    message("Job step has been finished previously.")
    return(TRUE)
  } else {
    return(FALSE)
  }
}

not_remote <- function(x) {
  x$set_remote <- FALSE
  return(x)
}

as_remote <- function(x) {
  x$set_remote <- TRUE
  return(x)
}

.install_basic_packages <- function() {
  packages <- c(
    "aplot", "bibtex", "cli", "ComplexHeatmap", "crayon", "data.table",
    "devtools", "dplyr", "flextable", "ggplotify", "ggraph", "ggrepel",
    "ggsci", "ggVennDiagram", "glmnet", "Hmisc", "knitr",
    "magick", "magrittr", "Matrix", "openxlsx", "igraph", "patchwork", "pbapply",
    "pdftools", "png", "qpdf", "R.utils", "RCurl", "readxl", "rlang",
    "rsvg", "stringr", "UpSetR", "usethis", "XML", "xml2"
  )
  packages <- lapply(packages, 
    function(x) {
      if (requireNamespace(x, quietly = TRUE)) {
        x
      } else NULL
    })
  packages <- lst_clear0(packages)
  if (length(packages)) {
    options("download.file.method" = "wget")
    BiocManager::install(unlist(packages))
  }
}

