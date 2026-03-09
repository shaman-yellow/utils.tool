# ==========================================================================
# huibang
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setup.hb <- function(project = guess_project(), ws = getRemoteWs(), 
  remote = "remote", path = "remote")
{
  if (!dir.exists(path)) {
   stop('!dir.exists("path").')
  }
  if (is_sshfs_mount(path)) {
    return(message("The directory has been mount."))
  }
  if (length(list.files(path, include.dirs = TRUE))) {
    stop('length(list.files(path, all.files = TRUE, include.dirs = TRUE)).')
  }
  # umount remote
  cdRun(glue::glue("nohup sshfs {remote}:{ws}/{project} {path} >/dev/null 2>&1 &"))
  repeat {
    Sys.sleep(1)
    if (is_sshfs_mount(path)) {
      return(TRUE)
    }
  }
}

.upd_pkg_to_remote <- function(...) {
  .send_pkg_to_remote(..., upd = TRUE)
}

.send_pkg_to_remote <- function(from = "~/utils.tool",
  exclude = ".git", to = "remote", upd = FALSE, remoteUntar = TRUE, remote = "remote")
{
  if (!is_sshfs_mount(to)) {
    stop('!is_sshfs_mount(to).')
  }
  pkg <- basename(from)
  archive_package <- paste0(pkg, ".tar.gz")
  if (!upd && file.exists(file.path(to, archive_package))) {
    stop('!upd && file.exists(file.path(to, archive_package))')
  }
  cdRun(
    "git ls-files -z | xargs -0 tar -czf ", archive_package,
    path = from
  )
  file.copy(
    file.path(from, archive_package), to, TRUE
  )
  exdir <- file.path(to, pkg)
  dir.create(exdir, FALSE)
  if (remoteUntar) {
    ws <- getRemoteWs()
    pr <- guess_project()
    dir_project <- paste0(ws, "/", pr)
    cmd <- glue::glue("cd {dir_project} && tar -xzvf {archive_package} -C {pkg}")
    cdRun("ssh ", remote, " '", cmd, "'")
  } else {
    untar(
      normalizePath(file.path(to, archive_package)), exdir = exdir
    )
  }
}

is_sshfs_mount <- function(path = "remote") {
  type <- system(
    glue::glue("findmnt -n -o FSTYPE --target {path}"), intern = TRUE
  )
  type == "fuse.sshfs"
}

getRemoteWs <- function() {
  path <- getOption("remote_working_space")
  if (is.null(path)) {
    stop('is.null(path).')
  }
  path
}

guess_project <- function(path = getwd()) {
  res <- stringr::str_extract(basename(path), "[0-9]+_project[0-9]+")
  if (is.na(res)) {
    stop('is.na(res).')
  }
  res
}

new_script.hb <- function(theme, num = "guess", project = guess_project(), 
  ws = getRemoteWs(), 
  path = "remote", exlibrary = getOption("remote_R_library", ""))
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  if (missing(theme)) {
   stop('missing(theme).')
  }
  dir_project <- paste0(ws, "/", project)
  pattern <- glue::glue(
    "r\\.[0-9]{2}_{{{theme}}}\\.r", .open = "{{{", .close = "}}}"
  )
  existFiles <- list.files(path, pattern)
  if (length(existFiles)) {
    stop('length(existFiles).')
  }
  if (num == "guess") {
    num <- guess_number.hb(path)
  } else if (is.numeric(num)) {
    num <- sprintf("%02d", as.integer(num))
  }
  pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
  dir_output <- glue::glue("{dir_project}/{num}_{theme}")
  script <- readLines(file.path(.expath, "job_templ", "script_setup_huibang.R"))
  script <- glue::glue(
    paste0(script, collapse = "\n"), 
    ORIGINAL_DIR = dir_project, output = dir_output,
    LIBRARY = exlibrary,
    .open = ".{{{", .close = "}}}."
  )
  writeLines(script, pathScript)
  return(pathScript)
}

guess_number.hb <- function(path = "remote", p.pattern = "r\\.[0-9]{2}",
  n.pattern = "[0-9]{2}", type = c("files", "dirs"))
{
  type <- match.arg(type)
  if (type == "dirs") {
    alls <- list.dirs(path, recursive = FALSE)
    alls <- alls[ grpl(alls, p.pattern) ]
  } else {
    alls <- list.files(path, p.pattern)
  }
  num <- as.integer(stringr::str_extract(alls, n.pattern))
  num <- num[!is.na(num)]
  if (length(num)) {
    max <- max(num)
  } else {
    max <- 0L
  }
  sprintf("%02d", max + 1)
}

spsv <- function(object, name = "figure", prefix = "tmp_") {
  write_graphics(object, name = name, mkdir = ".")
}

take_positions <- function(plots, envir = .GlobalEnv) {
  calls <- substitute(plots)
  if (as_label(calls[[1]]) != "{") {
    stop('as_label(calls[[1]]) != "{"')
  }
  for (i in rev(seq_along(calls)[-1])) {
    dev.new()
    print(eval(parse(text = as_label(calls[[i]])), envir = envir))
  }
}

loadJob <- function(path, env = .GlobalEnv, name = "AUTO") {
  if (identical(name, "AUTO")) {
    name <- basename(path)
    name <- tools::file_path_sans_ext(name)
    name <- tools::file_path_sans_ext(name)
  }
  object <- readRDS(path)
  assign(name, object, envir = env)
  writeJobSlotsAutoCompletion(name, envir = env)
}

setup.huibang <- function() {
  options(
    wd_prefix = "/data/nas1/huanglichuang_OD/project/",
    db_prefix = "/data/nas1/huanglichuang_OD/project/",
    op_prefix = "/data/nas1/huanglichuang_OD/project/"
  )
  options("download.file.method" = "wget")
}

run_in_project_nohup <- function(script = "", remote = "remote") {
  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)
  cmd <- glue::glue("cd {dir_project} && nohup Rscript {script} > task.log 2>&1 &")
  cdRun("ssh ", remote, " '", cmd, "'")
}

