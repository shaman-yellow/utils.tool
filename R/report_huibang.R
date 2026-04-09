# ==========================================================================
# huibang
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

create_remote_project.hb <- function(project = guess_project(), ws = getRemoteWs(), 
  remote = "remote")
{
  cmd <- glue::glue("cd {ws} && mkdir {project}")
  cdRun("ssh ", remote, " '", cmd, "'")
}

setup.sshfs <- function(project = guess_project(), ws = getRemoteWs(), 
  remote = "remote", path = "remote", mirror = "scripts_mirror")
{
  if (!dir.exists(mirror)) {
    stop('!dir.exists(mirror).')
  }
  if (!dir.exists(path)) {
   stop('!dir.exists("path").')
  }
  file_sync <- file.path(.expath, "scripts", "sync.sh")
  cmd <- glue::glue("sh -c 'nohup bash {file_sync} {mirror} {remote}:{ws}/{project} > sync.log 2>&1 &'")
  system(cmd, wait = FALSE)
  # system2("bash", c("-c", cmd), wait = FALSE)
  if (is_sshfs_mount(path)) {
    return(message("The directory has been mount."))
  }
  if (length(list.files(path, include.dirs = TRUE))) {
    stop('length(list.files(path, all.files = TRUE, include.dirs = TRUE)).')
  }
  # umount remote
  # cdRun(glue::glue("nohup sshfs {remote}:{ws} ../{path} >/dev/null 2>&1 &"))
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

table_qc.hb <- function() {
  ftibble(file.path(.expath, "report_qc.csv"))
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
  pathPkgFrom <- file.path(from, archive_package)
  # if (file.exists(pathPkgFrom)) {
  #   message(glue::glue('file.exists(pathPkgFrom), remove ...'))
  #   file.remove(pathPkgFrom)
  # }
  cdRun(
    "git ls-files -c -o --exclude-standard -z | xargs -0 tar -czf ", archive_package,
    path = from
  )
  file.copy(
    pathPkgFrom, to, TRUE
  )
  exdir <- file.path(to, pkg)
  if (dir.exists(exdir)) {
    unlink(exdir, TRUE)
  }
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
  if (grpl(path, "_mirror$")) {
    return(file.exists(path))
  }
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
  res <- stringr::str_extract(basename(path), "[0-9]+_[a-zA-Z]+[0-9]+")
  if (is.na(res)) {
    stop('is.na(res).')
  }
  res
}

push_script.hb <- function(..., .project = guess_project(), 
  .ws = getRemoteWs(), .path = "scripts_mirror", .exlibrary = getOption("remote_R_library", ""))
{
  project <- .project
  ws <- .ws
  path <- .path
  exlibrary <- .exlibrary
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  dir_project <- file.path(ws, project)
  maxNum <- 1L
  allFiles <- list.files(path)
  vapply(list(...), FUN.VALUE = character(1),
    function(theme) {
      pattern <- glue::glue(
        "^r\\.[0-9]{2}_{{{theme}}}\\.r$", .open = "{{{", .close = "}}}"
      )
      existFiles <- grpf(allFiles, pattern)
      num <- sprintf("%02d", maxNum)
      maxNum <<- maxNum + 1L
      pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
      if (length(existFiles)) {
        if (length(existFiles) > 1) {
          rlang::abort(glue::glue("Theme of {theme} found multiple files: {bind(existFiles)}"))
        } else {
          numReal <- strx(existFiles, "[0-9]+")
          if (numReal != num && sureThat("File exists: {existFiles}, rename to r.{num}_{theme}.r?")) {
            file.rename(file.path(path, existFiles), pathScript)
            return(file.path(pathScript))
          }
          return(file.path(path, existFiles))
        }
      }
      dir_output <- glue::glue("{num}_{theme}")
      script <- readLines(file.path(.expath, "job_templ", "script_setup_huibang.R"))
      script <- glue::glue(
        paste0(script, collapse = "\n"), 
        ORIGINAL_DIR = dir_project, output = dir_output,
        LIBRARY = exlibrary,
        .open = ".{{{", .close = "}}}."
      )
      writeLines(script, pathScript)
      pathScript
    })
}

project_packaging.hb <- function(file_report, overwrite = FALSE, 
  overwrite_report = overwrite,
  path = "./remote", remote = "remote", report_share_to = "~/.var/app/com.tencent.WeChat/xwechat_files",
  wait = TRUE)
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  if (!file.exists(file_report)) {
    stop('!file.exists(file_report).')
  }
  ws <- getRemoteWs()
  pr <- guess_project()
  message(glue::glue("Workspace: {ws}\nProject: {pr}"))
  prefix <- strx(pr, "(?<=_)[a-zA-Z]+")
  message(glue::glue("Prefix: {prefix}"))
  num_project <- strx(pr, glue::glue("(?<={prefix})[0-9]+"))
  message(glue::glue("Project number: {prefix}{num_project}"))
  # remote is linux!
  dir_project <- paste0(ws, "/", pr)
  time <- format(Sys.Date(), "%Y%m%d")
  types <- c("scripts", "results", "report")
  names <- setNames(
    as.list(glue::glue("{prefix}_{num_project}_{types}_{time}")), types
  )
  message(glue::glue("Send report file..."))
  toDocx <- file.path(glue::glue("{path}"), glue::glue("{names$report}.docx"))
  if (!file.exists(toDocx) || overwrite_report) {
    file.copy(file_report, toDocx, TRUE)
    cdRun(glue::glue("soffice --headless --convert-to pdf --outdir {path} {toDocx}"))
  }
  ## scripts and files.
  all_scripts <- list.files(path, "^r\\.[0-9]+.*\\.r$")
  all_results <- gs(all_scripts, "^r\\.|\\.r$", "")
  fun_bind <- function(x) paste(x, collapse = " ")
  cmd_sed <- glue::glue("sed -i \"/^ORIGINAL_DIR\\|^.libPaths/d\" {names$scripts}/*")
  cmd_packaging_scripts <- glue::glue(
    "cd {dir_project} && mkdir {names$scripts} && cp -r {fun_bind(all_scripts)} -t {names$scripts} && {cmd_sed} && zip -r {names$scripts}.zip {names$scripts}"
  )
  cmd_packaging_results <- glue::glue(
    "cd {dir_project} && mkdir {names$results} && cp -r {fun_bind(all_results)} -t {names$results} && zip -r {names$results}.zip {names$results}"
  )
  if (!file.exists(file.path(path, glue::glue("{names$scripts}.zip"))) || overwrite) {
    message(glue::glue("Pacakging scripts: {names$scripts} ..."))
    cdRun("ssh ", remote, " '", cmd_packaging_scripts, "'", wait = wait)
    message(glue::glue("Pacakging results: {names$results} ..."))
    cdRun("ssh ", remote, " '", cmd_packaging_results, "'", wait = wait)
  }
  text_reply <- glue::glue("已上传分析报告:{dir_project}/{names$report}.docx和对应pdf\n\n代码压缩包:{dir_project}/{names$scripts}.zip\n\n结果文件压缩包:{dir_project}/{names$results}.zip")
  browseURL(normalizePath(path), "xdg-open")
  gett(text_reply)
  if (!is.null(report_share_to)) {
    file.copy(toDocx, report_share_to, TRUE)
    file.copy(
      paste0(tools::file_path_sans_ext(toDocx), ".pdf"), report_share_to, TRUE
    )
  }
}

pkgVersion_remote <- function(pkgs, path = "remote",
  exlibrary = getOption("remote_R_library", ""), remote = "remote")
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  exlibrary <- getOption("remote_R_library", "")
  strs <- glue::glue("'{pkgs}'")
  cmd <- glue::glue("invisible(lapply(c({bind(strs)}), function(x) writeLines(as.character(packageVersion(x)))))")
  lines <- c(exlibrary, cmd)
  dir.create(file.path(path, "tmp"), FALSE)
  file_script <- file.path(path, "tmp", "getPkgInfo.R")
  writeLines(lines, file_script)
  map_file <- file.path("tmp", "getPkgInfo.R")
  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)
  cmd <- glue::glue("cd {dir_project} && Rscript {map_file} ")
  res <- system(paste0("ssh ", remote, " '", cmd, "'"), intern = TRUE)
  if (length(res) != length(pkgs)) {
    stop('length(res) != length(pkgs).')
  }
  res
}

run_remote_output.hb <- function(run = FALSE, skip = NULL,
  files = list.files(path, "^r\\.[0-9]+.*\\.r$", full.names = TRUE),
  order_by_number = TRUE,
  path = "remote", cl = NULL)
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  tmpdir <- file.path(path, "tmp")
  dir.create(tmpdir, FALSE)
  if (order_by_number) {
    nums <- as.integer(strx(files, "[0-9]+"))
    files <- files[order(nums)]
  }
  allCodes <- pbapply::pblapply(seq_along(files), cl = cl,
    function(n) {
      file <- files[n]
      message(glue::glue("In script ({n}): {file}"))
      if (n %in% skip) {
        return()
      }
      lines <- readLines(file)
      field_analysis <- grp(lines, "^# FIELD: analysis")
      field_output <- grp(lines, "^# FIELD: output")
      codes <- lines[-(field_analysis:field_output)]
      fileName <- basename(file)
      writeLines(codes, file.path(tmpdir, fileName))
      if (run) {
        run_in_project(glue::glue("tmp/{fileName}"))
      }
      codes
    })
}

push_overture_as_output.hb <- function(pull = FALSE, push = FALSE,
  ovLoc = getOption("overture_codes_and_location"), override_remote = FALSE,
  dir_check = "./remote_script/push_check", path = "scripts_mirror",
  replace = "take_positions")
{
  if (is.null(ovLoc)) {
    stop('is.null(ovLoc), has not run `project_publish.complex`?')
  }
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  dir.create(dir_check, FALSE)
  dir_all <- vapply(ovLoc, function(x) x$dir, character(1))
  file_codes <- paste0("r.", dir_all, ".r")
  # setup, be careful, maybe multiple overture into the same file.
  path_codes <- file.path(path, unique(file_codes))
  if (pull) {
    file.copy(path_codes, dir_check, overwrite = TRUE)
  }
  path_codes_local <- file.path(dir_check, unique(file_codes))
  allCodes <- lapply(path_codes_local, 
    function(file) {
      codes <- readLines(file)
      posMark <- grp(codes, "^# FIELD: output")
      c(codes[1:(posMark + 1)], "", "setup_counting_in_directory(output)")
    })
  names(allCodes) <- basename(path_codes_local)
  # append the codes
  lapply(seq_along(ovLoc), 
    function(n) {
      file_code <- file_codes[[n]]
      mainCodes <- allCodes[[file_code]]
      code_output <- ovLoc[[n]]$codes
      code_output <- s(code_output, replace, "output_with_counting_number", fixed = TRUE)
      codes <- c(mainCodes, "", code_output)
      allCodes[[file_code]] <<- codes
    })
  pbapply::pblapply(seq_along(path_codes_local), 
    function(n) {
      if (push) {
        path_remote <- file.path(path, basename(path_codes_local[n]))
        writeLines(c(allCodes[[n]], "", ""), path_remote)
      } else {
        writeLines(c(allCodes[[n]], "", ""), path_codes_local[n])
      }
    })
}

setup_counting_in_directory <- function(dir, pattern = "^[0-9]+") {
  unlink(list.files(dir, pattern, full.names = TRUE), force = TRUE)
  options(autor_counting_start_dir = dir)
}

output_with_counting_number <- function(plots, envir = .GlobalEnv, 
  fun_wrap = "autor", extra_cmd = NULL)
{
  if (is.null(output <- getOption("autor_counting_start_dir"))) {
    stop('is.null(getOption("autor_counting_start_dir")).')
  }
  if (!dir.exists(output)) {
    stop('!dir.exists(output).')
  }
  calls <- substitute(plots)
  if (!is(calls, "{")) {
    stop('!is(calls, "{").')
  }
  rapp_find_job_name <- function(x) {
    if (is(x, "call") || is(x, "{")) {
      rapp_find_job_name(x[[2]])
    } else if (is(x, "name")) {
      rlang::expr_text(x)
    } else {
      stop("The finally of the 'substitute' is: ", class(x))
    }
  }
  num <- as.integer(guess_number.hb(output, p.pattern = "^[0-9]{2}"))
  fun_num <- function(n) {
    sprintf("%02d", n)
  }
  lapply(calls[-1], 
    function(call) {
      name <- rapp_find_job_name(call)
      job <- try(get(name, envir = .GlobalEnv), TRUE)
      if (inherits(job, "try-error")) {
        .try_loadJob(name, FALSE)
      }
      object <- eval(parse(text = rlang::expr_text(call)))
      outputName <- paste0(fun_num(num), "_", label(object))
      message(glue::glue("Save as: {outputName}"))
      fun_save <- select_savefun(object)
      fun_save(object, name = outputName, mkdir = output)
      num <<- num + 1L
    })
}

convert_pdf_in_project <- function(path = "remote", skip = NULL)
{
  dirs <- gs(
    list.files(path, "^r\\.[0-9]+.*\\.r$", full.names = TRUE), 
    "(?<=/)r\\.|\\.r$", "", perl = TRUE
  )
  order <- as.integer(strx(basename(dirs), "[0-9]+"))
  dirs <- dirs[ order(order) ]
  if (!is.null(skip)) {
    message(glue::glue("Skip: \n{bind(dirs[skip], co = '\n')}"))
    dirs <- dirs[-skip]
  }
  targets <- list.files(dirs, "\\.pdf$", full.names = TRUE)
  pbapply::pblapply(targets,
    function(file) {
      newfile <- paste0(tools::file_path_sans_ext(file), ".png")
      res <- try(pdf_convert(file, filenames = newfile, dpi = 300, pages = 1))
      if (inherits(res, "try-error")) {
        sink()
        message(glue::glue("Failed to convert file: {file}"))
      }
    })
}

pull_jobs_from_script.hb <- function(files, override = FALSE, 
  test = !override, project = guess_project(), 
  ws = getRemoteWs(), 
  path = "scripts_mirror", dir_save = "remote_script",
  pattern_object = "(?<=\\bclear\\()[a-zA-Z0-9_.]+",
  pattern_level = "(?<=\\bstep)[0-9]+(?=\\()")
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  if (missing(files)) {
    stop('missing(files).')
  }
  if (any(!file.exists(files))) {
    stop('any(!file.exists(files)).')
  }
  dir.create(dir_save, FALSE)
  fileNames <- basename(files)
  localFiles <- vapply(fileNames, FUN.VALUE = character(1),
    function(name) {
      local <- file.path(dir_save, name)
      if (!override && file.exists(local)) {
        stop('file.exists(local): ', local)
      }
      local
    })
  if (!test) {
    file.copy(files, dir_save, overwrite = override)
  }
  scripts <- lapply(localFiles, readLines)
  belongs <- lapply(scripts, 
    function(script) {
      unlist(stringr::str_extract_all(script, pattern_object))
    })
  belongs <- as_df.lst(belongs)
  belongs <- split(belongs$type, belongs$name)
  belongs <- sapply(names(belongs), simplify = FALSE,
    function(oname) {
      file <- unique(belongs[[oname]])
      if (length(file) > 1) {
        message(glue::glue("Detected `{oname}` from multiple file: {bind(file)}"))
        levels <- lapply(file.path(dir_save, file),
          function(file) {
            lines <- grpf(readLines(file), oname, fixed = TRUE)
            max(unlist(stringr::str_extract_all(lines, pattern_level)))
          })
        file <- file[which.max(levels)]
        message(glue::glue("{crayon::yellow(oname)} -> {file}"))
      }
      file
    })
  belongs <- lapply(belongs, 
    function(file) {
      dir <- gs(file, "^r\\.|\\.r$", "")
      list(script = file, dir = dir)
    })
  saveRDS(belongs, ".job_locate_in_script.rds")
  belongs
}

# new_script.hb <- function(theme, num = "guess", project = guess_project(), 
#   ws = getRemoteWs(), 
#   path = "remote", exlibrary = getOption("remote_R_library", ""))
# {
#   if (!is_sshfs_mount(path)) {
#     stop('!is_sshfs_mount(path).')
#   }
#   if (missing(theme)) {
#     stop('missing(theme).')
#   }
#   dir_project <- paste0(ws, "/", project)
#   pattern <- glue::glue(
#     "r\\.[0-9]{2}_{{{theme}}}\\.r", .open = "{{{", .close = "}}}"
#   )
#   existFiles <- list.files(path, pattern)
#   if (length(existFiles)) {
#     stop('length(existFiles).')
#   }
#   if (num == "guess") {
#     num <- guess_number.hb(path)
#   } else if (is.numeric(num)) {
#     num <- sprintf("%02d", as.integer(num))
#   }
#   pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
#   dir_output <- glue::glue("{dir_project}/{num}_{theme}")
#   script <- readLines(file.path(.expath, "job_templ", "script_setup_huibang.R"))
#   script <- glue::glue(
#     paste0(script, collapse = "\n"), 
#     ORIGINAL_DIR = dir_project, output = dir_output,
#     LIBRARY = exlibrary,
#     .open = ".{{{", .close = "}}}."
#   )
#   writeLines(script, pathScript)
#   return(pathScript)
# }

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

spsv <- function(object, name = NULL, prefix = "tmp") {
  if (is.null(name)) {
    name <- formal_name(rlang::expr_text(substitute(object)))
  }
  fun <- select_savefun(object)
  fun(object, name = name, mkdir = prefix)
}

smart_wrap_expr <- function(plots, size = 3, ..., envir = .GlobalEnv)
{
  calls <- substitute(plots)
  if (as_label(calls[[1]]) != "{") {
    stop('as_label(calls[[1]]) != "{"')
  }
  plots <- lapply(calls[-1],
    function(call) {
      eval(parse(text = as_label(call)), envir = envir)
    })
  smart_wrap(plots, size = size, ...)
}

expect_package <- function(pkg, version, prio_lib = getOption("prio_lib")) {
  if (!requireNamespace(pkg)) {
    stop('!requireNamespace(pkg)')
  }
  if (packageVersion(pkg) >= version) {
    message("Pacakge ", pkg, " as expected.")
    return()
  }
  if (packageVersion(pkg) < version) {
    message(glue::glue("Detach the loaded namespace, search in preferred lib path."))
    unloadNamespace(asNamespace(pkg))
    loadNamespace(pkg, lib.loc = prio_lib)
  }
  if (packageVersion(pkg, lib.loc = prio_lib) < version) {
    stop('packageVersion(pkg, lib.loc = prio_lib) < version.')
  } else {
    message(glue::glue("Successfully loaded the latter R package"))
  }
}

clear_feature <- function(x, name = "key.rds", dir = ".", 
  file = file.path(dir, name))
{
  if (!is(x, "feature")) {
    stop('!is(x, "feature").')
  }
  saveRDS(x, file)
}

load_feature <- function(name = "key.rds", dir = ".", file = file.path(dir, name)) {
  readRDS(file)
}

save_small.huibang <- function(name, cutoff = 50, dir = "rdata_smallObject")
{
  dir.create(dir, FALSE)
  file <- file.path(dir, glue::glue("{name}.rdata"))
  message(glue::glue("Save rdata: {file}"))
  save_small(cutoff = cutoff, file = file)
}


setup.huibang <- function() {
  conflicted::conflict_prefer("map", "utils.tool")
  options(
    tibble.print_max = 100,
    pillar.width = 100,
    pillar.max_columns = 15,
    prio_lib = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/lib/R/library/",
    digits = 4,
    warning.length = 5000,
    max.print = 500L,
    path_jobSave = "rds_jobSave",
    future.globals.maxSize = 5e10,
    auto_convert_plots = TRUE,
    wd_prefix = "/data/nas1/huanglichuang_OD/project/",
    db_prefix = "/data/nas1/huanglichuang_OD/project/",
    op_prefix = "/data/nas1/huanglichuang_OD/project/",
    file_batman_compounds_info = "/data/nas2/database/graphban/db/BATMAN_TCM/cids_result.csv",
    path_jobLoadFrom = list(remote = "./rds_jobSave/", local = "./rds_jobSave/lite/"),
    pg_local_recode = list(
      db_scenic = "/data/nas1/huanglichuang_OD/project/SCENIC",
      # db_scenic = "/data/nas2/database/SCENIC",
      pyscenic = "conda run -n pyscenic pyscenic",
      compass = "conda run -n mebocost compass",
      cellchat_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
      rdkit_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
      conda = "/data/nas2/software/miniconda3/bin/conda",
      scsaEnv = "scsa",
      mebocostEnv = "mebocost",
      path_mebocost = "/data/nas1/huanglichuang_OD/MEBOCOST",
      scsa = "conda run -n scsa python3 /data/nas1/huanglichuang_OD/SCSA/SCSA.py",
      scsa_db = "/data/nas1/huanglichuang_OD/SCSA/whole_v2.db"
    )
  )
  options("download.file.method" = "wget", "download.file.extra" = "--no-check-certificate")
}

run_in_project_nohup <- function(script, ...) {
  run_in_project(
    script, ..., wait = FALSE, ex1 = "nohup", ex2 = "> task.log 2>&1 &"
  )
}

run_in_project <- function(script = "", remote = "remote", 
  fun_map = NULL, wait = TRUE, ex1 = "", ex2 = "")
{
  if (!is.null(fun_map)) {
    script <- fun_map(script)
  }
  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)
  cmd <- glue::glue("cd {dir_project} && {ex1} Rscript {script} {ex2}")
  cdRun("ssh ", remote, " '", cmd, "'", wait = wait)
}

mark_text <- function(string, color, bold = TRUE, ...) {
  string <- gs(string, "&lt;", "<")
  string <- gs(string, "&gt;", ">")
  ftext <- officer::ftext(
    string, officer::fp_text_lite(color = color, bold = bold, ...)
  )
  paste0("`", officer::to_wml(ftext), "`{=openxml}")
}

mark <- list()

mark$red <- function(string) {
  mark_text(string, color = "#C00000")
}

mark$sig <- mark$red

mark$blue <- function(string) {
  mark_text(string, color = "#2E75B5")
}

mark$th <- mark$blue

mark$green <- function(string) {
  mark_text(string, color = "green")
}

name.hb <- list()

name.hb$check <- function() {
  date <- format(Sys.Date(), "%m%d")
  glue::glue("{s(guess_project(), '[0-9]+_', '')}_关键节点核对_{date}")
}

get_file_with_format_name <- function(file, name) {
  filename <- paste0(name, ".", tools::file_ext(file))
  file_new <- file.path(dirname(file), filename)
  file.copy(file, file_new, TRUE)
  if (nchar(Sys.which("wl-copy"))) {
    cdRun(glue::glue("wl-copy < {file_new}"))
  } else {
    stop('nchar(Sys.which("wl-copy")).')
  }
}

