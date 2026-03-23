# ==========================================================================
# huibang
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setup.sshfs <- function(project = guess_project(), ws = getRemoteWs(), 
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
  cdRun(glue::glue("nohup sshfs {remote}:{ws} ../{path} >/dev/null 2>&1 &"))
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

spsv <- function(object, name = NULL, prefix = "tmp_") {
  if (is.null(name)) {
    name <- formal_name(rlang::expr_text(substitute(object)))
  }
  fun <- select_savefun(object)
  fun(object, name = name, mkdir = ".")
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
    cellchat_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
    rdkit_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
    path_jobLoadFrom = list(remote = "./rds_jobSave/"),
    pg_local_recode = list(
      # fusion = "~/fusion_twas",
      # ldscPython = "{conda}/bin/conda run -n ldsc python",
      # ldsc = "~/ldsc",
      # annovar = "~/disk_sda1/annovar",
      # vep = "~/ensembl-vep/vep",
      # vep_cache = "~/disk_sda1/.vep",
      # vina = "vina",
      # python = "{conda}/bin/python3",
      # conda = "{conda}/bin/conda",
      # conda_env = "{conda}/envs",
      # qiime = "{conda}/bin/conda run -n qiime2 qiime",
      # musitePython = "{conda}/bin/conda run -n musite python3",
      # musitePTM = "~/MusiteDeep_web/MusiteDeep/predict_multi_batch.py",
      # musitePTM2S = "~/MusiteDeep_web/PTM2S/ptm2Structure.py",
      # hobEnv = "hobpre",
      # hobPython = "{conda}/bin/conda run -n hobpre python",
      # hobPredict = "~/HOB/HOB_predict.py",
      # hobModel = "~/HOB/model",
      # hobExtra = "~/HOB/pca_hob.m",
      # dl = normalizePath("~/D-GCAN/DGCAN"),
      # dl_dataset = normalizePath("~/D-GCAN/dataset"),
      # dl_model = normalizePath("~/D-GCAN/DGCAN/model"),
      # scfeaPython = "{conda}/bin/conda run -n scFEA python",
      # scfea = "~/scFEA/src/scFEA.py",
      # scfea_db = "~/scFEA/data",
      # musiteModel = normalizePath("~/MusiteDeep_web/MusiteDeep/models"),
      # mk_prepare_ligand.py = "mk_prepare_ligand.py",
      # prepare_receptor = "prepare_receptor",
      # prepare_gpf.py = "prepare_gpf.py",
      # autogrid4 = "autogrid4",
      # pymol = "/usr/bin/python3 -m pymol",
      # sirius = .prefix("sirius/bin/sirius", "op"),
      # obgen = "obgen",
      conda = "/data/nas2/software/miniconda3/bin/conda",
      scsaEnv = "scsa",
      scsa = "conda run -n scsa python3 /data/nas1/huanglichuang_OD/SCSA/SCSA.py",
      scsa_db = "/data/nas1/huanglichuang_OD/SCSA/whole_v2.db"
    )
  )
  options("download.file.method" = "wget", "download.file.extra" = "--no-check-certificate")
}

run_in_project_nohup <- function(script = "", remote = "remote", fun_map = basename)
{
  script <- fun_map(script)
  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)
  cmd <- glue::glue("cd {dir_project} && nohup Rscript {script} > task.log 2>&1 &")
  cdRun("ssh ", remote, " '", cmd, "'")
}

