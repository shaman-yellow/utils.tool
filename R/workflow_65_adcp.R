# ==========================================================================
# workflow of adcp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_adcp <- setClass("job_adcp", 
  contains = c("job_vina"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "adcp",
    info = c("https://ccsb.scripps.edu/adcp/tutorial-redocking/"),
    cite = "[@DockingFlexiblZhang2019; @AutositeAnAuRavind2016; @AutodockfrAdvRavind2015]",
    method = "`ADCP` and `ADFR` software used for auto molecular docking (docking peptides)",
    tag = "dock:pep",
    analysis = "ADCP, ADFR peptides docking"
    ))

job_adcp <- function(cids, hgnc_symbols)
{
  .job_adcp(object = namel(cids, hgnc_symbols))
}

setMethod("step0", signature = c(x = "job_adcp"),
  function(x){
    step_message("Prepare your data with function `job_adcp`.")
  })

setMethod("step1", signature = c(x = "job_adcp"),
  function(x, bdb = NULL, each_target = 1, pdbs = NULL)
  {
    callNextMethod(x, bdb, each_target, pdbs)
  })

setMethod("step2", signature = c(x = "job_adcp"),
  function(x, dir = "ligand", use = c("this", "nextMethod"), cl = 10)
  {
    use <- match.arg(use)
    if (use == "nextMethod") {
      return(callNextMethod(x, cl))
    }
    step_message("Download sdf files and convert as pdbqt for ligands.")
    sdfFile <- query_sdfs(unique(names(x$dock_layout)), curl_cl = cl)
    dir.create(dir, F)
    cdRun("obabel -isdf ", sdfFile, " -opdb -m -O ", dir, "/obabel.pdb ")
    pdbs <- lapply(list.files(dir, full.names = T),
      function(file) {
        name <- strx(readLines(file, n = 1), "[0-9]+")
        if (is.na(name)) {
          stop("Can not recognize the compound id of: ", file)
        }
        file.rename(file, nfile <- paste0(dir, "/", name, ".pdb"))
        nl(name, nfile, F)
      })
    pdbqts <- lapply(unlist(pdbs),
      function(file) {
        file <- get_filename(file)
        cdRun("reduce ", file, " > ", fileH <- gs(file, "\\.pdb$", "_H.pdb"),
          path = dir)
        if (file.exists(paste0(dir, "/", fileH))) {
          cdRun("prepare_ligand -l ", fileH, " -o ", nfile <- gs(fileH, "\\.pdb", ".pdbqt"),
            path = dir)
          nfile <- paste0(dir, "/", nfile)
          if (file.exists(nfile)) {
            nfile
          } else {
            stop("The file not generated succeed: ", nfile)
          }
        } else {
          stop("The file not generated succeed: ", dir, "/", fileH)
        }
      })
    x$res.ligand <- pdbqts
    message("Got: ", length(x$res.ligand))
    alls <- as.character(object(x)$cid)
    message("Not got:", paste0(alls[!alls %in% names(x$res.ligand)], collapse = ", "))
    return(x)
  })

setMethod("step3", signature = c(x = "job_adcp"),
  function(x, cl = 10, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS",
    extra_pdb.files = NULL, extra_layouts = NULL, extra_symbols = NULL)
  {
    x <- callNextMethod()
    return(x)
  })

setMethod("step4", signature = c(x = "job_adcp"),
  function(x, time = 3600, savedir = "adcp_space", log = "/tmp/res.log", save.object = "ad3.rds")
  {
    step_message("Run ADFR and ADCP ...")
    runs <- tibble::tibble(
      Ligand = rep(names(x$dock_layout), lengths(x$dock_layout)),
      Receptor = tolower(unlist(x$dock_layout, use.names = F))
    )
    x$show_layout <- runs
    x$runs <- runs <- apply(runs, 1, unname, simplify = F)
    x$savedir <- savedir
    n <- 0
    if (file.exists(log)) {
      file.remove("/tmp/res.log")
    }
    saveRDS(x, save.object)
    pbapply::pblapply(x$runs,
      function(v) {
        lig <- x$res.ligand[[ v[1] ]]
        recep <- x$res.receptor[[ tolower(v[2]) ]]
        n <<- n + 1
        adcp(lig, recep, time, dir = savedir)
      }
    )
    return(x)
  })

adcp <- function(lig, recep, timeLimit = 3600, dir = "adcp_space", stout = "/tmp/res.log")
{
  if (!file.exists(dir)) {
    dir.create(dir, F)
  } 
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
  dir.create(wd, F)
  file.copy(c(recep, lig), wd)
  .message_info("Generate the target file containing the affinity maps: ", subdir)
  files <- get_filename(c(lig, recep))
  cdRun("agfr -r ", files[2],
    " -l ", files[1],
    " -o ", subdir,
    path = wd
  )
  if (file.exists(paste0(wd, "/", subdir, ".trg"))) {
    cat("\n$$$$\n", date(), "\n", subdir, "\n\n", file = stout, append = T)
    try(cdRun("timeout ", timeLimit, " adcp",
        " -t ", subdir, ".trg",
        " -s npisdvd",
        " -N 20",
        " -n 1000000",
        " -o redocking ",
        " -ref ", files[1],
        " >> ", stout, path = wd), T)
  } else {
    warning("File ", wd, "/", subdir, ".trg not exists.")
  }
}
