# ==========================================================================
# workflow of sirius
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_sirius <- setClass("job_sirius", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://boecker-lab.github.io/docs.sirius.github.io/"),
    cite = "[@Sirius4ARapDuhrko2019; @SearchingMolecDuhrko2015; @SystematicClasDuhrko2021; @AssigningConfiHoffma2021]",
    method = "`SIRIUS` 5 (`SIRIUS`, `CSI:FingerID`, `CANOPUS`, `COSMIC`) used for compounds identification and prediction with MASS spectra."
    ))

job_sirius <- function(mgf, ion = c("mix", "pos", "neg"), pg = "~/operation/sirius/bin/sirius")
{
  ion <- match.arg(ion)
  x <- .job_sirius()
  x$mgf <- mgf
  x$ion <- ion
  x@pg <- pg
  x
}

setMethod("step0", signature = c(x = "job_sirius"),
  function(x){
    step_message("Prepare your data with function `job_sirius`.")
  })

setMethod("step1", signature = c(x = "job_sirius"),
  function(x, tryLog = T, account = "202011113511016@zcmu.edu.cn", password = "Qiu!23224856")
  {
    step_message("Login, if not.")
    if (tryLog) {
      rem_run(pg(x), " login -u ", account, " -p <<EOF",
        "\n", password, "\nEOF")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_sirius"),
  function(x, savepath, dir = paste0(x$ion, "_sirius"), workers = 5, maxmz = 900)
  {
    step_message("Run sirius 5 and output.")
    if (!is.null(savepath)) {
      rem_dir.create(savepath, F)
      x$output <- paste0(savepath, "/", dir)
    } else {
      x$output <- dir
    }
    x$workers <- workers
    if (rem_file.exists(x$output)) {
      message("Previous directory exsits, use it as input.")
      Sys.sleep(3)
      input <- x$output
    } else {
      input <- x$mgf
    }
    x$maxmz <- maxmz
    rem_run(
      pg(x), " --input ", input,
      " --output ", x$output,
      " --cores ", x$workers,
      " formula -p orbitrap fingerprint structure compound-classes write-summaries",
      if (!is.null(x$maxmz)) paste0(" --maxmz=", x$maxmz) else NULL,
      " --output ", x$output
    )
    return(x)
  })

setMethod("step3", signature = c(x = "job_sirius"),
  function(x, local)
  {
    step_message("If is runing remote, copy the output to local.")
    if (is.remote(x)) {
      x$map_local <- local
      target <- rem_normalizePath(x$output)
      cdRun("scp -r ", x$remote, ":", target, " ", x$map_local)
    }
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_sirius"),
  function(x, wd, postfix = NULL, run_after_cd = NULL, tmpdir = NULL, remote = "remote")
  {
    x$wd <- wd
    x$set_remote <- T
    x$remote <- "remote"
    x$map_local <- "sirius_local"
    if (!rem_file.exists(wd)) {
      rem_dir.create(wd)
    }
    x$local_mgf <- x$mgf
    x$mgf <- paste0(x$wd, "/", get_filename(x$local_mgf))
    if (rem_file.exists(x$mgf)) {
      isThat <- usethis::ui_yeah("Mgf file exsists, continue?")
      if (!isThat) {
        stop("...")
      }
    }
    cdRun("scp ", x$local_mgf, " ", remote, ":", x$mgf)
    x@step <- 0L
    return(x)
  })

setGeneric("asjob_sirius", 
  function(x, ...) standardGeneric("asjob_sirius"))

setMethod("asjob_sirius", signature = c(x = "job_xcms"),
  function(x, pg = "~/operation/sirius/bin/sirius")
  {
    params <- x@params
    x <- .job_sirius(params = params, pg = pg)
    return(x)
  })

setMethod("ping", signature = c(x = "job_sirius"),
  function(x, logout = F){
    if (logout) {
      rem_run(pg(x), " login --logout")
    } else {
      rem_run(pg(x), " login --show")
    }
  })


