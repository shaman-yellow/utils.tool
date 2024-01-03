# ==========================================================================
# workflow of swiss
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_swiss <- setClass("job_swiss", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://www.swisstargetprediction.ch/index.php"),
    cite = "[@SwisstargetpredDaina2019]",
    method = "Web tool of `SwissTargetPrediction` used for drug-targets prediction"
    ))

job_swiss <- function(smiles)
{
  .job_swiss(object = smiles)
}

setMethod("step0", signature = c(x = "job_swiss"),
  function(x){
    step_message("Prepare your data with function `job_swiss`.")
  })

setMethod("step1", signature = c(x = "job_swiss"),
  function(x, db_file = "../swissTargetPrediction/targets.rds", tempdir = "download", sleep = 5)
  {
    step_message("Touch the online tools.")
    x$tempdir <- tempdir
    db <- new_db(db_file, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      link <- start_drive(download.dir = x$tempdir)
      Sys.sleep(3)
      link$open()
      lapply(db@query,
        function(query) {
          link$navigate("http://www.swisstargetprediction.ch/index.php")
          ele <- link$findElement("xpath", "//form//div//input[@id='smilesBox']")
          ele$sendKeysToElement(list(query))
          Sys.sleep(3)
          Sys.sleep(sleep)
          ele <- link$findElement("xpath", "//form//div//p//input[@id='submitButton']")
          ele$clickElement()
          Sys.sleep(3)
          Sys.sleep(sleep)
          ele <- F
          n <- 0L
          while ((is.logical(ele) | inherits(ele, "try-error")) & n < 20L) {
            Sys.sleep(1)
            n <- n + 1L
            ele <- try(link$findElement("xpath", "//div//button[@class='dt-button buttons-csv buttons-html5']"), T)
          }
          if (!inherits(ele, "try-error")) {
            ele$clickElement()
          }
          Sys.sleep(1)
          Sys.sleep(sleep)
        })
      link$close()
      end_drive()
      ids <- paste0("query", 1:length(db@query))
      files <- moveToDir(ids, "SwissTargetPrediction.*csv", from = x$tempdir, to = x$tempdir,
        suffix = ".csv")
      data <- ftibble(files)
      names(data) <- db@query
      data <- frbind(data, idcol = T, fill = T)
      db <- upd(db, data)
    }
    targets <- dplyr::filter(as_tibble(db@db), .id %in% object(x))
    targets <- dplyr::rename(targets, smiles = .id)
    x@tables[[ 1 ]] <- namel(targets)
    return(x)
  })
