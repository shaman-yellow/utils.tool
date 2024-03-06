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
  if (any(nchar(smiles) > 200)) {
    stop("any(nchar(smiles) > 200)")
  }
  .job_swiss(object = smiles)
}

setGeneric("asjob_swiss", 
  function(x, ...) standardGeneric("asjob_swiss"))

setMethod("asjob_swiss", signature = c(x = "job_tcmsp"),
  function(x){
    data <- x@tables$step2$ingredients
    .check_columns(data, c("Mol ID", "smiles"))
    job_swiss(nl(data$`Mol ID`, data$smiles, F))
  })

setMethod("asjob_swiss", signature = c(x = "job_pubchemr"),
  function(x){
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    job_swiss(unlist(x@params$smiles))
  })

setMethod("step0", signature = c(x = "job_swiss"),
  function(x){
    step_message("Prepare your data with function `job_swiss`.")
  })

setMethod("step1", signature = c(x = "job_swiss"),
  function(x, db_file = .prefix("swissTargetPrediction/targets.rds", "db"), tempdir = "download", sleep = 5, port = 4444)
  {
    step_message("Touch the online tools.")
    x$tempdir <- tempdir
    db <- new_db(db_file, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      link <- start_drive(download.dir = x$tempdir, port = port)
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
          Sys.sleep(20)
        })
      link$close()
      end_drive()
      ids <- paste0("query", 1:length(db@query))
      files <- collateFiles(ids, "SwissTargetPrediction.*csv", from = x$tempdir, to = x$tempdir,
        suffix = ".csv")
      data <- ftibble(files)
      names(data) <- db@query
      data <- frbind(data, idcol = T, fill = T)
      db <- upd(db, data)
    }
    targets <- dplyr::filter(as_tibble(db@db), .id %in% object(x))
    targets <- dplyr::rename(targets, smiles = .id, symbols = `Common name`)
    targets <- split_lapply_rbind(targets, 1:nrow(targets), args = list(fill = T),
      function(x) {
        if (grpl(x$symbols, " ")) {
          symbols <- strsplit(x$symbols, " ")[[1]]
          x$symbols <- NULL
          data.frame(x, symbols = symbols, check.names = F)
        } else x
      })
    x@tables[[ 1 ]] <- namel(targets)
    return(x)
  })

setMethod("map", signature = c(x = "job_tcmsp", ref = "job_swiss"),
  function(x, ref){
    refTargets <- ref@tables$step1$targets
    meta <- x@tables$step2$ingredients
    meta <- dplyr::select(meta, smiles, `Mol ID`, `Molecule Name`)
    refTargets <- tbmerge(refTargets, meta, by = "smiles", all.x = T)
    refTargets <- dplyr::select(refTargets, `Mol ID`, `Molecule Name`,
      `Target name` = Target, symbols, probability = `Probability*`)
    refTargets <- dplyr::filter(refTargets, `Mol ID` %in% x@tables$step2$ingredients$`Mol ID`)
    x@tables$step2$compounds_targets <- refTargets
    data <- split_lapply_rbind(refTargets, ~ `Molecule Name`,
      function(x) {
        x <- dplyr::filter(x, probability > 0)
        x <- tibble::add_row(x, symbols = "Others.",
          probability = 0, `Molecule Name` = x$`Molecule Name`[[1]])
        x
      })
    p.targets <- ggplot(data) +
      geom_col(aes(x = reorder(symbols, probability), y = probability, fill = probability)) +
      labs(y = "Targets", x = "Probability") +
      coord_flip() +
      guides(fill = "none") +
      facet_wrap(~ `Molecule Name`, scales = "free_y")
    p.targets <- wrap(p.targets)
    x$p.swissTargets <- .set_lab(p.targets, sig(x), "SwissTargetPrediction-results")
    return(x)
  })
