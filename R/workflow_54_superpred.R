# ==========================================================================
# workflow of superpred
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_superpred <- setClass("job_superpred", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://prediction.charite.de/subpages/target_prediction.php"),
    cite = "[@SuperpredUpdaNickel2014]",
    method = "Web tool of `Super-PRED` used for drug-targets prediction"
    ))

job_superpred <- function(smiles)
{
  .job_superpred(object = smiles)
}

setMethod("step0", signature = c(x = "job_superpred"),
  function(x){
    step_message("Prepare your data with function `job_superpred`.")
  })

setMethod("step1", signature = c(x = "job_superpred"),
  function(x, db_file = "../superPred/targets.rds", tempdir = "download", port = 4444)
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
          link$navigate("https://prediction.charite.de/subpages/target_prediction.php")
          ele <- link$findElement("xpath", "//form//div//input[@id='smiles_string']")
          ele$sendKeysToElement(list(query))
          Sys.sleep(3)
          ele <- link$findElement("xpath", "//form//div//input[@id='smiles_string']/../div[@class='input-group-append']/button")
          ele$sendKeysToElement(list("Search", key = "enter"))
          Sys.sleep(3)
          ele <- link$findElement("xpath", "//table//td/button[@type='submit']")
          ele$sendKeysToElement(list("Start Calculation", key = "enter"))
          ele <- F
          n <- 0L
          while ((is.logical(ele) | inherits(ele, "try-error")) & n < 20L) {
            Sys.sleep(1)
            n <- n + 1L
            ele <- try(link$findElement("xpath",
                "//div[@id='targets_wrapper']//div/button[@class='dt-button buttons-csv buttons-html5']"), T)
          }
          if (!inherits(ele, "try-error")) {
            ele$clickElement()
          }
          Sys.sleep(5)
        })
      link$close()
      end_drive()
      ids <- paste0("query", 1:length(db@query))
      files <- moveToDir(ids, "Targets.*csv", from = x$tempdir, to = x$tempdir,
        suffix = ".csv")
      data <- ftibble(files)
      names(data) <- db@query
      data <- frbind(data, idcol = T, fill = T)
      db <- upd(db, data)
    }
    targets <- dplyr::filter(as_tibble(db@db), .id %in% object(x))
    symbols <- e(UniProt.ws::mapUniProt(
        from = "UniProtKB_AC-ID", to = "Gene_Name",
        columns = c("accession", "id"),
        query = list(ids = unique(targets[[ 'UniProt ID' ]]))
        ))
    targets <- map(targets, "UniProt ID", symbols, "From", "To", col = "symbols")
    x@tables[[ 1 ]] <- namel(targets)
    return(x)
  })

setGeneric("asjob_superpred", 
  function(x, ...) standardGeneric("asjob_superpred"))

setMethod("asjob_superpred", signature = c(x = "job_pubchemr"),
  function(x){
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    job_superpred(unlist(x@params$smiles))
  })

setGeneric("do_herb", 
  function(x, ref, ...) standardGeneric("do_herb"))

setMethod("do_herb", signature = c(x = "job_pubchemr", ref = "job_superpred"),
  function(x, ref, disease = NULL, disease.score = 5, HLs = NULL, names = NULL)
  {
    if (is.null(names)) {
      names <- object(x)
      if (is.null(names)) {
        stop("is.null(names)")
      }
    } else {
      object(x) <- names
    }
    plots <- list()
    if (T) {
      targets <- dplyr::rename(ref@tables$step1$targets, smiles = .id)
      if (is.null(names(object(x)))) {
        stop("is.null(names(object(x)))")
      }
      dic <- nl(unlist(x$smiles, use.names = F), names(x$smiles))
      dic2 <- nl(unname(object(x)), names(object(x)))
      targets <- dplyr::mutate(targets, pubchem_id = dplyr::recode(smiles, !!!dic),
        name = dplyr::recode(pubchem_id, !!!dic2)
      )
      data <- dplyr::select(targets, name, symbols, probability = Probability, `Model accuracy`)
      data <- dplyr::mutate(data, probability = as.double(strx(probability, "^[0-9]+")),
        `Model accuracy` = as.double(strx(`Model accuracy`, "^[0-9]+")))
      data <- split_lapply_rbind(data, ~ name,
        function(x) {
          x <- head(x, n = 15)
          tibble::add_row(x, symbols = "Others.", probability = 0, name = x$name[[1]],
            `Model accuracy` = 0)
        })
      p.targets <- ggplot(data) +
        geom_col(aes(x = reorder(symbols, probability), y = probability, fill = `Model accuracy`)) +
        labs(y = "Targets", x = "Probability (%)", fill = "Model accuracy (%)") +
        coord_flip() +
        facet_wrap(~ name, scales = "free_y")
      p.targets <- wrap(p.targets)
      p.targets <- .set_lab(p.targets, sig(x), "SuperPred-results")
      plots <- c(plots, namel(p.targets))
    }
    data <- tibble::tibble(
      Herb_pinyin_name = "PseudoHerb",
      herb_id = "PseudoHerb_id",
      Herb_cn_name = "PseudoHerb_cn",
      Ingredient.id = unname(object(x)),
      Ingredient.name = names(object(x))
    )
    hb <- .job_herb(step = 2L)
    hb@tables$step1$herbs_compounds <- dplyr::select(
      data, herb_id, Ingredient.id, Ingredient.name
    )
    hb@tables$step2$compounds_targets <- dplyr::select(targets,
      Ingredient_id = pubchem_id, Target.name = symbols,
      Target.protein = `Target Name`
    )
    hb@params$herbs_info <- dplyr::select(
      dplyr::mutate(data, Herb_ = Herb_pinyin_name),
      Herb_, Herb_pinyin_name, Herb_cn_name
    )
    hb@object$herb <- hb@params$herbs_info
    hb <- suppressMessages(step3(hb, disease = disease, HLs = HLs))
    hb@plots[[ 3 ]] <- c(plots, hb@plots$step3)
    hb@tables[[ 3 ]] <- namel(disease_targets_annotation = hb@tables$step3$disease_targets_annotation)
    return(hb)
  })
