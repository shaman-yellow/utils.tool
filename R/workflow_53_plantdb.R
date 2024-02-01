# ==========================================================================
# workflow of plantdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_plantdb <- setClass("job_plantdb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://plantaedb.com/"),
    cite = "",
    method = "The Database `PlantaeDB` <https://plantaedb.com/> used for collating data of herbal ingredients"
    ))

job_plantdb <- function(query)
{
  .job_plantdb(object = query)
}

setMethod("step0", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Prepare your data with function `job_plantdb`.")
  })

setMethod("step1", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Query plants.")
    res <- sapply(object(x), simplify = F,
      function(query) {
        url <- paste0("https://plantaedb.com/search?src=", gs(query, " ", "%20"), "&tab=taxa")
        res <- RCurl::getURL(url)
        res <- try(get_table.html(res, elFun = tryGetLink.plantaedb), T)
        if (!inherits(res, "try-error")) {
          trySepCols.tcmsp(as_tibble(res[[1]]))
        } else {
          data.frame()
        }
      })
    res <- frbind(res, fill = T, idcol = T)
    x$herbs_info <- dplyr::relocate(res, .id, Result, Authority, Rank)
    return(x)
  })

setMethod("slice", signature = c(x = "job_plantdb"),
  function(x, ...){
    x@params$herbs_info <- dplyr::slice(x@params$herbs_info, ...)
    return(x)
  })

setMethod("step2", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Collating data from URLs.")
    t.data <- lapply(x$herbs_info$Result.link,
      function(x) {
        x <- paste0("https://plantaedb.com", x)
        get_plantaedb_data(x)$compounds
      })
    names(t.data) <- x$herbs_info$Result
    t.data <- frbind(t.data, idcol = T)
    x@tables[[ 2 ]] <- namel(t.data)
    return(x)
  })

setMethod("step3", signature = c(x = "job_plantdb"),
  function(x, hob_filter = F, ...)
  {
    step_message("Use PubChemR to obtain compounds information.")
    data <- x@tables$step2$t.data
    query <- nl(data$Name, data$`PubChem ID`, F)
    if (!identical(duplicated(query), duplicated(names(query)))) {
      warning("identical(duplicated(query), duplicated(names(query))) == F")
    }
    query <- query[ !duplicated(query) ]
    pr <- job_pubchemr(query)
    x$pr <- suppressMessages(step1(pr))
    tables <- list()
    plots <- list()
    if (hob_filter) {
      ho <- asjob_hob(x$pr)
      x$ho <- suppressMessages(step1(ho, ...))
      res <- x$ho@tables$step1$t.hob
      x$pr <- filter(x$pr, as.logical(res$prediction))
      .add_internal_job(ho)
    }
    .add_internal_job(pr)
    if (length(tables))
      x@tables[[ 3 ]] <- tables
    if (length(plots))
      x@plots[[ 3 ]] <- plots
    return(x)
  })

setMethod("step4", signature = c(x = "job_plantdb"),
  function(x, db_file = .prefix("superPred/targets.rds", "db"), tempdir = "download", port = 4444)
  {
    step_message("Use Super-Pred to get compounds targets.")
    sp <- asjob_superpred(x$pr)
    x$sp <- suppressMessages(step1(sp, db_file = db_file, tempdir = tempdir, port = port))
    .add_internal_job(sp)
    return(x)
  })

setMethod("step5", signature = c(x = "job_plantdb"),
  function(x, vis = F){
    step_message("Network pharmacology preparation.")
    metadata <- dplyr::distinct(x@tables[[2]]$t.data, herb = .id, cid = `PubChem ID`)
    x$hb <- do_herb(x$pr, x$sp, run_step3 = vis, metadata = metadata)
    return(x)
  })

setMethod("map", signature = c(x = "job_herb", ref = "job_plantdb"),
  function(x, ref){
    if (ref@step < 5L) {
      stop("ref@step < 5L")
    }
    fun_add <- function(x, y) {
      dplyr::bind_rows(x, y)
    }
    x@tables$step1$herbs_compounds %<>% fun_add(ref$hb@tables$step1$herbs_compounds)
    x@tables$step2$compounds_targets %<>% fun_add(ref$hb@tables$step2$compounds_targets)
    x@params$herbs_info %<>% fun_add(ref$hb@params$herbs_info)
    x@object$herb %<>% fun_add(ref$hb@object$herb)
    x@step <- 2L
    message("Set `step` to 2L")
    return(x)
  })

# setMethod("step2", signature = c(x = "job_plantdb"),
#   function(x){
#     step_message("Obtaining targets.")
#     # https://go.drugbank.com/releases/latest
#     return(x)
#   })

get_plantaedb_data <- function(url) {
  html <- RCurl::getURL(url)
  tables <- get_table.html(html)
  tables <- lapply(tables,
    function(x) {
      if (any(colnames(x) %in% c("Internal ID", "PubChem ID"))) {
        as_tibble(x)
      } else NULL
    })
  tables <- lst_clear0(tables)
  if (length(tables) == 2)
    names(tables) <- c("info", "compounds")
  else
    names(tables) <- "info"
  fun_format <- function(x) {
    seps <- grp(x$Name, "^>")
    group <- 0L
    index <- 0L
    seps <- vapply(1:nrow(x), FUN.VALUE = double(1),
      function(n) {
        if (is.na(seps[ index + 1L ]) || n > seps[ index + 1L ]) {
          group <<- group + 1L
          index <<- index + 1L
        } else if (n == seps[ index + 1L ]) {
          group <- 0L
        }
        group
      })
    x <- split(x, factor(seps, levels = unique(seps)))
    names <- x[[ "0" ]]$Name
    x <- x[-1]
    names(x) <- names
    x <- dplyr::rename(frbind(x, idcol = T), classes = .id)
    x <- dplyr::mutate(x, `Canonical SMILES` = gs(`Canonical SMILES`, "^Click to see\n\\s*", ""))
    x
  }
  if (nrow(tables$compounds))
    tables$compounds <- fun_format(tables$compounds)
  tables
}

tryGetLink.plantaedb <- function(x, sep = " ### ", ...) {
  val <- XML::xmlValue(x)
  children <- XML::xmlChildren(x)
  if (any(names(children) == "a")) {
    attrs <- XML::xmlAttrs(children$a)
    if (any(names(attrs) == "href")) {
      paste0(val, sep, attrs[[ "href" ]])
    } else {
      val
    }
  } else {
    val
  }
}

get_bindingdb_data <- function(url = "https://www.bindingdb.org/bind/downloads/BindingDB_All_202401_tsv.zip",
  save = paste0(.prefix(name = "db"), get_filename(url)),
  file_unzip = gs(gs(save, ".zip$", ""), "_tsv$", ".tsv"))
{
  if (!file.exists(file_unzip)) {
    if (!file.exists(save)) {
      cdRun("wget ", url, " -O ", save)
    }
    utils::unzip(save, exdir = get_path(save))
  }
}

get_drugbank_data <- function(url = "https://go.drugbank.com/releases/5-1-11/downloads/all-full-database",
  user = c("202011113511016@zcmu.edu.cn", "qiu23224856"),
  save = .prefix("drugbank/drugbank_all_full_database.xml.zip", "db"),
  file_unzip = paste0(get_path(save), "/full database.xml"))
{
  if (!file.exists(file_unzip)) {
    if (!file.exists(save)) {
      dir.create(get_path(save), F)
      cdRun("wget ", " --http-user=", user[1], " --http-passwd=", user[2], " ", url, " -O ", save)
    }
    utils::unzip(save, exdir = get_path(save))
  }
}

