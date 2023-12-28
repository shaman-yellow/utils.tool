# ==========================================================================
# workflow of tcmsp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tcmsp <- setClass("job_tcmsp", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://tcmsp-e.com/tcmsp.php"),
    cite = "[@TcmspADatabaRuJi2014]",
    method = "Website `TCMSP` <https://tcmsp-e.com/tcmsp.php> used for data source"
    ))

job_tcmsp <- function(herbs, db = get_tcmsp_data())
{
  db$herbs <- dplyr::mutate(db$herbs, Herb_cn_name = gs(`Chinese name`, "^([^ ]+).*", "\\1"),
    Herb_pinyin_name = gs(`Chinese name`, ".*\\(([a-zA-Z]+)\\)", "\\1")
  )
  x <- .job_tcmsp(object = db)
  x@params$herbs <- herbs
  herbs_info <- filter(object(x)$herb, Herb_cn_name %in% !!params(x)$herbs)
  x@params$herbs_info <- herbs_info
  print(herbs_info)
  message("Got the herbs:\n\n\t", paste0(herbs_info$Herb_cn_name, collapse = ", "),
    "\n\n", "\tTotal: ", colSum(herbs_info$Herb_cn_name))
  x
}

setMethod("step0", signature = c(x = "job_tcmsp"),
  function(x){
    step_message("Prepare your data with function `job_tcmsp`.")
  })

setMethod("step1", signature = c(x = "job_tcmsp"),
  function(x, savedir = "../tcmsp/largedata"){
    step_message("Dowload detail information of herbs.")
    if (!dir.exists(savedir)) {
      dir.create(savedir)
    }
    items <- c(
      "ingredients" = "Ingredients",
      "targets" = "Related Targets",
      "disease" = "Related Diseases"
    )
    isDriveStarted <- F
    link <- F
    db_file <- as.list(.list_allfiles(savedir, "^herb_"))
    n <- 0L
    lst <- sapply(x$herbs_info$Herb_pinyin_name, simplify = F,
      function(name) {
        n <<- n + 1L
        fileid <- paste0("herb_", name, "_", names(items))
        files <- paste0(savedir, "/", fileid, ".tsv")
        fun <- function(n) is.null(db_file[[ fileid[n] ]])
        if (fun(1) | fun(2) | fun(3)) {
          link <<- start_drive(browser = "firefox")
          Sys.sleep(3)
          isDriveStarted <<- T
          link$open()
          url <- x$herbs_info$`Latin name.link`[n]
          status <- try(link$navigate(url), T)
          if (inherits(status, "try-error")) {
            res <- usethis::ui_yeah("Do you want to continue?")
            if (!res) {
              stop("Stop.")
            }
          }
          sets <- lapply(1:length(items),
            function(i) {
              path <- paste0("//div//ul//li//a[text()='", items[i], "']")
              ele <- try(link$findElement(using = "xpath", value = path), T)
              if (inherits(ele, "try-error")) {
                stop("In URL: ", url, "\n\tCan not find the element of path: ", path)
              }
              ele$clickElement()
              Sys.sleep(.1)
              ###########
              ###########
              if (items[i] == "Related Targets") {
                obj <- .get_folded_tables.tcmsp(link, "[@id='grid2']", 3:4)
              } else if (items[i] == "Related Diseases") {
                obj <- .get_folded_tables.tcmsp(link, "[@id='grid3']", 7:8)
              } else {
                obj <- .get_folded_tables.tcmsp(link)
              }
              ###########
              ###########
              obj <- .try_format_folded_tables.tcmsp(obj)
              write_tsv(obj, files[i])
              obj
            })
        } else {
          sets <- ftibble(files)
        }
        names(sets) <- names(items)
        sets
      })
    if (isDriveStarted) {
      link$close()
      end_drive()
    }
    tables <- sapply(names(items), simplify = F,
      function(name) {
        data <- frbind(lapply(lst, function(x) x[[name]]), idcol = T, fill = T)
        dplyr::relocate(data, Herb_pinyin_name = .id)
      })
    x@tables[[ 1 ]] <- tables
    x$savedir <- savedir
    return(x)
  })

setMethod("step2", signature = c(x = "job_tcmsp"),
  function(x, mode = c("tcmsp", "pubchem")){
    step_message("Dowload detail information of compounds")
    data <- x@tables$step1$ingredients
    mode <- match.arg(mode)
    if (mode == "tcmsp") {
      compounds_info <- .get_compounds_info.tcmsp(x, data)
    } else {
      stop("...")
    }
    compounds_info <- dplyr::select(compounds_info, `Mol ID` = .id, InChIKey)
    data <- tbmerge(data, compounds_info, by = "Mol ID")
    targets <- dplyr::select(x@tables$step1$targets, `Mol ID`, `Target name`)
    targets <- dplyr::distinct(targets)
    compounds_targets <- tbmerge(data, targets, by = "Mol ID", all.x = T)
    x@tables[[ 2 ]] <- namel(compounds_info, ingredients = data, compounds_targets)
    return(x)
  })

setMethod("map", signature = c(x = "job_tcmsp", ref = "job_classyfire"),
  function(x, ref) {
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    if (ref@step < 2L) {
      stop("ref@step < 2L")
    }
    if (is.null(ref$match)) {
      stop("is.null(ref$match)")
    }
    if (is.null(x$step2.compounds_targets.raw)) {
      data <- x$step2.compounds_targets.raw <- x@tables$step2$compounds_targets
    } else {
      data <- x$step2.compounds_targets.raw
    }
    data <- dplyr::mutate(data, inchikey2d = stringr::str_extract(InChIKey, "^[A-Z]+"))
    data <- dplyr::filter(data, inchikey2d %in% !!ref$match$inchikey2d)
    x@tables$step2$compounds_targets <- data
    return(x)
  })

setMethod("asjob_classyfire", signature = c(x = "job_tcmsp"),
  function(x){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    allcpds <- x@tables$step2$ingredients
    allcpds <- dplyr::filter(allcpds, !is.na(InChIKey))
    allcpds <- dplyr::mutate(allcpds, inchikey2d = gs(InChIKey, "^([A-Z]+).*", "\\1"))
    x <- job_classyfire(allcpds$inchikey2d, type = "inchikey")
    x$tcmsp_cpds <- allcpds
    return(x)
  })

.get_compounds_info.tcmsp <- function(x, data) {
  ids <- unique(data[[ "Mol ID" ]])
  if (!length(ids)) {
    message("No `Mol ID` to query.")
    return()
  }
  db <- new_db(paste0(x$savedir, "/compounds_id.rdata"), ".id")
  db <- not(db, ids)
  query <- db@query
  if (length(query)) {
    fun <- function() {
      urls <- data[[ "Molecule Name.link" ]]
      urls <- urls[ match(query, ids) ]
      nl(query, urls)
    }
    link <- start_drive()
    Sys.sleep(3)
    link$open()
    res <- pbapply::pblapply(fun(),
      function(url) {
        res <- try(link$navigate(url), T)
        if (inherits(res, "try-error")) {
          usethis::ui_yeah("Continue?")
        }
        table <- as_tibble(.get_current_tables(link)[[1]])
        dplyr::filter(dplyr::select(table, V1, V2),
          V1 %in% c("Molecule ID", "InChIKey"))
      })
    link$close()
    end_drive()
    res <- frbind(res, idcol = T)
    res <- tidyr::spread(res, V1, V2)
    db <- upd(db, res)
  }
  dplyr::filter(db@db, .id %in% ids)
}

.list_allfiles <- function(dir, pattern = ".") {
  locals <- list.files(dir, pattern, full.names = T, recursive = T)
  names(locals) <- get_realname(locals)
  locals
}
# UniProt.ws::mapUniProt

get_tcmsp_data <- function(update = F, savedir = "../tcmsp") {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  url_base <- "https://tcmsp-e.com/browse.php?qc="
  used <- c("herbs", "ingredients", "targets", "diseases")
  isDriveStarted <- F
  link <- F
  token <- NULL
  db <- sapply(used, simplify = F,
    function(name) {
      file <- paste0(savedir, "/", name, ".tsv")
      if (!isDriveStarted) {
        link <<- start_drive(browser = "firefox")
        Sys.sleep(3)
        isDriveStarted <<- T
        link$open()
        # link$setTimeout(type = "page load", milliseconds = 3000000)
      }
      url <- paste0(url_base, name)
      if (name == "herbs") {
        res <- try(link$navigate(url), T)
        if (inherits(res, "try-error")) {
          res <- usethis::ui_yeah("The loading time out. Continue?")
          if (!res) {
            stop("Stop.")
          }
        }
        isUrlNavigated <- T
        fun_tryGetToken <- function(link) {
          table <- tibble::as_tibble(.get_current_tables.tcmsp(link)[[2]])
          table <- trySepCols.tcmsp(table)
          token <- unique(gs(table$V2.link, ".*?([a-z0-9]+)$", "\\1"))
          if (length(token) == 1) {
            message("Get token:\n\t", token)
            Sys.sleep(1)
          } else {
            token <- NULL
          }
          token
        }
        token <- fun_tryGetToken(link)
        if (!is.null(token)) {
          token <<- token
        }
      } else {
        isUrlNavigated <- F
      }
      if (!file.exists(file) || update) {
        if (!isUrlNavigated) {
          res <- try(link$navigate(url), T)
        }
        nextFun <- function() {
          obj <- .get_folded_tables.tcmsp(link)
          obj <- .try_format_folded_tables.tcmsp(obj)
          write_tsv(obj, file)
          obj
        }
        if (inherits(res, "try-error")) {
          options(nextFun = nextFun)
          stop("RSelenium Timeout")
        } else {
          nextFun()
        }
      } else {
        ftibble(file)
      }
    })
  if (isDriveStarted) {
    link$close()
    end_drive()
  }
  if (!is.null(token)) {
    db <- lapply(db,
      function(x) {
        dplyr::mutate(x, dplyr::across(dplyr::ends_with(".link"),
            function(col) {
              gs(col, "token=[a-z0-9]+$", paste0("token=", token))
            }))
      })
    message("\n\tToken set successful: ", token, "\n")
  }
  db
}

trySepCols.tcmsp <- function(obj, sep = " ### ") {
  testCols <- function(obj) {
    vals <- unlist(obj[1, ], use.names = F)
    grp(vals, sep)
  }
  colsLink <- testCols(obj)
  colnames <- colnames(obj)
  if (length(colsLink)) {
    for (i in colsLink) {
      obj <- tidyr::separate(obj, !!rlang::sym(colnames[i]),
        into = paste0(colnames[i], c("", ".link")), sep = sep)
    }
  }
  obj
}

tryGetLink.tcmsp <- function(x, sep = " ### ", ...) {
  val <- XML::xmlValue(x)
  if (!is.null(role <- XML::xmlAttrs(x)[[ "role" ]])) {
    if (role != "gridcell") {
      return(val)
    }
  }
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

.try_format_folded_tables.tcmsp <- function(obj) {
  colnames <- colnames(obj[[1]][[1]])
  obj <- lapply(obj, function(x) x[[2]])
  obj <- frbind(obj, fill = T)
  colnames(obj) <- colnames
  obj <- trySepCols.tcmsp(obj)
  obj <- dplyr::mutate(obj,
    dplyr::across(dplyr::ends_with(".link"),
      function(x) paste0("https://tcmsp-e.com/", gs(x, " ", "%20"))))
  obj
}

.get_current_tables.tcmsp <- function(link) {
  .get_current_tables(link, tryGetLink.tcmsp)
}

.get_current_tables <- function(link, elFun = XML::xmlValue) {
  html <- link$getPageSource()[[1]]
  html <- XML::htmlParse(html)
  XML::readHTMLTable(html, elFun = elFun)
}

.get_folded_tables.tcmsp <- function(link, div = NULL, which = NULL)
{
  ## div is paramater for `xpath` find the specific div, e.g., '[@id="test"]'
  end <- F
  lst <- list()
  n <- 0L
  while (!end) {
    table <- .get_current_tables.tcmsp(link)
    if (!is.null(which)) {
      table <- table[ which ]
    }
    n <- n + 1L
    lst[[ n ]] <- table
    path <- paste0("//div", div, "//a[@title='Go to the next page']")
    ele <- try(link$findElement(using = "xpath", value = path), T)
    if (inherits(ele, "try-error")) {
      stop("Can not find the element of path: ", path)
    }
    class <- ele$getElementAttribute("class")
    if (class[[1]] == "k-link k-state-disabled") {
      end <- T
    } else {
      ele$sendKeysToElement(list("Go to the next page", key = "enter"))
      Sys.sleep(.1)
    }
  }
  lst
}
