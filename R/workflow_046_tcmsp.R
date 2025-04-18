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
    method = "Website `TCMSP` <https://tcmsp-e.com/tcmsp.php> used for data source",
    tag = "db:tcmsp",
    analysis = "TCMSP 网络药理学"
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
  function(x, savedir = .prefix("tcmsp/largedata", "db"))
  {
    step_message("Dowload detail information of herbs.")
    if (!dir.exists(savedir)) {
      dir.create(savedir)
    }
    items <- c(
      "ingredients" = "Ingredients",
      "targets" = "Related Targets",
      "disease" = "Related Diseases"
    )
    isDriveStarted <- FALSE
    link <- FALSE
    db_file <- as.list(.list_allfiles(savedir, "^herb_"))
    n <- 0L
    lst <- sapply(x$herbs_info$Herb_pinyin_name, simplify = FALSE,
      function(name) {
        n <<- n + 1L
        fileid <- paste0("herb_", name, "_", names(items))
        files <- paste0(savedir, "/", fileid, ".tsv")
        fun <- function(n) is.null(db_file[[ fileid[n] ]])
        if (fun(1) | fun(2) | fun(3)) {
          link <<- start_drive(browser = "firefox")
          isDriveStarted <<- TRUE
          link$open()
          url <- x$herbs_info$`Latin name.link`[n]
          status <- try(link$navigate(url), TRUE)
          if (inherits(status, "try-error")) {
            res <- usethis::ui_yeah("Do you want to continue?")
            if (!res) {
              stop("Stop.")
            }
          }
          sets <- lapply(seq_along(items),
            function(i) {
              path <- paste0("//div//ul//li//a[text()='", items[i], "']")
              ele <- try(link$findElement(using = "xpath", value = path), TRUE)
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
    tables <- sapply(names(items), simplify = FALSE,
      function(name) {
        data <- frbind(lapply(lst, function(x) x[[name]]), idcol = TRUE, fill = TRUE)
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
    ## there is a terreble name rule break !!!  `Molecule Name` <-> `Molecule name`
    message("For `compounds_targets`: switch colname of `Molecule name` to `Molecule Name`.")
    targets <- dplyr::select(x@tables$step1$targets, `Mol ID`,
      `Molecule Name` = `Molecule name`, `Target name`)
    targets <- dplyr::distinct(targets)
    x@tables[[ 2 ]] <- namel(compounds_info, ingredients = data, compounds_targets = targets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_tcmsp"),
  function(x, disease = NULL, disease.score = 5, HLs = NULL,
    db_uniprot = .prefix("tcmsp/largedata/uniprotkb.rds", "db"), ...)
  {
    step_message("Query genes (symbol) for target proteins, and do step3 similar to `job_herb`...")
    ########################
    ########################
    compounds_targets <- x@tables$step2$compounds_targets
    if (!any(colnames(compounds_targets) == "symbols")) {
      cli::cli_alert_info("UniprotKB")
      kb <- job_uniprotkb(compounds_targets[[ "Target name" ]], db_uniprot)
      kb <- suppressMessages(step1(kb, ...))
      .add_internal_job(kb)
      res <- kb@tables$step1$format_results
      compounds_targets <- tbmerge(compounds_targets, res, by.x = "Target name", by.y = "query",
        all.x = TRUE, allow.cartesian = TRUE)
      x@tables$step2$compounds_targets <- compounds_targets
    }
    ########################
    ########################
    cli::cli_alert_info("step3: job_herb")
    hb <- .job_herb(step = 2L)
    hb@tables$step1$herbs_compounds <- dplyr::select(
      x@tables$step2$ingredients,
      herb_id = Herb_pinyin_name, Ingredient.id = `Mol ID`,
      Ingredient.name = `Molecule Name`
    )
    hb@tables$step2$compounds_targets <- dplyr::select(compounds_targets,
      Ingredient_id = `Mol ID`, Target.name = symbols,
      Target.protein = `Target name`
    )
    hb@params$herbs_info <- dplyr::select(
      dplyr::mutate(x@params$herbs_info, Herb_ = Herb_pinyin_name),
      Herb_, Herb_pinyin_name, Herb_cn_name
    )
    hb@object$herb <- hb@params$herbs_info
    hb <- suppressMessages(step3(hb, disease = disease, HLs = HLs))
    x@plots[[ 3 ]] <- c(kb@plots$step1, hb@plots$step3)
    easyRead <- tbmerge(
      dplyr::select(x@tables$step2$ingredients, `Mol ID`, Herb_pinyin_name),
      compounds_targets, by = "Mol ID",
      allow.cartesian = TRUE
    )
    x$easyRead <- dplyr::select(easyRead, - `Mol ID`)
    x@tables[[ 3 ]] <- namel(disease_targets_annotation = hb@tables$step3$disease_targets_annotation)
    x$data.allu <- hb$data.allu
    return(x)
  })

setMethod("map", signature = c(x = "job_tcmsp", ref = "job_classyfire"),
  function(x, ref, ...) {
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    if (is.null(x$step2.compounds_targets.raw)) {
      data <- x$step2.compounds_targets.raw <- x@tables$step2$compounds_targets
    } else {
      data <- x$step2.compounds_targets.raw
    }
    data <- dplyr::mutate(data, inchikey2d = ink2d(InChIKey))
    if (ref@step < 2L) {
      data <- dplyr::filter(data, inchikey2d %in% ref@params$db_inchi$inchikey2d, ...)
    } else {
      if (is.null(ref$match)) {
        stop("is.null(ref$match)")
      }
      data <- dplyr::filter(data, inchikey2d %in% !!ref$match$inchikey2d, ...)
    }
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
    link$open()
    res <- pbapply::pblapply(fun(),
      function(url) {
        res <- try(link$navigate(url), TRUE)
        if (inherits(res, "try-error")) {
          usethis::ui_yeah("Continue?")
        }
        table <- as_tibble(.get_current_tables(link)[[1]])
        dplyr::filter(dplyr::select(table, V1, V2),
          V1 %in% c("Molecule ID", "InChIKey"))
      })
    link$close()
    end_drive()
    res <- frbind(res, idcol = TRUE)
    res <- tidyr::spread(res, V1, V2)
    db <- upd(db, res)
  }
  dplyr::filter(db@db, .id %in% ids)
}

.list_allfiles <- function(dir, pattern = ".") {
  locals <- list.files(dir, pattern, full.names = TRUE, recursive = TRUE)
  names(locals) <- get_realname(locals)
  locals
}
# UniProt.ws::mapUniProt

get_tcmsp_data <- function(update = FALSE, savedir = .prefix("tcmsp", "db")) {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  url_base <- "https://tcmsp-e.com/browse.php?qc="
  used <- c("herbs", "ingredients", "targets", "diseases")
  isDriveStarted <- FALSE
  link <- FALSE
  token <- NULL
  db <- sapply(used, simplify = FALSE,
    function(name) {
      file <- paste0(savedir, "/", name, ".tsv")
      if (!isDriveStarted) {
        link <<- start_drive(browser = "firefox")
        isDriveStarted <<- TRUE
        link$open()
        # link$setTimeout(type = "page load", milliseconds = 3000000)
      }
      url <- paste0(url_base, name)
      if (name == "herbs") {
        res <- try(link$navigate(url), TRUE)
        if (inherits(res, "try-error")) {
          res <- usethis::ui_yeah("The loading time out. Continue?")
          if (!res) {
            stop("Stop.")
          }
        }
        isUrlNavigated <- TRUE
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
        isUrlNavigated <- FALSE
      }
      if (!file.exists(file) || update) {
        if (!isUrlNavigated) {
          res <- try(link$navigate(url), TRUE)
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
    vals <- unlist(obj[1, ], use.names = FALSE)
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
  obj <- frbind(obj, fill = TRUE)
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
  end <- FALSE
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
    ele <- try(link$findElement(using = "xpath", value = path), TRUE)
    if (inherits(ele, "try-error")) {
      stop("Can not find the element of path: ", path)
    }
    class <- ele$getElementAttribute("class")
    if (class[[1]] == "k-link k-state-disabled") {
      end <- TRUE
    } else {
      ele$sendKeysToElement(list("Go to the next page", key = "enter"))
      Sys.sleep(.1)
    }
  }
  lst
}

setMethod("filter", signature = c(x = "job_tcmsp"),
  function(x, ..., cutoff = 30){
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    message("Only `ingredients` data in `step1` was modified.")
    data <- x@tables$step2$ingredients <- dplyr::filter(
      x@tables$step2$ingredients, ...
    )
    x@tables$step2$compounds_targets <- dplyr::filter(
      x@tables$step2$compounds_targets, `Mol ID` %in% !!data$`Mol ID`
    )
    p <- ggplot(data) +
      geom_col(aes(x = reorder(`Molecule Name`, `OB (%)`), y = `OB (%)`, fill = DL), width = .7) +
      geom_hline(yintercept = cutoff, linetype = 4) +
      coord_flip() +
      labs(x = "Molecule Name", y = "OB (%)", fill = "DL") +
      theme_minimal()
    p <- wrap(p, 7, height = 2 + nrow(data) * .3)
    p <- .set_lab(p, sig(x), "Filterd ingredients")
    x$p.filtered <- p
    x@step <- 2L
    x
  })
