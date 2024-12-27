# ==========================================================================
# workflow of herb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_herb <- setClass("job_herb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@HerbAHighThFang2021]",
    method = "Website `HERB` <http://herb.ac.cn/> used for TCM data source",
    tag = "pharm, herb:HERB",
    analysis = "HERB 网络药理学"
    ))

job_herb <- function(herbs, db = get_herb_data())
{
  x <- .job_herb(object = db)
  x@params$herbs <- herbs
  herbs_info <- filter(object(x)$herb, Herb_cn_name %in% !!params(x)$herbs)
  herbs_info <- .set_lab(herbs_info, "herbs information")
  x@params$herbs_info <- herbs_info
  print(herbs_info)
  message("Got the herbs:\n\n\t", paste0(herbs_info$Herb_cn_name, collapse = ", "),
      "\n\n", "\tTotal: ", colSum(herbs_info$Herb_cn_name))
  message("By modify `x@params$herbs_info` to filter the results.")
  x
}

setMethod("step0", signature = c(x = "job_herb"),
  function(x){
    step_message("Prepare your data with function `job_herb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_herb"),
  function(x, filter.hob = TRUE, filter.dl = TRUE, ...,
    tempdir = "download", db = .prefix("herb/herbs_ingredient.rds", "db"))
  {
    step_message("Dowload compounds of herbs.")
    ids <- params(x)$herbs_info$Herb_
    db <- new_db(db, "herb_id")
    db <- not(db, ids)
    if (length(db@query)) {
      x$tempdir <- tempdir
      link <- start_drive(browser = "firefox", download.dir = x$tempdir)
      Sys.sleep(3)
      tryCatch(
        download_herbCompounds(link, db@query, tempdir = x$tempdir),
        finally = end_drive()
      )
      fun_clear <- function() {
        if (file.exists("herbs_ingredient")) {
          unlink("herbs_ingredient", TRUE, TRUE)
        }
      }
      fun_clear()
      herbs_compounds <- collateFiles_herbs(db@query, from = x$tempdir, to = "herbs_ingredient")
      db <- upd(db, herbs_compounds, db@query)
      fun_clear()
      unlink(tempdir, TRUE, TRUE)
    }
    herbs_compounds <- dplyr::filter(db@db, herb_id %in% ids)
    herbs_compounds <- .set_lab(herbs_compounds, sig(x), "Components of Herbs")
    if (filter.hob || filter.dl) {
      smiles <- dplyr::filter(object(x)$component,
        Ingredient_id %in% !!unique(herbs_compounds$Ingredient.id),
        Ingredient_Smile != "Not Available")
      smiles <- dplyr::select(smiles, Ingredient_id, Ingredient_Smile)
      smiles <- nl(smiles$Ingredient_id, smiles$Ingredient_Smile, FALSE)
      ## Check smiles available
      message("Use rcdk to check smiles validity.")
      notValide <- vapply(rcdk::parse.smiles(unname(smiles)), is.null, logical(1))
      if (any(notValide)) {
        message("Some smiles was not valid: ", length(which(notValide)))
      }
      smiles <- smiles[ !notValide ]
      lst <- filter_hob_dl(smiles, filter.hob, filter.dl, use.cids = FALSE, ...)
      if (length(lst)) {
        p.upset <- .set_lab(lst$p.upset, sig(x), "Prediction of HOB and Drug-Likeness")
        x@plots[[ 1 ]] <- namel(p.upset)
        herbs_compounds <- dplyr::filter(herbs_compounds, Ingredient.id %in% !!lst$ids)
      }
    }
    x@tables[[ 1 ]] <- namel(herbs_compounds)
    return(x)
  })

setMethod("step2", signature = c(x = "job_herb"),
  function(x, group_number = 50, group_sleep = 3,
    searchDir = .prefix(), db = .prefix("herb/compounds_target.rds", "db"),
    removeExistsTemp = FALSE)
  {
    step_message("Dowload targets of compounds")
    all_ids <- unique(x@tables$step1$herbs_compounds$Ingredient.id)
    db <- new_db(db, "Ingredient_id")
    db <- not(db, all_ids)
    if (length(db@query)) {
      queries <- db@query
      merge_data <- FALSE
      if (!is.null(searchDir)) {
        message("The feature of `searchDir` was deprecated.")
        dbfile <- list.files(searchDir, "^HBIN.*.xlsx$", full.names = TRUE, recursive = TRUE)
        names(dbfile) <- get_realname(dbfile)
        dbfile <- dbfile[ !duplicated(dbfile) ]
        has_ids <- queries[ queries %in% names(dbfile) ]
        if (length(has_ids)) {
          files <- dbfile[ names(dbfile) %in% has_ids ]
          data <- collateFiles_herbs(readFrom = files, .id = "Ingredient_id", from = x$tempdir)
          queries <- queries[ !queries %in% has_ids ]
          merge_data <- TRUE
        }
      }
      if (length(queries)) {
        queries <- grouping_vec2list(queries, group_number, TRUE)
        link <- start_drive(browser = "firefox")
        Sys.sleep(3)
        lst <- pbapply::pblapply(queries,
          function(ids) {
            .dir <- "compounds_target"
            to <- paste0(.dir, "/", attr(ids, "name"))
            if (file.exists(to)) {
              .dir <- timeName("compounds_target")
              to <- paste0(.dir, "/", attr(ids, "name"))
            }
            download_compoundTargets(link, ids, tempdir = x$tempdir)
            data <- collateFiles_herbs(ids, to = to, .id = "Ingredient_id", from = x$tempdir)
            cli::cli_alert_info("Sleep between groups ...")
            Sys.sleep(group_sleep)
            data
          })
        lst <- data.table::rbindlist(lst, fill = TRUE)
        end_drive()
      } else {
        lst <- data.frame()
      }
      if (merge_data) {
        if (!is(data, "list")) {
          data <- list(data)
        }
        if (nrow(lst)) {
          lst <- c(list(lst), data)
        } else {
          lst <- data
        }
      }
      if (is(lst, "list")) {
        data <- data.table::rbindlist(lst, fill = TRUE)
      } else {
        data <- lst
      }
      compounds_targets <- as_tibble(data)
      db <- upd(db, compounds_targets, db@query)
    }
    compounds_targets <- dplyr::filter(db@db, Ingredient_id %in% !!all_ids)
    x@tables[[ 2 ]] <- namel(compounds_targets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_herb"),
  function(x, disease = NULL, disease.score = 5,
    mart_dataset = "hsapiens_gene_ensembl", HLs = NULL)
  {
    step_message("Get disease targets from genecards and annotate targets with biomaRt.
      As well, plot the intersection of herbs targets.
      "
    )
    if (is.null(x@params$mart)) {
      mart <- new_biomart(mart_dataset)
      x@params$mart <- mart
    } else {
      mart <- x@params$mart
    }
    ## annotate compounds targets
    message("Annotate compounds targets")
    compounds_targets <- x@tables$step2$compounds_targets
    targets <- unique(compounds_targets$Target.name)
    targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol", targets)
    x@params$targets_annotation <- targets_annotation
    compounds_targets_annotation <- tbmerge(
      compounds_targets, targets_annotation,
      by.x = "Target.name", by.y = "hgnc_symbol", all.x = TRUE,
      allow.cartesian = TRUE
    )
    ## create data: herbs_target
    message("Create data: herbs_target")
    herbs_compounds <- x@tables$step1$herbs_compounds
    if (!identical(class(herbs_compounds$Ingredient.id), class(compounds_targets$Ingredient_id))) {
      message("'ID' not identical in class, convert as character.")
      herbs_compounds$Ingredient.id %<>% as.character()
      compounds_targets$Ingredient_id %<>% as.character()
    }
    herbs_targets <- tbmerge(
      herbs_compounds, compounds_targets,
      by.x = "Ingredient.id", by.y = "Ingredient_id", all.x = TRUE,
      allow.cartesian = TRUE
    )
    x@params$herbs_targets <- herbs_targets
    easyRead <- map(x@params$herbs_targets, "herb_id", object(x)$herb, "Herb_", "Herb_pinyin_name")
    easyRead <- .set_lab(easyRead, sig(x), "Herbs compounds and targets")
    x@params$easyRead <- easyRead
    if (length(unique(herbs_targets$herb_id)) > 1) {
      ## plot upset of herbs targets
      sets <- dplyr::distinct(herbs_targets, herb_id, Target.name)
      sets <- split(sets$Target.name, sets$herb_id)
      fun <- function(sets) x@params$herbs_info$Herb_pinyin_name[match(names(sets), x@params$herbs_info$Herb_)]
      names(sets) <- fun(sets)
      p.herbs_targets <- new_upset(lst = sets)
      p.herbs_targets <- .set_lab(p.herbs_targets, sig(x), "Intersection of herbs all targets")
      p.herbs_targets <- setLegend(p.herbs_targets, "UpSet 图展示了中药的所有靶点 (中药名为拼音) 。该靶点是通过加和所有中药成分的靶点实现的。{try_summary(sets)}")
      ## plot upset of herbs ingredient
      sets <- split(herbs_compounds$Ingredient.id, herbs_compounds$herb_id)
      names(sets) <- fun(sets)
      Terror <<- sets
      p.herbs_compounds <- new_upset(lst = sets)
      p.herbs_compounds <- .set_lab(p.herbs_compounds, sig(x), "Intersection of herbs compounds")
      p.herbs_compounds <- setLegend(p.herbs_compounds, "UpSet 图展示了中药各组成成分之间的交集 (中药名为拼音)。左侧横柱状图展示了各中药对应化合物数量。右侧展示了交集数目。各中药的化合物数目为数据库包含的化合物数量。{try_summary(sets)}")
    } else {
      p.herbs_targets <- NULL
      p.herbs_compounds <- NULL
    }
    ## set plots and tables
    plots <- namel(p.herbs_targets, p.herbs_compounds)
    tables <- namel(compounds_targets_annotation)
    if (!is.null(disease)) {
      x$disease <- disease
      disease_targets <- get_from_genecards(disease, score = disease.score)
      .add_internal_job(.job_genecard())
      x@params$disease_targets <- disease_targets
      disease_targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol",
        disease_targets$Symbol)
      .add_internal_job(.job_biomart())
      tables <- c(tables, namel(disease_targets_annotation))
      lst <- list(
        targets = targets,
        targets_annotation = targets_annotation$hgnc_symbol,
        disease_targets = disease_targets$Symbol,
        disease_targets_annotation = disease_targets_annotation$hgnc_symbol
      )
      p.targets_disease <- new_upset(lst = lst)
      plots <- c(plots, namel(p.targets_disease))
      x@params$ppi_used <- dplyr::filter(
        targets_annotation,
        hgnc_symbol %in% dplyr::all_of(disease_targets_annotation$hgnc_symbol)
      )
    }
    data.allu <- dplyr::select(easyRead, Herb_pinyin_name, Ingredient.name, Target.name)
    data.allu <- dplyr::mutate(data.allu, Target.name = ifelse(is.na(Target.name), "Unkown", Target.name))
    data.allu <- dplyr::distinct(data.allu)
    x$data.allu <- data.allu
    p.pharm <- plot_network.pharm(data.allu, seed = x$seed, HLs = HLs)
    p.pharm <- .set_lab(p.pharm, sig(x), "network pharmacology visualization")
    plots <- c(plots, namel(p.pharm))
    x@plots[[ 3 ]] <- plots
    x@tables[[ 3 ]] <- tables
    return(x)
  })

# setMethod("merge", signature = c(x = "job_herb", y = "list"),
#   function(x, y){
#   })


rstyle <- function(get = "pal", seed = NULL, n = 1L) {
  if (get == "pal") {
    set <- lapply(ggsci:::ggsci_db[1:19], function(x) unname(x[[1]]))
  } else if (get == "shape.c") {
    set <- as.list(15:19)
  } else if (get == "shape.f") {
    set <- as.list(21:25)
  } else if (get == "theme") {
    set <- list(theme_minimal(), theme_bw(), theme_light())
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  res <- sample(set, n)
  if (n == 1) {
    return(res[[ 1 ]])
  }
  res
}

get_herb_data <- function(herb = .prefix("HERB_herb_info.txt", "db"),
  component = .prefix("HERB_ingredient_info.txt", "db"),
  target = .prefix("HERB_target_info.txt", "db"))
{
  db <- lapply(namel(herb, component, target),
    function(file) {
      tibble::as_tibble(suppressWarnings(data.table::fread(file)))
    })
  names <- colnames(db$component)
  colnames(db$component) <- c(names[-1], "extra_id")
  db
}

download_herbTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Gene Targets", ...)
{
  download_herbCompounds(link, ids, heading = heading, ...)
}

download_compoundTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Ingredient"),
  heading = "Related Gene Targets", ...)
{
  download_herbCompounds(link, ids, urls, heading = heading, ...)
}

download_herbCompounds <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Ingredients", tempdir = "download")
{
  get_all <- function() list.files(tempdir, "\\.xlsx$", full.names = TRUE)
  checks <- get_all()
  if (length(checks) > 0)
    file.remove(checks)
  link$open()
  jobn <- 0
  res <- lapply(urls,
    function(url) {
      link$navigate(url)
      jobn <<- jobn + 1
      if (jobn > 1) {
        link$refresh()
      }
      Sys.sleep(.5)
      ## check whether all the web page is all loaded
      isOK <- FALSE
      if (grpl(url, "Herb$")) {
        checkWhat <- "Herb id"
      } else if (grpl(url, "Ingredient$")) {
        checkWhat <- "Ingredient id"
      }
      while (!isOK) {
        status <- try(link$findElement(
            using = "xpath",
            value = paste0("//main//div//ul//span//b[text()='", checkWhat, "']")
            ), TRUE)
        if (!inherits(status, "try-error")) {
          isOK <- TRUE
        } else {
          message("Guess the page is loading...")
          Sys.sleep(1)
        }
      }
      res <- try(silent = TRUE, {
        ele <- link$findElement(
          using = "xpath",
          value = paste0("//main//div//h4[text()='", heading, "']/../div//button")
        )
        Sys.sleep(.5)
        ele$sendKeysToElement(list("Download", key = "enter"))
        })
      Sys.sleep(.5)
      if (inherits(res, "try-error")) {
        writeLines("", paste0(tempdir, "/empty_", jobn, ".xlsx"))
      }
      checks <- get_all()
      cli::cli_alert_info(paste0("Job number ", jobn))
      ## sometimes these is no response
      while (length(checks) < jobn) {
        cli::cli_alert_info("Retry remotes response ...")
        ele$sendKeysToElement(list("Download", key = "enter"))
        Sys.sleep(1)
        checks <- get_all()
      }
    })
  link$close()
}

start_drive <- function(command = paste0("java -jar ", .prefix("/selenium.jar", "op")),
  port = 4444, extra = NULL, browser = c("firefox", "chrome"), download.dir = "download", ...)
{
  system(paste(command, "-port", port, extra), wait = FALSE)
  new_link(port = port, browser = match.arg(browser), download.dir = download.dir, ...)
}

end_drive <- function(pattern = "[0-9]\\s*java.*-jar") {
  if (.Platform$OS.type == "unix") {
    if (getOption("end_drive", FALSE)) {
      tmp <- tempfile(fileext = ".txt")
      system(paste0("ps aux > ", tmp))
      text <- readLines(tmp)
      text <- text[ grepl(pattern, text) ]
      pid <- stringr::str_extract(text, "(?<=\\s)[0-9]{1,}(?=\\s)")
      system(paste("kill", paste0(pid, collapse = " ")))
    }
  }
}

get_table.html <- function(x, ...) {
  if (file.exists(x)) {
    str <- readLines(x)
  } else {
    str <- x
  }
  ht <- e(XML::htmlParse(str))
  XML::readHTMLTable(ht, ...)
}

collateFiles_herbs <- function(ids,
  file.pattern = "\\.xlsx$", index.pfun = file_seq.by_time,
  from = "download", to = "herbs_ingredient", suffix = ".xlsx", .id = "herb_id", readFrom = NULL)
{
  if (is.null(readFrom)) {
    if (!file.exists(to)) {
      args <- as.list(environment())
      args$readFrom <- NULL
      files <- do.call(collateFiles, args)
      isOrdered <- TRUE
    } else {
      files <- list.files(to, file.pattern, full.names = TRUE)
      isOrdered <- FALSE
    }
  } else {
    files <- readFrom
    isOrdered <- FALSE
  }
  data <- lapply(files,
    function(file) {
      res <- try(openxlsx::read.xlsx(file), TRUE)
      if (inherits(res, "try-error")) {
        tibble::tibble()
      } else res
    })
  if (isOrdered) {
    names(data) <- ids
  } else {
    names(data) <- get_realname(files)
  }
  data <- frbind(data, idcol = .id, fill = TRUE)
  if (nrow(data)) {
    data <- dplyr::relocate(data, !!rlang::sym(.id))
  }
  data
}

collateFiles <- function(ids,
  file.pattern, index.pfun = file_seq.by_time,
  from, to, suffix, ...)
{
  files <- list.files(from, file.pattern, full.names = TRUE)
  index <- index.pfun(files)
  files <- files[order(index)][seq_along(ids)]
  files <- rev(files)
  dir.create(to, FALSE, TRUE)
  files <- lapply(seq_along(ids),
    function(n) {
      file.copy(files[ n ], file <- paste0(to, "/", ids[n], suffix), TRUE)
      file
    })
  files
}

file_seq.by_time <- function(files) {
  info <- fs::file_info(files)
  info <- dplyr::mutate(info, .INDEX = seq_len(nrow(info)))
  info <- dplyr::arrange(info, dplyr::desc(modification_time))
  info <- dplyr::mutate(info, .SEQ = seq_len(nrow(info)))
  info <- dplyr::arrange(info, .INDEX)
  info$.SEQ
}

format_index.by_num <- function(index) {
  index <- gsub("\\(|\\)", "", index)
  index <- ifelse(is.na(index), "0", index)
  as.integer(index)
}

new_link <- function(port = 4444L, browser = c("firefox", "chrome"), addr = "localhost", download.dir = "download")
{
  if (!dir.exists(download.dir)) {
    dir.create(download.dir)
  }
  browser <- match.arg(browser)
  if (browser == "firefox") {
    download.dir <- normalizePath(download.dir)
    # about:config
    if (!is.null(path_prof <- getOption("FirefoxProfile"))) {
      prof <- RSelenium::getFirefoxProfile(path_prof)
    } else {
      prof <- RSelenium::makeFirefoxProfile(
        list('permissions.default.image' = 2L,
          'browser.download.folderList' = 2L,
          'browser.download.manager.showWhenStarting' = FALSE,
          'browser.download.lastDir' = download.dir,
          'browser.download.dir' = download.dir
        )
      )
    }
    link <- RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = prof
    )
    link
  } else if (browser == "chrome") {
    RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = list()
    )
  }
}

merge.componentsGenes <- function(data, genes.lst) {
  names(genes.lst) <- data[[1]]
  genes.lst <- lapply(genes.lst,
    function(lst) {
      if (nrow(lst$ids) == 0) {
        return(NULL)
      } else {
        lst$ids$genes <- lst$data
        return(lst$ids)
      }
    })
}

get_tcm.components <- function(id, key = "Components",
  link_prefix = "http://www.tcmip.cn/ETCM/index.php/Home/Index/yc_details.html?id=",
  src = NULL)
{
  if (is.null(src)) {
    src <- get_tcm.base(id, link_prefix)
  }
  pos <- grep(paste0("^\\s*\\{\"ID*.*", key), src)
  data <- src[pos]
  ids.data <- extract_id(data)
  data <- gsub("<[^<^>]*>", "##", data)
  data <- gsub("^.*?####|####[^#]*$", "", data)
  data <- strsplit(data, "##, ##")[[1]]
  list(data = data, ids = ids.data, src = src)
}

extract_id <- function(ch) {
  links <- stringr::str_extract_all(ch, "(?<=href=\').*?(?=\')")[[1]]
  ids <- stringr::str_extract(links, "[^=]*$")
  data.frame(links = links, ids = ids)
}

get_tcm.componentToGenes <- function(ids, key = "Candidate Target Genes", dep = TRUE) 
{
  link_prefix <- "http://www.tcmip.cn/ETCM/index.php/Home/Index/cf_details.html?id="
  lst <- pbapply::pblapply(ids,
    function(id) {
      data <- get_tcm.components(id, key, link_prefix)
      if (dep) {
        data$src <- NULL
      }
      data
    })
  lst
}

get_tcm.base <- function(id, link_prefix) {
  link <- paste0(link_prefix, id)
  src <- xml2::read_html(link)
  data <- rvest::html_text(src)
  data <- strsplit(data, "\r\n")[[1]]
  data
}

get_coords.sin <- function(num = 100, nWave = 4) {
  x <- seq(0, nWave * 2 * pi, length.out = num)
  tibble::tibble(x = x, y = sin(x))
}

get_coords.spiral <- function(num = 100, nCir = 3, minRad = 1, maxRad = 3) {
  ang <- seq(0, by = 2 * pi * nCir / num, length.out = num)
  rad <- seq(minRad, maxRad, length.out = num)
  tibble::tibble(x = sin(ang) * rad, y = cos(ang) * rad)
}

