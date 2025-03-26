# ==========================================================================
# workflow of tcmsp2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tcmsp2 <- setClass("job_tcmsp2", 
  contains = c("job_tcmsp"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "tcmsp2",
    info = c("Tutorial: https://www.tcmsp-e.com/#/home"),
    cite = "[@TcmspADatabaRuJi2014]",
    method = "Website `TCMSP` <https://tcmsp-e.com/> used for data source",
    tag = "db:tcmsp, pharm",
    analysis = "TCMSP 网络药理学"
    ))

job_tcmsp2 <- function(herbs)
{
  .job_tcmsp2(params = list(herbs = rm.no(herbs)))
}

setMethod("step0", signature = c(x = "job_tcmsp2"),
  function(x){
    step_message("Prepare your data with function `job_tcmsp2`.")
  })

setMethod("step1", signature = c(x = "job_tcmsp2"),
  function(x, forceGetHerbInfo = FALSE, dir = .prefix("tcmsp2", "db"), port = 8888, ...){
    step_message("Check local database.")
    dir.create(dir, FALSE)
    x$dir <- dir
    items <- c(
      # "disease" = "Related Diseases",
      "ingredients" = "Ingredients",
      "targets" = "Related Targets"
    )
    dbs <- sapply(names(items), simplify = FALSE,
      function(name) {
        db <- new_db(file.path(dir, paste0(name, ".rds")), "herb")
        db <- not(db, x@params$herbs)
      })
    queries <- dbs[[ "ingredients" ]]@query
    if (length(queries) || forceGetHerbInfo) {
      x$link <- start_drive(port = port, ...)
      x$link$open()
      x$url_home <- "https://www.tcmsp-e.com/#/"
      x$link$navigate(x$url_home)
    }
    dbHerbs <- new_db(file.path(dir, "herbs_info.rds"), "herb")
    x$dbHerbs <- not(dbHerbs, x$herbs)
    x$dbs <- dbs
    x$queries <- queries
    x$items <- items
    return(x)
  })

setMethod("step2", signature = c(x = "job_tcmsp2"),
  function(x, forceGetHerbInfo = FALSE)
  {
    step_message("Download data.")
    queries <- x$queries
    queries <- queries[ !queries %in% unique(x$dbs[[ "ingredients" ]]@db$herb) ]
    if (length(queries) || forceGetHerbInfo) {
      if (forceGetHerbInfo) {
        queries <- x$herbs
      }
      link <- x$link
      items <- x$items
      isLogin <- try(link$findElement("xpath", "//div//div[@id='login-btn']"), TRUE)
      if (!inherits(isLogin, "try-error")) {
        stop("You should login for TCMSP first.")
      } else {
        message("TCMSP has logined.")
      }
      lstHerbItems <- pbapply::pblapply(queries,
        function(herb) {
          ## download
          # >>>
          fun_input <- function() {
            link$findElement("xpath", "//div//input[@class='el-input__inner']")
          }
          fun_check <- function(info) {
            ele <- try(link$findElement("xpath", "//div//div[@class='mspTips']//div[@class='left']"), TRUE)
            if (inherits(ele, "try-error")) {
              return(FALSE)
            }
            text <- ele$getElementText()
            grpl(text, info$V3)
          }
          fun_pageSize <- function() {
            # change page size
            message("Switch page size.")
            ele <- link$findElement("xpath",
              paste0("//span[@class='el-pagination__sizes']",
                "//i[@class='el-icon el-select__caret el-select__icon']"))
            ele$clickElement()
            Sys.sleep(1)
            # span[contains(text(), '200条/页')]
            ele <- link$findElement("xpath", "//ul[@role='listbox']//li//span[text()='200条/页']")
            ele$clickElement()
          }
          fun_collate <- function() {
            lst <- get_table.html(link$getPageSource()[[1]])
            data <- lst[[2]]
            colnames(data) <- colnames(lst[[1]])
            data
          }
          fun_next <- function() {
            ## switch to next page
            ele <- link$findElement("xpath", "//div[@class='el-pagination']//button[@class='btn-next']")
            attr <- ele$getElementAttribute("aria-disabled")[[1]]
            if (attr == "true") {
              return(FALSE)
            } else {
              ele$clickElement()
              return(TRUE)
            }
          }
          fun_get <- function() {
            get_table.html(link$getPageSource()[[1]])
          }
          # >>>
          message("Try Dowload data of ", herb)
          ele <- try(fun_input(), TRUE)
          if (inherits(ele, "try-error")) {
            link$navigate(x$url_home)
            ele <- try(fun_input(), TRUE)
          }
          ele$clearElement()
          ele$sendKeysToElement(list(paste0(" ", herb, " ")))
          Sys.sleep(.1)
          ele <- link$findElement("xpath", "//div//span//img[@class='btn-search']")
          ele$clickElement()
          Sys.sleep(3)
          # check search page status
          info <- fun_get()
          n <- 0L
          empty <- FALSE
          ele <- try(link$findElement("xpath", "//div//div[@class='el-table__empty-block']"), TRUE)
          if (!inherits(ele, "try-error")) {
            empty <- TRUE
            warning("No results for matching: ", herb)
            return()
          }
          while ((is.null(info[[2]]) || !length(nrow(info[[2]]))) && (!empty)) {
            n <- n + 1L
            if (n > 5) {
              stop("Search failed.")
            }
            Sys.sleep(1)
            info <- fun_get()
          }
          # check search results
          infoNames <- colnames(info[[1]])
          info <- dplyr::filter(info[[2]], V1 == herb)
          if (!nrow(info)) {
            warning("No results for matching: ", herb)
            Sys.sleep(3)
            return()
          } else if (nrow(info) > 1) {
            message("Too many results found for: ", herb)
            print(info)
            which <- menu(info$V3, title = "Use which?")
            info <- info[which, ]
          }
          if (forceGetHerbInfo) {
            colnames(info) <- infoNames
            return(list(herbs_info = info))
          }
          ele <- link$findElement("xpath", paste0("//span[text()='", info$V3, "']"))
          ele$clickElement()
          Sys.sleep(3)
          # check whether content loaded
          n <- 0L
          while (!fun_check(info)) {
            if (n > 10) {
              warning("Something error for ", herb, ": ", info$V3)
              Sys.sleep(3)
              return()
            }
            n <- n + 1L
            Sys.sleep(1)
          }
          ## collate results 
          message("In frame of ", herb, ": ", info$V3)
          lstItems <- lapply(items,
            function(type) {
              message("In collating: ", type)
              # check whether focused
              active <- FALSE
              n <- 0L
              while (!active && n < 10L) {
                n <- n + 1L
                xpath <- paste0("//div[@class='el-radio-group']//span[text()='", type, "']")
                ele <- link$findElement("xpath", xpath)
                ele$clickElement()
                Sys.sleep(1)
                if (n > 5) {
                  system("notify-send 'Crawl TCMSP' 'Can not switch frame automaticly! Please manualy switch!'")
                  hasSwitch <- usethis::ui_yeah("Has the windows been switched?")
                  if (!hasSwitch) {
                    stop("Can not switch frame automaticly! Please manualy switch!")
                  }
                }
                # the parent frame is 'label'
                ele <- link$findElement("xpath", paste0(xpath, "/.."))
                attr <- ele$getElementAttribute("class")[[1]]
                active <- grpl(attr, "is-active")
              }
              if (!active) {
                stop("Can not focused.")
              }
              message("Is focused: ", type)
              if (FALSE) {
                Sys.sleep(1)
                fun_pageSize()
                Sys.sleep(2)
                page.size <- 200
              } else {
                page.size <- 50
              }
              ele <- link$findElement("xpath", "//div[@class='el-pagination']//span[@class='el-pagination__total is-first']")
              total <- as.integer(strx(ele$getElementText()[[1]], "[0-9]+"))
              len <- total %/% page.size + 1
              lstContents <- lapply(1:len,
                function(x) {
                  data <- fun_collate()
                  canNext <- fun_next()
                  Sys.sleep(2)
                  if (x == len && canNext) {
                    stop("Total Page number error.")
                  }
                  return(data)
                })
              dplyr::mutate(frbind(lstContents), herb = herb)
            })
          colnames(info) <- infoNames
          c(lstItems, list(herbs_info = info))
        })
      ## update local database
      herbs_info <- frbind(lapply(lstHerbItems, function(x) x$herbs_info))
      if (nrow(herbs_info)) {
        herbs_info <- dplyr::mutate(herbs_info, herb = `Chinese name`)
        if (forceGetHerbInfo) {
          x$herbs_info <- herbs_info
          return(x)
        }
        x$dbHerbs <- upd(x$dbHerbs, herbs_info)
        x$dbs <- sapply(names(x$items), simplify = FALSE,
          function(name) {
            data <- frbind(lapply(lstHerbItems, function(x) x[[ name ]]))
            upd(x$dbs[[ name ]], data)
          })
      }
    }
    lst <- lapply(x$dbs,
      function(db) {
        as_tibble(dplyr::filter(db@db, herb %in% !!x@params$herbs))
      })
    herbs_info <- dplyr::filter(x$dbHerbs@db, herb %in% !!x$herbs)
    tables <- list()
    tables$ingredients <- map(lst$ingredients,
      "herb", herbs_info, "herb", "Chinese Pinyin name",
      col = "Herb_pinyin_name"
    )
    tables$compounds_targets <- dplyr::distinct(lst$targets,
      `Mol ID`, `Molecule name`, `Target name`
    )
    x@tables[[ 2 ]] <- tables
    x@params$herbs_info <- dplyr::select(herbs_info,
      Herb_pinyin_name = `Chinese Pinyin name`,
      Herb_cn_name = `Chinese name`
    )
    ## remove NA
    x$dbs[[ "ingredients" ]]@db <- dplyr::filter(x$dbs[[ "ingredients" ]]@db, !is.na(`Mol ID`))
    x$dbHerbs@db <- dplyr::filter(x$dbHerbs@db, !is.na(`Chinese name`))
    x$herbs_notGet <- x$herbs[ !x$herbs %in% x$herbs_info$Herb_cn_name ]
    return(x)
  })

setMethod("step3", signature = c(x = "job_tcmsp2"),
  function(x, filter = TRUE, cutoff.ob = 30, cutoff.dl = 0.18, ...)
  {
    step_message("Call next method.")
    object(x) <- list()
    if (filter) {
      data <- x@tables$step2$ingredients
      data <- dplyr::mutate(data,
        `OB (%)` = as.double(`OB (%)`),
        DL = as.double(DL)
      )
      data <- dplyr::filter(data, `OB (%)` > cutoff.ob, DL > cutoff.dl)
      message("Filter: before (", nrow(x@tables$step2$ingredients), ") vs after (", nrow(data), ")")
      x@tables$step2$ingredients <- data
      x@tables$step2$compounds_targets <- dplyr::filter(
        x@tables$step2$compounds_targets,
        `Mol ID` %in% !!data[[ "Mol ID" ]]
      )
    }
    x <- callNextMethod(x, ...)
    data <- .set_lab(data, sig(x), "compounds filtered by OB and DL")
    attr(data, "lich") <- new_lich(list("OB (\\%) and DL cut-off" = "OB >= 30\\%; DL >= 0.18"))
    x@tables[[ 3 ]]$t.ob_dl_filtered <- data
    x$easyRead <- dplyr::select(x$easyRead, Herb_pinyin_name, `Molecule name`, symbols, protein.names)
    return(x)
  })
