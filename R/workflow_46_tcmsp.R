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

job_tcmsp <- function()
{
  .job_tcmsp()
}

setMethod("step0", signature = c(x = "job_tcmsp"),
  function(x){
    step_message("Prepare your data with function `job_tcmsp`.")
  })

setMethod("step1", signature = c(x = "job_tcmsp"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

# UniProt.ws::mapUniProt

get_tcmsp_data <- function(update = F, savedir = "../tcmsp") {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  url_base <- "https://tcmsp-e.com/browse.php?qc="
  used <- c("herbs", "ingredients", "targets", "diseases")
  isDriveStarted <- F
  link <- F
  db <- sapply(used, simplify = F,
    function(name) {
      file <- paste0(savedir, "/", name, ".tsv")
      if (!file.exists(file) || update) {
        if (!isDriveStarted) {
          link <<- start_drive(browser = "firefox")
          Sys.sleep(3)
          isDriveStarted <<- T
          link$open()
        }
        url <- paste0(url_base, name)
        link$navigate(url)
        lst <- .get_folded_tables.tcmsp(link)
        lst
      } else {
        ftibble(file)
      }
    })
  if (isDriveStarted) {
    end_drive()
  }
  db
}

.get_folded_tables.tcmsp <- function(link) {
  end <- F
  lst <- list()
  n <- 1L
  while (!end) {
    html <- link$getPageSource()[[1]]
    html <- XML::htmlParse(html)
    table <- XML::readHTMLTable(html)
    n <- n + 1L
    lst[[ n ]] <- table
    path <- "//div//a[@title='Go to the next page']"
    ele <- try(link$findElement(using = "xpath", value = path), T)
    if (inherits(ele, "try-error")) {
      stop("Can not find the element of path: ", path)
    }
    class <- ele$getElementAttribute("class")
    if (class[[1]] == "k-link k-state-disabled") {
      end <- T
    } else {
      ele$sendKeysToElement(list("Go to the next page", key = "enter"))
      Sys.sleep(.2)
    }
  }
  lst
}
