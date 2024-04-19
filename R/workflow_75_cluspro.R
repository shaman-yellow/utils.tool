# ==========================================================================
# workflow of cluspro
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_cluspro <- setClass("job_cluspro", 
  contains = c("job_hawkdock"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "cluspro",
    info = c("https://cluspro.bu.edu/help.php"),
    cite = "[@TheClusproWebKozako2017]",
    method = "The `ClusPro` server used for Protein-Protein docking"
    ))

job_cluspro <- function(symbols, .layout = NULL)
{
  if (!is.null(.layout)) {
    message("Use columns of `.layout` as 'from' and 'to'.")
    .layout <- .layout[, 1:2]
    colnames(.layout) <- c("from", "to")
    symbols <- unlist(.layout, use.names = F)
  }
  .job_cluspro(object = rm.no(symbols), params = list(.layout = .layout))
}

setMethod("step0", signature = c(x = "job_cluspro"),
  function(x){
    step_message("Prepare your data with function `job_cluspro`.")
  })

setMethod("step1", signature = c(x = "job_cluspro"),
  function(x, extra_pdbs, ...){
    step_message("Found PDB.")
    x <- callNextMethod(x, extra_pdbs, ...)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cluspro"),
  function(x, port = 9999, ..., usr = getOption("cluspro_user"), psw = getOption("cluspro_password"))
  {
    step_message("Login.")
    x <- login(x, port, ..., usr = usr, psw = psw)
    return(x)
  })

setMethod("step3", signature = c(x = "job_cluspro"),
  function(x, n = 10){
    step_message("Prepare groups for upload.")
    layout <- x$layouts
    groups <- grouping_vec2list(names(layout), n, T)
    groups_status <- rep(F, length(groups))
    x$groups <- groups
    x$groups_status <- groups_status
    x$in_queue <- F
    message("Use `upload` to start the first jobs.")
    return(x)
  })

setMethod("login", signature = c(x = "job_cluspro"),
  function(x, port = 9999, ..., usr = getOption("cluspro_user"), psw = getOption("cluspro_password")){
    message("Login")
    x$link <- start_drive(port = port, ...)
    Sys.sleep(3)
    x$link$open()
    x$url_home <- "https://cluspro.bu.edu/"
    x$link$navigate(x$url_home)
    if (is.null(usr) || is.null(psw)) {
      message("Please enter the user name and password to login.")
    } else {
      ele <- x$link$findElement("xpath", "//div//form//input[@id='username']")
      ele$sendKeysToElement(list(usr))
      ele <- x$link$findElement("xpath", "//div//form//input[@id='password']")
      ele$sendKeysToElement(list(psw))
      ele <- x$link$findElement("xpath", "//div//form//input[@type='submit']")
      ele$clickElement()
      Sys.sleep(1)
    }
    return(x)
  })

setMethod("upload", signature = c(x = "job_cluspro"),
  function(x){
    notFinish <- which(!x$groups_status)
    if (!length(notFinish)) {
      message("All job has been finished.")
      return(x)
    }
    if (x$in_queue) {
      message("The server maybe overload, Please wait for completion of previous jobs.")
      return(x)
    }
    whichGroup <- notFinish[1]
    ## check link
    if (is.null(x$link)) {
      stop("is.null(x$link)")
    }
    check <- try(x@params$link$getCurrentUrl(), T)
    url <- "https://cluspro.bu.edu/home.php"
    if (inherits(check, "try-error")) {
      stop("Please check the web server is login.")
    } else if (check[[1]] != url) {
      x$link$navigate(url)
      Sys.sleep(1)
      check <- try(x@params$link$getCurrentUrl(), T)
      if (check[[1]] != url) {
        stop("Maybe not login.")
      }
    }
    message("Start upload.")

    x$groups_status[ whichGroup ] <- T
    return(x)
  })
