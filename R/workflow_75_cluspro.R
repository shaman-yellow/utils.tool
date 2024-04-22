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
    x$usr <- usr
    return(x)
  })

setMethod("pull", signature = c(x = "job_cluspro"),
  function(x, ref = NULL, dir = "cluspro", dir_download = "download")
  {
    message("Get the models and scores for finished job.")
    if (is.null(x$link)) {
      stop("is.null(x$link)")
    }
    link <- x$link
    url <- "https://cluspro.bu.edu/results.php"
    link$navigate(url)
    fun_format <- function(x) {
      colnames(x) <- unlist(x[1, ])
      x[-1, ]
    }
    history <- fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
    if (!is.null(x$groups) && !is.null(x$whichGroup)) {
      group <- x$groups[[ x$whichGroup ]]
      group <- group[ !group %in% x$excludes ]
      status <- dplyr::filter(history, Name %in% !!group)
      message("Current finished: ", nrow(status), "/", length(group))
      if (is.null(ref)) {
        ref <- status$Name
      }
    }
    if (length(ref)) {
      status <- dplyr::filter(history, Name %in% !!ref)
    } else {
      status <- data.frame()
    }
    if (nrow(status)) {
      dir.create(dir, F)
      fun_format <- function(x) {
        cluster <- 0L
        members <- 0L
        names <- colnames(x)
        x <- apply(x, 1, simplify = T,
          function(v) {
            if (grpl(v[[ "Cluster" ]], "^[0-9]+$")) {
              cluster <<- v[[ "Cluster" ]]
              members <<- v[[ "Members" ]]
              type <- T
            } else {
              type <- F
            }
            if (type) {
              c(cluster, members, v[3:4])
            } else {
              c(cluster, members, v[1:2])
            }
          })
        x <- data.frame(t(x))
        colnames(x) <- names
        x <- tidyr::spread(x, Representative, `Weighted Score`)
        x <- dplyr::mutate_all(x, as.double)
        x <- dplyr::arrange(x, Cluster)
      }
      if (!dir.exists(dir_download)) {
        stop("dir.exists(dir_download)")
      }
      res <- pbapply::pbapply(status, 1, simplify = F,
        function(v) {
          id <- v[[ "Id" ]]
          message("\nPull results: \n\t", v[[ "Name" ]])
          ## download model: 0
          if (file.exists(dfile <- paste0(dir_download, "/model.000.00.pdb"))) {
            file.remove(dfile)
          }
          url <- paste0("https://cluspro.bu.edu/models.php?job=", id)
          message("Navigate:\n\t", url)
          link$navigate(url)
          ele <- link$findElement("xpath",
            paste0("//td//a[@href='file.php?jobid=", id, "&coeffi=0&model=0&filetype=model_file']"))
          ele$clickElement()
          Sys.sleep(1)
          # R.utils::withTimeout(link$navigate(url), timeout = 1, onTimeout = "silent")
          subdir <- paste0(dir, "/", v[[ "Name" ]])
          dir.create(subdir, F)
          exist <- F
          file_model <- paste0(subdir, "/model.000.00.pdb")
          n <- 0L
          while (!exist && n < 10L) {
            n <- n + 1L
            Sys.sleep(1)
            exist <- file.copy(dfile, file_model)
          }
          ## score
          link$navigate(paste0("https://cluspro.bu.edu/scores.php?job=", id, "&coeffi=0"))
          score <- fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
          namel(score, file_model)
        })
      names(res) <- status$Name
      scores <- lapply(res, function(x) x$score)
      scores <- frbind(scores, idcol = "Name")
      scores <- dplyr::rename_all(scores, make.names)
      x$tmp_scores <- scores
      x$tmp_models <- vapply(res, function(x) x$file_model, character(1))
    } else {
      if (length(ref)) {
        message("Job maybe not finished: ", ref)
      } else {
        message("No job has been finished.")
      }
    }
    return(x)
  })

setMethod("upload", signature = c(x = "job_cluspro"),
  function(x){
    message("Collate current status and upload jobs.")
    notFinish <- which(!x$groups_status)
    if (!length(notFinish)) {
      message("All job has been finished.")
      return(x)
    }
    if (x$in_queue) {
      message("Please wait for completion of previous jobs.")
      return(x)
    }
    x$whichGroup <- whichGroup <- notFinish[1]
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
    link <- x$link
    group <- rawGroup <- x$groups[[ whichGroup ]]
    fun_get <- function(url) {
      link$navigate(url)
      fun_format <- function(x) {
        colnames(x) <- unlist(x[1, ])
        x[-1, ]
      }
      fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
    }
    ## check if got results
    url <- "https://cluspro.bu.edu/results.php"
    finished <- dplyr::filter(fun_get(url), Name %in% !!group)
    message("Has got results: ", nrow(finished), "/", length(group), "/", length(rawGroup))
    group <- group[ !group %in% finished$Name ]
    ## check now running
    url <- "https://cluspro.bu.edu/queue.php"
    queues <- dplyr::filter(fun_get(url), User == x$usr)
    runs <- group[ group %in% queues$Name ]
    message("Is now in queue: ", length(runs), "/", length(group), "/", length(rawGroup))
    group <- group[ !group %in% queues$Name ]
    ####################################################
    ####################################################
    ## upload
    if (length(group)) {
      message("Start upload.")
      n <- 0L
      excludes <- pbapply::pblapply(group,
        function(i) {
          message("Upload jobname: ", i)
          n <<- n + 1L
          if (n >= 1) {
            link$navigate(url)
          }
          layout <- x$layouts[[ i ]]
          if (any(is.na(layout))) {
            warning("No valid PDB id (NA), scape combining of ", i)
            Sys.sleep(3)
            return(i)
          }
          ele <- link$findElement("xpath", "//div//input[@name='jobname']")
          ele$sendKeysToElement(list(i))
          ele <- link$findElement("xpath", "//form//div//select[@name='server']//option[@value='gpu']")
          ele$clickElement()
          ele <- link$findElement("xpath", "//div//span//input[@id='ligpdb']")
          ele$sendKeysToElement(list(layout[[ "from" ]]))
          ele <- link$findElement("xpath", "//div//span//input[@id='recpdb']")
          ele$sendKeysToElement(list(layout[[ "to" ]]))
          ele <- link$findElement("xpath", "//div//input[@type='submit']")
          ele$clickElement()
          message("JOB has successfully upload: ", i)
          Sys.sleep(3)
        })
      x$excludes <- unlist(excludes)
    }
    x$groups_status[ whichGroup ] <- T
    x$in_queue <- T
    return(x)
  })


