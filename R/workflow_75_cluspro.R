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
    method = "The `ClusPro` server used for Protein-Protein docking",
    tag = "dock:protein+protein",
    analysis = "ClusPro 蛋白质-蛋白质对接预测"
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
    message("Use `upload` to start the first jobs.")
    return(x)
  })

setMethod("step4", signature = c(x = "job_cluspro"),
  function(x, top = 3)
  {
    step_message("Plot scores and visualize models.")
    if (!is.null(x$tmp_scores)) {
      data <- dplyr::filter(x$tmp_scores, Cluster == 0)
      data <- tidyr::separate(data, "Name", c("pro1", "pro2"), "_", remove = F)
      recodes <- as.list(x$pdb_from)
      data <- dplyr::mutate(data, value = Lowest.Energy,
        pdb1 = dplyr::recode(pro1, !!!recodes),
        pdb2 = dplyr::recode(pro2, !!!recodes),
        var = paste0(pro1, " (", pdb1, ") + ", pro2, " (", pdb2, ")")
      )
      t.data <- data <- dplyr::arrange(data, Lowest.Energy)
      p.score <- ggplot(data, aes(x = reorder(var, value, decreasing = TRUE), y = value, fill = value)) +
        geom_col(width = .5) +
        geom_text(aes(x = var, y = value + max(value) * .02, label = value), hjust = 1, size = 3) +
        ylim(c(min(data$value) * 1.2, 0)) +
        coord_flip() +
        labs(y = "Lowest Energy of Top", x = "Proteins (PDB)") +
        rstyle("theme") +
        theme(legend.position = "")
      p.score <- wrap(p.score, 9, nrow(data) * .3 + .5)
      p.score <- .set_lab(p.score, sig(x), "Overview of protein docking results.")
      plots <- namel(p.score)
      if (!is.null(top)) {
        data <- head(data, top)
        n <- 0L
        figs <- pbapply::pbapply(data, 1, simplify = F,
          function(v) {
            n <<- n + 1L
            name <- v[[ "Name" ]]
            file <- vis_pdb.cluspro(
              x$tmp_models[[ name ]],
              paste0(v[[ "pro2" ]], " (", v[[ "pdb2" ]], ")"),
              paste0(v[[ "pro1" ]], " (", v[[ "pdb1" ]], ")")
            )
            fig <- file_fig(paste0("Top_", n, "_", name), file)
            lab(fig[[1]]) <- paste0("protein docking of ", name)
            return(fig)
          })
        figs <- unlist(figs, recursive = F)
        plots <- c(plots, list(top = figs))
      }
      x@plots[[ 4 ]] <- plots
      t.data <- .set_lab(t.data, sig(x), "Overview of protein docking results data.")
      x@tables[[ 4 ]] <- namel(t.data)
    } else {
      message("No score information, please pull results first.")
    }
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
    history <- .get_finished.cluspro(link)
    if (!is.null(x$groups)) {
      group <- unlist(x$groups)
      group <- group[ !group %in% x$excludes ]
      finished <- dplyr::filter(history, Name %in% !!group)
      message("Current finished: ", nrow(finished), "/", length(group))
      if (is.null(ref)) {
        ref <- finished$Name
      }
    }
    if (length(ref)) {
      finished <- dplyr::filter(history, Name %in% !!ref)
    } else {
      finished <- data.frame()
    }
    if (!is.null(x$tmp_scores)) {
      finishedNotGet <- dplyr::filter(finished, !Name %in% !!unique(x$tmp_scores$Name))
    } else {
      finishedNotGet <- finished
    }
    if (nrow(finishedNotGet)) {
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
      res <- pbapply::pbapply(finishedNotGet, 1, simplify = F,
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
          Sys.sleep(1)
          score <- fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
          namel(score, file_model)
        })
      names(res) <- finishedNotGet$Name
      scores <- lapply(res, function(x) x$score)
      scores <- frbind(scores, idcol = "Name")
      scores <- dplyr::rename_all(scores, make.names)
      models <- vapply(res, function(x) x$file_model, character(1))
      if (!is.null(x$tmp_scores)) {
        x$tmp_scores <- dplyr::bind_rows(x$tmp_scores, scores)
        x$tmp_models <- c(x$tmp_models, models)
      } else {
        x$tmp_scores <- scores
        x$tmp_models <- models
      }
    } else {
      if (length(ref)) {
        message("Job maybe not finished: ", paste0(ref, collapse = ", "))
      } else {
        message("No job has been finished.")
      }
    }
    x@step <- 3L
    return(x)
  })

setMethod("upload", signature = c(x = "job_cluspro"),
  function(x, whichGroup = NULL, maxAllow = 15){
    message("Collate current status and upload jobs.")
    if (is.null(whichGroup)) {
      notFinish <- which(!x$groups_status)
      if (!length(notFinish)) {
        message("All job has been upload.")
        return(x)
      }
      x$whichGroup <- whichGroup <- notFinish[1]
    } else {
      x$whichGroup <- whichGroup
    }
    message("Arrange for uploading: Group ", whichGroup, " of ", length(x$groups))
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
    ## 
    link <- x$link
    group <- rawGroup <- x$groups[[ whichGroup ]]
    ## check if got results
    url <- "https://cluspro.bu.edu/results.php"
    link$navigate(url)
    finished <- .get_finished.cluspro(link)
    finished <- dplyr::filter(finished, Name %in% !!group)
    message("Has got results: ", nrow(finished), "/", length(group), "/", length(rawGroup))
    group <- group[ !group %in% finished$Name ]
    ## check now running
    url <- "https://cluspro.bu.edu/queue.php"
    fun_get <- function(url) {
      fun_format <- function(x) {
        colnames(x) <- unlist(x[1, ])
        x[-1, ]
      }
      fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
    }
    queues <- dplyr::filter(fun_get(url), User == x$usr)
    if (nrow(queues) >= maxAllow) {
      stop("Too many jobs in server: ", nrow(queues))
    } else {
      message("Job Server in server: ", nrow(queues))
    }
    runs <- group[ group %in% queues$Name ]
    message("The group in queue: ", length(runs), "/", length(group), "/", length(rawGroup))
    group <- group[ !group %in% queues$Name ]
    if (length(group) + nrow(queues) > maxAllow) {
      stop("Upload will cause too many jobs: ", nrow(queues), " + ", length(group), " > ", maxAllow)
    } else {
      message("Enough space for further upload.")
    }
    ####################################################
    ####################################################
    ## upload
    if (length(group)) {
      message("Start upload.")
      url <- "https://cluspro.bu.edu/home.php"
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
          ## ligand
          ele <- link$findElement("xpath", "//div//span//input[@id='ligpdb']")
          ele$sendKeysToElement(list(layout[[ "from" ]]))
          ## rec
          ele <- link$findElement("xpath", "//div//span//input[@id='recpdb']")
          ele$sendKeysToElement(list(layout[[ "to" ]]))
          ## submit
          ele <- link$findElement("xpath", "//div//input[@type='submit']")
          ele$clickElement()
          message("JOB has successfully upload: ", i)
          Sys.sleep(3)
        })
      x$excludes <- unique(c(x$excludes, unlist(excludes)))
    }
    x$groups_status[ whichGroup ] <- T
    return(x)
  })


vis_pdb.cluspro <- function(pdb, label.a, label.b,
  save = paste0(get_savedir("figs"), "/", make.names(label.a), "_with_", make.names(label.b), ".png"),
  command = pg("pymol"))
{
  # label.a, rec
  # label.b, lig
  data <- ftibble(pdb, fill = T)
  sig <- which(data$V1 == "HEADER")
  coa <- dplyr::slice(data, sig[1]:(sig[2] - 1))
  cob <- dplyr::slice(data, sig[2]:nrow(data))
  fun_pos <- function(x) {
    c(max(x$V7, na.rm = T), median(x$V8, na.rm = T), median(x$V9, na.rm = T))
  }
  pos.a <- fun_pos(coa)
  pos.b <- fun_pos(cob)
  expr <- .expr_pretty_pdb(label.a, pos.a, label.b, pos.b,
    a = "rec.pdb", b = "lig.000.00.pdb"
  )
  vis_pdb(pdb, expr, save)
  return(save)
}

.get_finished.cluspro <- function(link) {
  fun_format <- function(x) {
    colnames(x) <- unlist(x[1, ])
    x[-1, ]
  }
  canGet <- T
  res <- list()
  n <- 0L
  while (canGet) {
    n <- n + 1L
    res[[ n ]] <- fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
    ele <- try(link$findElement("xpath", "//a[text()='next ->']"), T)
    if (inherits(ele, "try-error")) {
      canGet <- F
    } else {
      ele$clickElement()
      Sys.sleep(1)
    }
  }
  frbind(res)
}
