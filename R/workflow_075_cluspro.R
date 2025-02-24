# ==========================================================================
# workflow of cluspro
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_cluspro <- setClass("job_cluspro", 
  contains = c("job"),
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
    symbols <- unlist(.layout, use.names = FALSE)
  } else {
    symbols <- rm.no(symbols)
    .layout <- data.frame(t(combn(symbols, 2)))
    names(.layout) <- c("from", "to")
  }
  .job_cluspro(object = list(hgnc_symbols = rm.no(symbols)), params = list(.layout = .layout))
}

setMethod("step0", signature = c(x = "job_cluspro"),
  function(x){
    step_message("Prepare your data with function `job_cluspro`.")
  })

setMethod("step1", signature = c(x = "job_cluspro"),
  function(x, order = TRUE, each_target = 1, custom_pdbs = NULL){
    step_message("Found PDB.")
    x <- .find_proper_pdb(x, order, each_target, custom_pdbs)
    if (TRUE) {
      fun_set_layout <- function(anno) {
        layouts <- x$.layout
        layouts <- map(layouts, "from", anno, "hgnc_symbol", "pdb", rename = FALSE)
        layouts <- map(layouts, "to", anno, "hgnc_symbol", "pdb", rename = FALSE)
        layouts <- apply(layouts, 1, tolower, simplify = FALSE)
        names(layouts) <- paste0(x$.layout$from, "_", x$.layout$to)
        layouts
      }
      anno <- dplyr::distinct(x$targets_annotation, hgnc_symbol, .keep_all = TRUE)
      x$layouts <- fun_set_layout(anno)
    }
    x$dock_layout <- x$layouts
    res <- .get_pdb_files(x)
    x <- res$x
    if (any(is.na(unlist(x$dock_layout)))) {
      x$layouts <- fun_set_layout(
        data.frame(hgnc_symbol = names(x$used_pdbs), pdb = unname(x$used_pdbs))
      )
    }
    keep <- vapply(x$layouts, 
      function(x) {
        all(!is.na(x))
      }, logical(1))
    x$layouts <- x$layouts[ keep ]
    t.docking_layouts <- tibble::tibble(
      name = names(x$layouts), 
      from_id = toupper(vapply(x$layouts, function(x) x[[ "from" ]], character(1))),
      to_id = toupper(vapply(x$layouts, function(x) x[[ "to" ]], character(1)))
    )
    t.docking_layouts <- setLegend(t.docking_layouts, "为蛋白质对接所使用的 PDB 来源附表 (如果是 PDB ID，则 pdb 文件直接来源于 `PDB` 数据库；如果是 `UniProtKB-Swiss-Prot` ID，则根据该 ID，从 `AlphaFold` 获取预测的结构文件。")
    x <- tablesAdd(x, t.docking_layouts)
    x$pdb.files <- res$pdb.files
    x$pdb_from <- x$used_pdbs
    return(x)
  })

setMethod("step2", signature = c(x = "job_cluspro"),
  function(x, port = 9999, ..., usr = getOption("cluspro_user"),
    psw = getOption("cluspro_password"))
  {
    step_message("Login.")
    x <- login(x, port, ..., usr = usr, psw = psw)
    x <- methodAdd(x, "登录 `ClusPro` ({cite_show('TheClusproWebKozako2017')}) (<https://cluspro.bu.edu/home.php>)，上传需要对接的蛋白结构 (PDB) 。")
    x <- snapAdd(x, "将 PDB 上传至 `ClusPro` 进行对接。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_cluspro"),
  function(x, n = 10){
    step_message("Prepare groups for upload.")
    layout <- x$layouts
    groups <- grouping_vec2list(names(layout), n, TRUE)
    groups_status <- rep(FALSE, length(groups))
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
      data <- tidyr::separate(data, "Name", c("pro1", "pro2"), "_", remove = FALSE)
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
      p.score <- setLegend(p.score, "为每组对接结果的最小能量柱状图 (请参考 <https://cluspro.bu.edu/help.php>)。")
      plots <- namel(p.score)
      if (!is.null(top)) {
        data <- head(data, top)
        n <- 0L
        figs <- pbapply::pbapply(data, 1, simplify = FALSE,
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
        figs <- unlist(figs, recursive = FALSE)
        figs <- setLegend(figs, glue::glue("以 `pymol` 将蛋白质 ({names(figs)}) 对接结果可视化。"))
        plots <- c(plots, list(top = figs))
      }
      x@plots[[ 4 ]] <- plots
      t.data <- .set_lab(t.data, sig(x), "Overview of protein docking results data.")
      x@tables[[ 4 ]] <- namel(t.data)
    } else {
      stop("No score information, please pull results first.")
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
      dir.create(dir, FALSE)
      fun_format <- function(x) {
        cluster <- 0L
        members <- 0L
        names <- colnames(x)
        x <- apply(x, 1, simplify = TRUE,
          function(v) {
            if (grpl(v[[ "Cluster" ]], "^[0-9]+$")) {
              cluster <<- v[[ "Cluster" ]]
              members <<- v[[ "Members" ]]
              type <- TRUE
            } else {
              type <- FALSE
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
      res <- pbapply::pbapply(finishedNotGet, 1, simplify = FALSE,
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
          dir.create(subdir, FALSE)
          exist <- FALSE
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
    check <- try(x@params$link$getCurrentUrl(), TRUE)
    url <- "https://cluspro.bu.edu/home.php"
    if (inherits(check, "try-error")) {
      stop("Please check the web server is login.")
    } else if (check[[1]] != url) {
      x$link$navigate(url)
      Sys.sleep(1)
      check <- try(x@params$link$getCurrentUrl(), TRUE)
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
    link$navigate(url)
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
      fun_need_upload_pdbFile <- function(id) {
        id_needFile <- x$pdb_notGot_uniprot$Entry
        if (!is.null(id_needFile)) {
          any(tolower(id_needFile) == id)
        } else {
          FALSE
        }
      }
      if (any(grpl(names(x$pdb.files), "[A-Z]"))) {
        message('any(grpl(names(x$pdb.files), "[A-Z]")), convert to lower case.')
        names(x$pdb.files) <- tolower(names(x$pdb.files) )
      }
      pdb.files <- as.list(x$pdb.files)
      fun_get_pdb_file <- function(id) {
        file <- pdb.files[[ id ]]
        if (is.null(file)) {
          stop('is.null(file), do not known the filepath.')
        }
        normalizePath(file)
      }
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
          if (fun_need_upload_pdbFile(layout[[ "from" ]])) {
            ele_switch <- link$findElement("xpath", "//div//span[@id='showligfile']")
            ele_switch$clickElement()
            Sys.sleep(.5)
            ele <- link$findElement("xpath", "//div//input[@id='lig' and @type='file']")
            file <- fun_get_pdb_file(layout[[ "from" ]])
            ele$sendKeysToElement(list(file))
          } else {
            ele <- link$findElement("xpath", "//div//span//input[@id='ligpdb']")
            ele$sendKeysToElement(list(layout[[ "from" ]]))
          }
          ## rec
          if (fun_need_upload_pdbFile(layout[[ "to" ]])) {
            ele_switch <- link$findElement("xpath", "//div//span[@id='showrecfile']")
            ele_switch$clickElement()
            Sys.sleep(.5)
            ele <- link$findElement("xpath", "//div//input[@id='rec' and @type='file']")
            file <- fun_get_pdb_file(layout[[ "to" ]])
            ele$sendKeysToElement(list(file))
          } else {
            ele <- link$findElement("xpath", "//div//span//input[@id='recpdb']")
            ele$sendKeysToElement(list(layout[[ "to" ]]))
          }
          ## submit
          ele <- link$findElement("xpath", "//div//input[@type='submit']")
          ele$clickElement()
          message("JOB has successfully upload: ", i)
          Sys.sleep(3)
        })
      x$excludes <- unique(c(x$excludes, unlist(excludes)))
    }
    x$groups_status[ whichGroup ] <- TRUE
    return(x)
  })


vis_pdb.cluspro <- function(pdb, label.a, label.b,
  save = paste0(get_savedir("figs"), "/", make.names(label.a), "_with_", make.names(label.b), ".png"),
  command = pg("pymol"))
{
  # label.a, rec
  # label.b, lig
  data <- ftibble(pdb, fill = TRUE)
  sig <- which(data$V1 == "HEADER")
  coa <- dplyr::slice(data, sig[1]:(sig[2] - 1))
  cob <- dplyr::slice(data, sig[2]:nrow(data))
  fun_pos <- function(x) {
    c(max(x$V7, na.rm = TRUE), median(x$V8, na.rm = TRUE), median(x$V9, na.rm = TRUE))
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
  canGet <- TRUE
  res <- list()
  n <- 0L
  while (canGet) {
    n <- n + 1L
    res[[ n ]] <- fun_format(get_table.html(link$getPageSource()[[1]])[[1]])
    ele <- try(link$findElement("xpath", "//a[text()='next ->']"), TRUE)
    if (inherits(ele, "try-error")) {
      canGet <- FALSE
    } else {
      ele$clickElement()
      Sys.sleep(1)
    }
  }
  frbind(res)
}

.expr_pretty_pdb <- function(label.a, pos.a, label.b, pos.b,
  a = "chain A", b = "chain B")
{
  pos <- function(x) paste0("[", paste0(x, collapse = ","), "]")
  paste0(
    " -d \"",
    "color pink, ", a, "; ",
    "color yellow, ", b, "; ",
    "pseudoatom pA, pos = ", pos(pos.a), "; ",
    "pseudoatom pB, pos = ", pos(pos.b), "; ",
    "set label_size, 25",
    "\" ",
    " -d \"", "label pA, '", label.a, "'", "\" ",
    " -d \"", "label pB, '", label.b, "'", "\" "
    )
}
