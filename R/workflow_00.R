# ==========================================================================
# for wrap various of analysis workflow contains many steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("virtual_job", "VIRTUAL")

.meth <- setClass("meth", 
  contains = c("list"),
  representation = representation(params = "list"),
  prototype = NULL)

#' @exportClass job
.job <- setClass("job", 
  contains = c("virtual_job"),
  representation = representation(
    object = "ANY",
    params = "ANY",
    plots = "ANY",
    tables = "ANY",
    others = "ANY",
    step = "integer",
    info = "ANY",
    pg = "character",
    cite = "character",
    method = "character",
    meth = "list",
    sig = "character",
    tag = "character",
    analysis = "character"
    ),
  prototype = prototype(step = 0L,
    tag = character(1),
    analysis = character(1),
    meth = .meth(),
    params = list(wd = ".", remote = "remote"),
    others = list()
    ))

setGeneric("snap", 
  function(x, ...) standardGeneric("snap"))

setGeneric("meth", 
  function(x, ...) standardGeneric("meth"))
setMethod("meth", 
  signature = c(x = "ANY"),
  function(x){
    x@meth
  })

setGeneric("meth<-", 
  function(x, value) standardGeneric("meth<-"))
setReplaceMethod("meth", signature = c(x = "ANY"),
  function(x, value){
    initialize(x, meth = value)
  })

collate_details <- function(tag = "meth", envir = parent.frame(1)) {
  expr <- stringr::str_extract(unlist(knitr::knit_code$get()),
    paste0("(?<=^#' @", tag, " ).*"))
  expr <- expr[ !is.na(expr) ]
  sys <- Sys.info()
  expr <- c("## 数据分析平台", "",
    glue::glue("在 {paste(sys[[ 'sysname' ]], sys[[ 'nodename' ]], sys[[ 'machine' ]])} ({sys[[ 'release' ]]}) 上，使用 {R.Version()$version.string} (https://www.r-project.org/) 对数据统计分析与整合分析。"),
    "", expr)
  if (length(expr)) {
    text <- lapply(expr, glue::glue, .sep = "\n", .envir = envir)
    writeLines(unlist(text))
  }
}

get_meth <- function(x, setHeader = T) {
  header <- NULL
  if (is(x, "job")) {
    if (length(x@analysis)) {
      if (length(x@sig)) {
        sig <- glue::glue("(Dataset: {x@sig})")
      } else {
        sig <- ""
      }
      header <- paste("## ", x@analysis, sig)
    }
  }
  if (setHeader) {
    c("", header, "", meth(x))
  } else {
    c("", meth(x), "")
  }
}

setMethod("show", 
  signature = c(object = "meth"),
  function(object){
    lapply(object, writeLines)
    writeLines("\n")
  })

# setValidity("job", 
#   function(object){
#     class <- class(object)
#     options(method_name = class)
#     if (is(object, "job")) T else F
#   })

setClass("hclust")
## cutree, the tree higher, the sort number lower

.marker_list <- setClass("marker_list", 
  contains = c(),
  representation = representation(level = "numeric", marker = "list"),
  prototype = NULL)

setMethod("show", 
  signature = c(object = "marker_list"),
  function(object){
    .show(object)
  })

#' @exportMethod show
setMethod("show", 
  signature = c(object = "job"),
  function(object){
    message("A workflow of '", class(object), "' in step (done): ", object@step)
    message("Object size: ", obj.size(object))
    if (!is.null(object@params$set_remote))
      message("Remote: ", object@params$remote)
  })

setGeneric("params", 
  function(x, ...) standardGeneric("params"))
setGeneric("plots", 
  function(x, step, ...) standardGeneric("plots"))
setGeneric("tables", 
  function(x, step, ...) standardGeneric("tables"))
setGeneric("others", 
  function(x, ...) standardGeneric("others"))
setGeneric("object", 
  function(x, ...) standardGeneric("object"))

setGeneric("pg", 
  function(x, ...) standardGeneric("pg"))
setGeneric("relative", 
  function(x, ...) standardGeneric("relative"))
setGeneric("object<-", 
  function(x, value) standardGeneric("object<-"))
setGeneric("plots<-", 
  function(x, value) standardGeneric("plots<-"))
setGeneric("tables<-", 
  function(x, value) standardGeneric("tables<-"))
setGeneric("params<-", 
  function(x, value) standardGeneric("params<-"))
setGeneric("others<-", 
  function(x, value) standardGeneric("others<-"))
setGeneric("ids", 
  function(x, ...) standardGeneric("ids"))
setGeneric("sig", 
  function(x, ...) standardGeneric("sig"))
setGeneric("sig<-", 
  function(x, value) standardGeneric("sig<-"))
setGeneric("lab", 
  function(x, ...) standardGeneric("lab"))
setGeneric("lab<-", 
  function(x, value) standardGeneric("lab<-"))
setGeneric("pattern", 
  function(x, ...) standardGeneric("pattern"))
setGeneric("ping", 
  function(x, ...) standardGeneric("ping"))
setGeneric("login", 
  function(x, ...) standardGeneric("login"))
setGeneric("ref", 
  function(x, ...) standardGeneric("ref"))

.character_ref <- setClass("character_ref", 
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

setMethod("ref", signature = c(x = "character_ref"),
  function(x, add_internal_job = T, ...){
    ref <- strx(x, "[A-Z][A-Za-z]{10,}[0-9]{4}")
    message("Extracted reference: ", ref)
    if (add_internal_job) {
      job <- .job_publish(cite = paste0("[@", ref, "]"))
      .add_internal_job(job)
    }
    return(x)
  })

setMethod("ref", signature = c(x = "character"),
  function(x, ...){
    codes_list <- knitr::knit_code$get()
    if (!length(codes_list)) {
      stop("!length(codes_list), no in rendering circumstance?")
    }
    codes <- codes_list[[x]] 
    if (is.null(codes)) {
      stop('is.null(codes), can not found ref target, please check the character.')
    }
    objName <- stringr::str_extract(codes, "(?<=autor\\().*(?=\\)$)")
    objName <- objName[ !is.na(objName) ]
    if (!length(objName)) {
      stop('!length(objName), can not match object name in "autor()".')
    } else if (length(objName) > 1) {
      stop("Too many object name matched in 'autor()'.")
    }
    obj <- eval(parse(text = objName), parent.frame(2))
    if (is(obj, "df")) {
      glue::glue("Tab. \\@ref(tab:{x})")
    } else if (is(obj, "can_not_be_draw") || is(obj, "can_be_draw") || is(obj, "fig")) {
      glue::glue("Fig. \\@ref(fig:{x})")
    }
  })

setMethod("lab", signature = c(x = "ANY"),
  function(x, ...){
    attr(x, ".LABEL")
  })

setReplaceMethod("lab", signature = c(x = "ANY", value = "character"),
  function(x, value){
    if (is.null(x)) {
      return(x)
    }
    attr(x, ".LABEL") <- value
    return(x)
  })

.lab_out <- function(x) {
  if (requireNamespace("nvimcom")) {
    x <- lab(x)
    if (is.null(x)) {
      x <- getOption("autor_unnamed_number", 1)
      options(autor_unnamed_number = x + 1)
      label <- paste0("#| Unnamed-", x)
    } else {
      label <- as_chunk_label(x)
      label <- paste0("#| ", label)
    }
    .C("nvimcom_msg_to_nvim", glue::glue('append(line(".") - 1, "{label}")'), PACKAGE = "nvimcom")
  }
}

as_chunk_label <- function(x) {
  Hmisc::capitalize(gs(gs(x, "(?<![a-zA-Z])ids[ \\-]?", "", perl = T), " |_|\\.", "-"))
}

.set_lab <- function(x, sig, group = NULL, body = NULL, suffix = NULL) {
  if (identical(class(x), "list")) {
    n <- 0L
    x <- lapply(x,
      function(obj) {
        n <<- n + 1L
        group <- group[n]
        if (is.na(group)) {
          stop("is.na(group)")
        }
        if (is(obj, "can_be_draw")) {
          obj <- wrap(obj)
        }
        lab(obj) <- gs(gs(paste(sig, group, body, suffix), "[ ]+", " "), "^[ \\-]|[ \\-]$", "")
        return(obj)
      })
  } else {
    if (is(x, "can_be_draw")) {
      x <- wrap(x)
    }
    lab(x) <- gs(gs(paste(sig, group, body, suffix), "[ ]+", " "), "^[ \\-]|[ \\-]$", "")
  }
  return(x)
}

setMethod("sig", signature = c(x = "virtual_job"),
  function(x){
    x@sig
  })
setReplaceMethod("sig", signature = c(x = "virtual_job"),
  function(x, value){
    initialize(x, sig = value)
  })

#' @exportMethod pg
setMethod("pg", signature = c(x = "job"),
  function(x,
    recode.remote = getOption("pg_remote_recode", pg_remote_recode()),
    recode.local = getOption("pg_local_recode", pg_local_recode()))
  {
    pg(x@pg, is.remote(x), recode.remote, recode.local)
  })

setMethod("pg", signature = c(x = "character"),
  function(x, is.remote = F,
    recode.remote = getOption("pg_remote_recode", pg_remote_recode()),
    recode.local = getOption("pg_local_recode", pg_local_recode()))
  {
    cli::cli_alert_info(x)
    if (is.remote) {
      recode <- recode.remote[[ x ]]
    } else {
      recode <- recode.local[[ x ]]
    }
    if (is.null(recode)) {
      x
    } else {
      recode
    }
  })

pg_local_recode <- function() {
  conda <- getOption("conda", "~/miniconda3")
  lst <- list(
    python = "{conda}/bin/python3",
    conda = "{conda}/bin/conda",
    conda_env = "{conda}/envs",
    qiime = "{conda}/bin/conda run -n qiime2 qiime",
    musitePython = "{conda}/bin/conda run -n musite python3",
    musitePTM = "~/MusiteDeep_web/MusiteDeep/predict_multi_batch.py",
    musitePTM2S = "~/MusiteDeep_web/PTM2S/ptm2Structure.py",
    dl = normalizePath("~/D-GCAN/DGCAN"),
    dl_dataset = normalizePath("~/D-GCAN/dataset"),
    dl_model = normalizePath("~/D-GCAN/DGCAN/model"),
    scfeaPython = "{conda}/bin/conda run -n scFEA python",
    scfea = "~/scFEA/src/scFEA.py",
    scfea_db = "~/scFEA/data",
    musiteModel = normalizePath("~/MusiteDeep_web/MusiteDeep/models"),
    mk_prepare_ligand.py = "mk_prepare_ligand.py",
    prepare_gpf.py = "prepare_gpf.py",
    autogrid4 = "autogrid4",
    scsa = "python3 ~/SCSA/SCSA.py",
    scsa_db = "~/SCSA/whole_v2.db",
    pymol = "pymol",
    # sirius = .prefix("sirius/bin/sirius", "op"),
    obgen = "obgen"
  )
  envir <- environment()
  lapply(lst, glue::glue, .envir = envir)
}

pg_remote_recode <- function() {
  conda_remote <- getOption("conda_remote", "~/miniconda3")
  lst <- list(
    qiime = "{conda_remote}/bin/conda run -n qiime2 qiime",
    fastp = "{conda_remote}/bin/conda run -n base fastp",
    bcftools = "{conda_remote}/bin/conda run -n base bcftools",
    elprep = "{conda_remote}/bin/conda run -n base elprep",
    # biobakery_workflows = "{conda_remote}/bin/conda run -n biobakery biobakery_workflows",
    bowtie2 = "{conda_remote}/bin/conda run -n base bowtie2",
    samtools = "{conda_remote}/bin/conda run -n base samtools",
    metaphlan = "{conda_remote}/bin/conda run -n mpa metaphlan",
    merge_metaphlan_tables.py = "{conda_remote}/bin/conda run -n mpa merge_metaphlan_tables.py",
    sirius = "~/operation/sirius/bin/sirius",
    scfeaPython = "{conda_remote}/bin/conda run -n scFEA python",
    scfea = "~/scFEA/src/scFEA.py",
    scfea_db = "~/scFEA/data"
  )
  envir <- environment()
  lapply(lst, glue::glue, .envir = envir)
}

#' @exportMethod object
setMethod("object", signature = c(x = "job"),
  function(x){
    x@object
  })
setReplaceMethod("object", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })

#' @exportMethod plots
setMethod("plots", signature = c(x = "job"),
  function(x){
    x@plots[[ x@step ]]
  })

#' @exportMethod plots
setMethod("plots", signature = c(x = "job", step = "numeric"),
  function(x, step){
    x@plots[[ step ]]
  })

setReplaceMethod("plots", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })

plotsAdd <- function(x, ...) {
  jobSlotAdd(x, "plots", ...)
}

tablesAdd <- function(x, ...) {
  jobSlotAdd(x, "tables", ...)
}

jobSlotAdd <- function(x, name, ...) {
  calls <- getNamesOfCall(substitute(list(...)))
  objs <- list(...)
  names(objs) <- calls
  labs <- formatNames(calls)
  db <- if (length(slot(x, name)) >= x@step) {
    slot(x, name)[[ x@step ]]
  } else {
    list()
  }
  slot(x, name)[[ x@step ]] <- c(db, .set_lab_batch(objs, labs, sig(x)))
  return(x)
}

methodAdd <- function(x, glueString, env = parent.frame(1)) {
	meth(x)[[ paste0("step", x@step) ]] <- glue::glue(glueString, .envir = env)
  return(x)
}

.set_lab_batch <- function(objs, labs, sig) {
  mapply(objs, labs, SIMPLIFY = F,
    FUN = function(obj, lab) {
      if (is.null(lab(obj))) {
        obj <- .set_lab(obj, sig, lab)
      }
      return(obj)
    })
}

formatNames <- function(names) {
  names <- gsub("\\.|_", " ", sub("^[a-z]\\.", "", names))
  gsub("([a-z0-9])([A-Z])", "\\1 \\2", names)
}

#' @exportMethod tables
setMethod("tables", signature = c(x = "job"),
  function(x){
    x@tables[[ x@step ]]
  })

#' @exportMethod tables
setMethod("tables", signature = c(x = "job", step = "double"),
  function(x, step){
    x@tables[[ step ]]
  })

#' @exportMethod tables
setMethod("tables", signature = c(x = "job", step = "double"),
  function(x, step){
    x@tables[[ step ]]
  })
setReplaceMethod("tables", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })


#' @exportMethod params
setMethod("params", signature = c(x = "job"),
  function(x){
    x@params
  })
setReplaceMethod("params", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })


#' @exportMethod others
setMethod("others", signature = c(x = "job"),
  function(x){
    x@others[[ x@step ]]
  })
setReplaceMethod("others", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })

# ==========================================================================
# remote or local 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setGeneric("set_remote", 
  function(x, ...) {
    x <- set_remote.default(x, tmpdir = getOption("remote_tmpdir", NULL),
      map_local = paste0(gs(class(x), "^job_", ""), "_local"),
      remote = "remote"
    )
    standardGeneric("set_remote")
  })

set_remote.default <- function(x, tmpdir, map_local, remote) {
  if (is.null(x$set_remote))
    x$set_remote <- T
  if (is.null(x$remote))
    x$remote <- remote
  if (is.null(x$map_local))
    x$map_local <- map_local
  if (is.null(x$tmpdir))
    x$tmpdir <- tmpdir
  return(x)
}

setMethod("set_remote", signature = c(x = "job"),
  function(x, wd){
    x$wd <- wd
    return(x)
  })

# ==========================================================================
# methods of step series
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

available_signatures <- function(f, exclude = c("ANY", "missing", "job")) {
  signatures <- findMethodSignatures(f)[, 1]
  signatures <- signatures[ !signatures %in% exclude ]
  signatures
}

setGeneric("step0", 
  function(x, ...) {
    if (!missing(x)) {
      if (!is(x, "numeric") & !is(x, "character")) {
        message("This step do nothing but description.")
        message(crayon::silver("Use `step0()` to show all available workflow.\n"))
      }
    }
    standardGeneric("step0")
  })

setMethod("step0", signature = c(x = "missing"),
  function(x){
    message(crayon::blue("All available workflow for class:\n"))
    signatures <- available_signatures("step1")
    mess <- paste0(paste0(1:length(signatures), ": ", signatures), collapse = "\n")
    step_message(mess)
    message(crayon::silver("Example of show description: use `step0(1)`."))
  })

setMethod("step0", signature = c(x = "numeric"),
  function(x){
    signatures <- available_signatures("step1")
    fun <- get_fun(paste0(".", signatures[ x ]))
    step0(fun())
  })

setMethod("step0", signature = c(x = "character"),
  function(x){
    signatures <- available_signatures("step1")
    sig <- match.arg(x, sub("job_", "", signatures))
    fun <- get_fun(paste0(".job_", sig))
    step0(fun())
    options(method_name = paste0("job_", sig))
  })

set_sig <- function(x) {
  name <- substitute(x)
  sig <- toupper(gs(gs(gs(name, "^[a-zA-Z0-9]+", ""), "\\.", " "), "^[ ]*", ""))
  attr(sig, "name") <- as.character(name)
  sig(x) <- sig
  return(x)
}

setGeneric("step1",
  function(x, ...) {
    # if (identical(sig(x), character(0))) {
    if (identical(parent.frame(1), .GlobalEnv)) {
      name <- substitute(x)
      sig <- toupper(gs(gs(gs(name, "^[a-zA-Z0-9]+", ""), "\\.", " "), "^[ ]*", ""))
      attr(sig, "name") <- as.character(name)
      sig(x) <- sig  
    }
    # }
    x$seed <- sample(1:100, 1)
    x <- checkAddStep(x, 1L)
    x <- standardGeneric("step1")
    x <- stepPostModify(x, 1)
    if (identical(parent.frame(1), .GlobalEnv)) {
      job_append_heading(x)
    }
    x
  })

job_append_heading <- function (x, mutate = T) {
  if (getOption("job_appending", F)) {
    if (is(x, "job")) {
      heading <- x@analysis
      if (length(x@sig)) {
        heading <- paste0(heading, " (", x@sig, ")")
      }
      .append_heading(heading, mutate)
    }
  }
}

.append_heading <- function(x, mutate = F) {
  args <- paste0("CheckAndAppendHeading('", x, "'", if (mutate) ", 1" else ", 0", ")")
  if (requireNamespace("nvimcom")) {
    if (nchar(x)) {
      .C("nvimcom_msg_to_nvim", args, PACKAGE = "nvimcom")
    }
  }
}

.append_class <- function(expr, short = T, envir = .GlobalEnv) {
  if (requireNamespace("nvimcom")) {
    class <- class(expr)
    if (length(class) > 1) {
      len <- nchar(class)
      class <- class[ which(len == min(len))[1] ]
    }
    args <- glue::glue("InsertText_I('{class}', -1)")
    message("Send to vim: ", args)
    .C("nvimcom_msg_to_nvim", args, PACKAGE = "nvimcom")
  }
}

setGeneric("step2",
  function(x, ...) {
    x <- checkAddStep(x, 2L)
    x <- standardGeneric("step2")
    stepPostModify(x, 2)
  })

setGeneric("step3",
  function(x, ...) {
    x <- checkAddStep(x, 3L)
    x <- standardGeneric("step3")
    stepPostModify(x, 3)
  })

setGeneric("step4",
  function(x, ...) {
    x <- checkAddStep(x, 4L)
    x <- standardGeneric("step4")
    stepPostModify(x, 4)
  })

setGeneric("step5",
  function(x, ...) {
    x <- checkAddStep(x, 5L)
    x <- standardGeneric("step5")
    stepPostModify(x, 5)
  })

setGeneric("step6",
  function(x, ...) {
    x <- checkAddStep(x, 6L)
    x <- standardGeneric("step6")
    stepPostModify(x, 6)
  })

setGeneric("step7",
  function(x, ...) {
    x <- checkAddStep(x, 7L)
    x <- standardGeneric("step7")
    stepPostModify(x, 7)
  })

setGeneric("step8",
  function(x, ...) {
    x <- checkAddStep(x, 8L)
    x <- standardGeneric("step8")
    stepPostModify(x, 8)
  })

setGeneric("step9",
  function(x, ...) {
    x <- checkAddStep(x, 9L)
    x <- standardGeneric("step9")
    stepPostModify(x, 9)
  })

setGeneric("step10",
  function(x, ...) {
    x <- checkAddStep(x, 10L)
    x <- standardGeneric("step10")
    stepPostModify(x, 10)
  })

setGeneric("step11",
  function(x, ...) {
    x <- checkAddStep(x, 11L)
    x <- standardGeneric("step11")
    stepPostModify(x, 11)
  })

setGeneric("step12",
  function(x, ...) {
    x <- checkAddStep(x, 12L)
    x <- standardGeneric("step12")
    stepPostModify(x, 12)
  })

stepPostModify <- function(x, n) {
  cli::cli_h1("Job finished & Start post modify")
  x <- convertPlots(x, n)
  xname <- attr(x@sig, "name")
  news <- c()
  if (length(x@plots) >= n) {
    if (!is.null(x@plots[[ n ]]))
      names(x@plots)[ n ] <- paste0("step", n)
    new_plots <- paste0(xname, "@plots$step", x@step, "$", names(x@plots[[ x@step ]]))
    cli::cli_alert_info("Plot news:")
    lapply(new_plots, cli::cli_alert_success)
    news <- c(news, new_plots)
  }
  if (length(x@tables) >= n) {
    if (!is.null(x@tables[[ n ]]))
      names(x@tables)[ n ] <- paste0("step", n)
    new_tables <- paste0(xname, "@tables$step", x@step, "$", names(x@tables[[ x@step ]]))
    cli::cli_alert_info("Table news:")
    lapply(new_tables, cli::cli_alert_success)
    news <- c(news, new_tables)
  }
  isNewParams <- !names(x@params) %in% x@others$.oldParams
  if (any(isNewParams)) {
    cli::cli_alert_info("Params news:")
    new_params <- paste0(xname, "@params$", names(x@params)[ isNewParams ])
    lapply(new_params, cli::cli_alert_success)
    news <- c(news, new_params)
  }
  if (!is.null(xname) && length(xname)) {
    assign(xname, x)
    writeJobSlotsAutoCompletion(xname, environment())
  }
  x
}

writeAllCompletion <- function(jobs = .get_job_list(extra = NULL)) {
  updAlls(jobs)
  lapply(jobs,
    function(x) {
      writeJobSlotsAutoCompletion(x)
    })
  invisible()
}

updAlls <- function(names = .get_job_list(extra = NULL), envir = .GlobalEnv) {
  lapply(names,
    function(x) {
      obj <- get(x, envir = envir)
      if (is(obj, "job")) {
        if (inherits(try(validObject(obj), T), "try-error")) {
          obj <- upd(obj)
          assign(x, obj, envir = envir)
        }
      }
    })
}

writeJobSlotsAutoCompletion <- function(x, envir = .GlobalEnv) {
  ## x is character(1)
  if (!is.character(x)) {
    stop("`x` must be character(1).")
  }
  dir <- getOption("SelectCompletion")
  if (is.null(dir)) {
    return()
  }
  dir.create(dir, F, recursive = T)
  obj <- try(get(x, envir = envir), T)
  if (inherits(obj, "try-error")) {
    message("May not run in .GlobalEnv.")
    return()
  }
  if (!is.null(dir) || !is(obj, "job")) {
    filename <- file.path(dir, paste0(class(obj), "__", x, ".R"))
    content <- treeObj(x, envir = envir)
    writeLines(content, filename)
    message("Setup completion: ", filename)
  }
}

convertPlots <- function(x, n) {
  if (length(x@plots) >= n) {
    x@plots[[ n ]] <- rapply2_asWrap(x@plots[[ n ]], "gg.obj")
  }
  gc()
  x
}

rapply2_asWrap <- function(x, class = "gg.obj") {
  fun <- function(x) {
    if (is.null(x))
      return()
    if (is(x, "gg.obj"))
      wrap(as_grob(x))
    else {
      wrap(x)
    }
  }
  if (is(x, class)) {
    name <- attr(x, ".name")
    x <- fun(x)
    cli::cli_alert_success(paste("Convert", paste0("`", name, "`"), "as wrap"))
    x
  } else if (is(x, "wrap")) {
    if (is(x@data, class)) {
      x@data <- as_grob(x@data)
    }
    x
  } else if (is(x, "can_not_be_draw")) {
    x
  } else if (is(x, "list")) {
    attrs <- attributes(x)
    names <- names(x)
    x <- lapply(1:length(x),
      function(n) {
        name <- names[[ n ]]
        obj <- x[[ n ]]
        if (!is.null(obj))
          attr(obj, ".name") <- name
        obj
      })
    x <- lapply(x, rapply2_asWrap, class = class)
    attributes(x) <- attrs
    x
  } else {
    x
  }
}

rapply2 <- function(x, fun, class, ...) {
  if (is(x, class)) {
    fun(x, ...)
  } else if (is(x, "list")) {
    lapply(x, rapply2, fun = fun, class = class, ...)
  } else {
    x
  }
}

checkAddStep <- function(x, n, clear_tables_plots = T) {
  STEP <- paste0("step", n)
  if (getOption("step_check", T)) {
    if (x@step >= n)
      stop("x@step >= ", n, ". This step has been executed.")
    if (x@step < n - 1)
      stop("x@step <", n - 1, ". The previous step has not been executed.")
    if (STEP == "step1") {
      if (!is.null(x@info)) {
        message(crayon::yellow("Workflow information:"))
        step_message(x@info, show_end = NULL)
      }
    }
  }
  if (clear_tables_plots) {
    x@plots[[ n ]] <- NULL
    x@tables[[ n ]] <- NULL
  }
  message(crayon::green("Running", STEP), ":")
  x@step <- n
  x@others$.oldParams <- names(x@params)
  x
}

step_message <- function(..., show_end = "Job Status") {
  str <- paste0(unlist(list(...)), collapse = "")
  str <- gs(str, "red\\{\\{(.*?)\\}\\}", "\033[31m\\1\033[39m")
  str <- gs(str, "grey\\{\\{(.*?)\\}\\}", "\033[90m\\1\033[39m")
  str <- gs(str, "yellow\\{\\{(.*?)\\}\\}", "\033[33m\\1\033[39m")
  str <- gs(str, "blue\\{\\{(.*?)\\}\\}", "\033[34m\\1\033[39m")
  textSh(str, pre_trunc = F)
  if (!is.null(show_end))
    cli::cli_h1(show_end)
}

.argEnv <- new.env()

E <- function(expr, text = NULL, internal = F) {
  e(expr, text, internal, n = 1)
}

e <- function(expr, text = NULL, internal = T, n = NULL) {
  if (is.null(n)) {
    expr <- substitute(expr)
    n <- 1
  } else {
    expr <- substitute(expr, parent.frame(n))
    n <- 2
  }
  if (internal)
    names <- stringr::str_extract(deparse(expr), "[a-zA-Z0-9_.]*:::?[^\\(]*")
  else {
    names <- gs(deparse(expr), "^.*cdRun\\(([^ ]+ [^ ]+ [^ ^\"]+).*$", "\\1")
    names <- gs(names, "\"|,", "")
  }
  names <- names[ !is.na(names) ]
  if (!length(names))
    name <- "FUNCTION"
  else name <- names[1]
  if (is.null(text))
    cli::cli_alert_info(name)
  else
    cli::cli_alert_info(paste0(name, " (", text, ")" ))
  suppressMessages(eval(expr, envir = parent.frame(n)))
}

parallel <- function(x, fun, workers = 3) {
  # if (isNamespaceLoaded("future")) {
    # e(pkgload::unload("future"))
  # }
  options(future.globals.maxSize = Inf)
  e(future::plan(future::multicore, workers = workers))
  x <- fun(x)
  x
}

parallel_collect <- function() {
  future::value(future::future(Sys.getpid()))
}

show_nonstandardGenericFunction <- selectMethod(
  "show", "nonstandardGenericFunction"
)

show_standardGeneric <- selectMethod(
  "show", "standardGeneric"
)

.show_method <- function(object, default, filter_by_options){
  if (!hasMethods(object, where = topenv())) {
    return(default(object))
  }
  methods <- findMethods(object)
  len <- length(methods@names)
  if (filter_by_options) {
    if (!is.null(ms <- getOption("method_name"))) {
      nums <- which(methods@names %in% ms)
      if (!length(nums))
        return(cli::cli_h1(paste0("No methods for ",
              paste0(ms, collapse = ", "))))
    } else {
      nums <- 1:len
    }
  } else {
    nums <- 1:len
  }
  lapply(nums,
    function(n) {
      args <- "test"
      str <- deparse(methods@.Data[[n]])
      posLocal <- grep("\\.local <- function \\(", str)
      if (length(posLocal) > 0)
        str <- str[ posLocal[1]:length(str) ]
      else
        str[1] <- sub("^.*?(function)", "\\1", str[1])
      str <- paste0(str, collapse = " ")
      args <- match_pair(str)
      message(crayon::green(gsub("#", ", ", methods@names[[n]])), ":")
      message(strwrap(args, indent = 4))
    })
  cli::cli_h1("Methods parameters")
}

match_pair <- function(str, left = "(", right = ")") {
  str <- strsplit(str, "")[[1]]
  sig <- 0L
  pair <- 0L
  char <- c()
  for (i in str) {
    if (i == left) {
      pair <- 1L
      sig <- sig + 1L
    } else if (i == right) {
      sig <- sig - 1L
    }
    if (sig == 0L & pair)
      break
    else if (pair) {
      if (i == left & sig == 1L)
        next
    char <- c(char, i)
    }
  }
  paste0(char, collapse = "")
}

setMethod("show", 
  signature = c(object = "standardGeneric"),
  function(object){
    .show_method(object, show_standardGeneric, F)
  })

setMethod("show", 
  signature = c(object = "nonstandardGenericFunction"),
  function(object) {
    .show_method(object, show_nonstandardGenericFunction, T)
  }
)

setMethod("print", signature = c(x = "nonstandardGenericFunction"),
  function(x){
    show(x)
  })

setGeneric("clear", 
  function(x, ...) standardGeneric("clear"))

setMethod("clear", signature = c(x = "job"),
  function(x, save = T, suffix = NULL, name = substitute(x, parent.frame(1))){
    if (save)
      saveRDS(x, paste0(name, ".", x@step, suffix, ".rds"))
    object(x) <- NULL
    return(x)
  })

setGeneric("via_symbol", 
  function(x, ...) standardGeneric("via_symbol"))

setGeneric("vis",
  function(x, ...) standardGeneric("vis"))

setGeneric("focus", 
  function(x, ...) standardGeneric("focus"))

setGeneric("meta", 
  function(x, ...) standardGeneric("meta"))

setGeneric("anno", 
  function(x, ...) standardGeneric("anno"))

setGeneric("map", 
  function(x, ref, ...) standardGeneric("map"))

setGeneric("gname", 
  function(x, ...) standardGeneric("gname"))

setGeneric("tops", 
  function(x, ...) standardGeneric("tops"))

setGeneric("regroup", 
  function(x, ref, ...) standardGeneric("regroup"))

setGeneric("diff", 
  function(x, ...) standardGeneric("diff"))

setGeneric("infer", 
  function(x, ...) standardGeneric("infer"))

setGeneric("getsub", 
  function(x, ...) standardGeneric("getsub"))

setGeneric("active", 
  function(x, ...) standardGeneric("active"))

setGeneric("skel", 
  function(x, ...) standardGeneric("skel"))

setGeneric("res", 
  function(x, ref, ...) standardGeneric("res"))

setGeneric("intersect", 
  function(x, y, ...) standardGeneric("intersect"))

setGeneric("fill", 
  function(x, ...) standardGeneric("fill"))

setGeneric("upd", 
  function(x, ...) standardGeneric("upd"))

setGeneric("upload", 
  function(x, ...) standardGeneric("upload"))

setGeneric("pull", 
  function(x, ...) standardGeneric("pull"))

setGeneric("not", 
  function(x, ...) standardGeneric("not"))

setMethod("intersect", signature = c(x = "ANY", y = "ANY"),
  function(x, y){
    base::intersect(x, y)
  })

setMethod("not", signature = c(x = "job"),
  function(x){
    options(method_name = class(x))
    invisible(x)
  })

setGeneric("merge", 
  function(x, y, ...) standardGeneric("merge"))

setMethod("upd", signature = c(x = "ANY"),
  function(x){
    new <- new(class(x))
    new@object <- x@object
    new@params <- x@params
    new@plots <- x@plots
    new@tables <- x@tables
    new@step <- x@step
    if (!inherits(try(x@sig, T), "try-error")) {
      new@sig <- x@sig
    }
    new
  })

setMethod("fill", signature = c(x = "df"),
  function(x, ref, fill, lst){
    for (i in 1:length(lst)) {
      x[x[[ref]] == names(lst[i]), ][[ fill ]] <- lst[[ i ]]
    }
    return(x)
  })

ins <- function(..., lst = NULL) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  if (!is.null(lst) && length(list(...))) {
    lst <- c(lst, list(...))
  }
  x <- lst[[ 1 ]]
  for (i in lst[2:length(lst)]) {
    x <- intersect(x, i)
  }
  x
}

svr <- function(object) {
  saveRDS(object, paste0(substitute(object, parent.frame(1)), ".rds"))
}

ldr <- function(object) {
  x <- loadRDS(obj <- paste0(substitute(object, parent.frame(1)), ".rds"))
  assign(obj, x, envir = parent.frame(2))
}

set_palette <- function(x, values = color_set()) {
  x + scale_color_manual(values = values) +
    scale_fill_manual(values = values)
}

setGeneric("is.remote", 
  function(x, ...) standardGeneric("is.remote"))

setMethod("is.remote", signature = c(x = "job"),
  function(x){
    if (!is.null(x@params$set_remote)) {
      if (x$set_remote) T else F
    } else F
  })

setGeneric("is_workflow_object_exists", 
  function(object, ...) standardGeneric("is_workflow_object_exists"))

setMethod("[[", signature = c(x = "job"),
  function(x, i, ...){
    x@params[[ i ]]
  })

setMethod("[[<-", signature = c(x = "job"),
  function(x, i, ..., value){
    x@params[[ i ]] <- value
    return(x)
  })

setMethod("$", signature = c(x = "job"),
  function(x, name){
    x[[ name ]]
  })

setMethod("$<-", signature = c(x = "job"),
  function(x, name, value){
    x[[ name ]] <- value
    return(x)
  })

setMethod("map", signature = c(x = "df"),
  function(x, ref, y, y.ref, y.get, rename = T, col = NULL)
  {
    if (any(!c(y.ref, y.get) %in% colnames(y))) {
      stop("Can not found columns of '", y.ref, "' and '", y.get, "'")
    }
    if (is.null(col)) {
      x[[ ref ]] <- y[[ y.get ]][match(x[[ ref ]], y[[ y.ref ]])]
      if (rename) {
        colnames(x)[ colnames(x) == ref ] <- y.get
      }
    } else {
      x[[ col ]] <- y[[ y.get ]][match(x[[ ref ]], y[[ y.ref ]])]
    }
    x
  })

# setMethod("map", signature = c(x = "list"),
  # function(x, y, y.ref, y.get, ...)
  # {
  #   if (is(x, "df")) {
  #     return(callNextMethod(x, y, y.ref, y.get, ...))
  #   }
  #   lapply(x,
  #     function(x) {
  #       y[[ y.get ]][match(x, y[[ y.ref ]])]
  #     })
#   })

setMethod("gname", signature = c(x = "character"),
  function(x){
    gs(x, "\\.[0-9]", "")
  })

activate_celldancer <- function(env_pattern = "cellDancer", conda = "~/miniconda3/bin/conda") {
  activate_base(env_pattern, conda = conda)
}

activate_base <- function(env_pattern = "base",
  env_path = pg("conda_env"), conda = pg("conda"))
{
  meta <- dplyr::filter(e(reticulate::conda_list()), grepl(env_pattern, name))
  conda_env <- meta$name[1]
  python <- meta$python[1]
  e(base::Sys.setenv(RETICULATE_PYTHON = python))
  ## e(reticulate::py_config())
  e(reticulate::use_condaenv(conda_env, conda, required = TRUE))
  e(reticulate::import("platform"))
}

treeObj <- function(obj, stops = c("gg", "wrap"), maxdepth = 10L, envir = parent.frame(1)) {
  depth <- 0L
  subsL1 <- ""
  form <- function(x) {
    ifelse(grpl(x, "^[a-zA-Z.][a-zA-Z0-9_.]*$"), x, paste0("`", x, "`"))
  }
  rapp <- function(names) {
    if (depth > maxdepth) {
      return(names)
    }
    depth <<- depth + 1L
    lst <- lapply(names,
      function(oname) {
        obj <- eval(parse(text = oname), envir = envir)
        isThats <- vapply(stops, function(x) is(obj, x), logical(1))
        if (any(isThats)) {
          names <- NULL
        } else if (isS4(obj)) {
          subs <- form(slotNames(obj))
          names <- paste0(oname, "@", subs)
          if (depth == 1L) {
            subsL1 <<- names
          }
          nowDepth <- depth
          names <- rapp(names)
          depth <<- nowDepth
        } else if (is.list(obj)) {
          subs <- form(names(obj))
          names <- ifelse(names(obj) == "",
            paste0(oname, "[[", 1:length(obj), "]]"),
            paste0(oname, "$", subs))
          if (depth == 1L) {
            subsL1 <<- names
          }
          nowDepth <- depth
          names <- rapp(names)
          depth <<- nowDepth
        } else {
          names <- NULL
        }
        c(oname, names)
      })
    unlist(lst)
  }
  res <- rapp(obj)
  res <- res[!is.na(res)]
  unlist(lapply(res,
      function(x) {
        if (any(x == subsL1)) c("", x)
        else x
      }))
}

.view_obj <- function(x, howto = "vsplit") {
  nvimcom:::nvim_viewobj(x, howto = howto)
}

view_obj_for_vim <- function(x, y, view = T) {
  if (isGeneric(x)) {
    if (nchar(y)) {
      classY <- class(eval(parse(text = y), envir = .GlobalEnv))
      if (length(classY) > 1) {
        classY <- classY[1]
      }
      message("Class of target: ", classY)
      obj <- try(selectMethod(x, classY), silent = T)
    } else {
      obj <- NULL
      classY <- NULL
    }
    if (inherits(obj, "try-error") || is.null(classY)) {
      alls <- findMethodSignatures(x)
      if (!nrow(alls)) {
        stop("No methods signatures.")
      } else if (nrow(alls) > 1) {
        if (!is.null(classY)) {
          alls <- alls[alls[, 1] == classY, ]
        }
        if (is.matrix(alls)) {
          useWhich <- menu(apply(alls, 1, paste0, collapse = ", "), title = "Show which?")
        } else if (is.character(alls)) {
          useWhich <- integer(0)
        }
      } else {
        useWhich <- 1L
      }
      if (length(useWhich)) {
        sigs <- alls[useWhich, ]
      } else {
        sigs <- alls
      }
      obj <- selectMethod(x, sigs)
    }
    assign(x <- ".tmpfun", obj, envir = .GlobalEnv)
  }
  if (view) {
    .view_obj(x)
  }
}

.vimShow_complete <- function(obj, ...) {
  .vimShow(treeObj(obj), ...)
}

.vimShow <- function(text, name = "SelectCompletion", howto = "vsplit") {
  if (is.null(text)) {
    message("The `text` is NULL.")
  }
  .prepare_vimShow <- function(text) {
    file <- paste0(Sys.getenv("NVIMR_TMPDIR"), "/Rinsert")
    writeLines(text, file)
  }
  .vimPre <- function(name = "Content", howto = "vsplit") {
    msg <- paste0("ShowRObj(\"", howto, "\", \"", name, "\", \"r\")")
    .C("nvimcom_msg_to_nvim", msg, PACKAGE = "nvimcom")
  }
  .prepare_vimShow(text)
  invisible(.vimPre(name, howto))
}


get_workflow_news <- function(..., data = get_orders()) {
  lst <- as.Date(unlist(list(...)))
  later <- dplyr::filter(data, belong <= max(!!lst))
  former <- dplyr::filter(data, belong < min(!!lst))
  fun <- function(x) unique(x$name)
  x <- fun(plot_workflow_summary(later, T))
  y <- fun(plot_workflow_summary(former, T))
  x[ !x %in% y ]
}

plot_workflow_summary <- function(data, get_jobs_used = F) {
  fun.co_job <- function() {
    alls <- available_signatures("step1")
    name <- gs(alls, "^job_", "")
    nodes <- data.frame(name = name)
    asjobs <- paste0("asjob_", name)
    asjobs <- sapply(asjobs, simplify = F,
      function(x) {
        res <- try(get_fun(x), T)
        if (!inherits(res, "try-error")) {
          available_signatures(x)
        } else {
          NULL
        }
      })
    edges <- as_df.lst(lst_clear0(asjobs))
    edges <- dplyr::mutate(edges, from = gs(type, "^asjob_", ""), to = gs(name, "^job_", ""))
    edges <- dplyr::relocate(edges, from, to)
    edges <- dplyr::filter(edges, from %in% nodes$name, to %in% nodes$name)
    graph <- fast_layout(edges, "linear", nodes = nodes, circular = T)
    p <- ggraph(graph) +
      geom_edge_arc(color = "lightblue", alpha = .5) +
      geom_node_point(aes(size = centrality_degree,
          fill = centrality_degree), alpha = .7, shape = 21) +
      geom_node_text(aes(x = x * 1.1,  y = y * 1.1, 
          label = Hmisc::capitalize(name), angle = -((-node_angle(x,  y) + 90) %% 180) + 90), 
        size = 3, hjust = 'outward', family = "Times") +
      scale_size(range = c(2, 12)) +
      scale_fill_gradientn(colors = color_set()[1:3]) +
      coord_cartesian(xlim = zoRange(graph$x, 1.2), ylim = zoRange(graph$y, 1.2)) +
      theme_graph() +
      theme(text = element_text(family = "Times"))
    wrap(p, 8, 6)
  }
  ###################################
  ###################################
  fun.use_freq <- function(alls) {
    used <- lapply(alls$.dir,
      function(dir) {
        file <- paste0(dir, "/index.Rmd")
        if (file.exists(file)) {
          x <- stringr::str_extract(readLines(file), "job_[a-zA-Z0-9._]+|new_[a-zA-Z0-9._]+|^pb\\.[a-zA-Z0-9._]+")
          isPb <- grpl(x, "^pb")
          x[ isPb ] <- ifelse(!duplicated(x[ isPb ]), x[ isPb ], NA)
          x <- x[!is.na(x)]
          x[ !x %in% c("new_upset", "new_venn", "new_lich", "job_esearch") ]
        }
      })
    alls <- dplyr::mutate(alls,
      .type = dplyr::recode(type,
        `备单业务` = "BI", `固定业务` = "IN", "其他业务" = "Others"))
    names(used) <- paste0(alls$receive_date, "_", 1:length(used), "##", alls$.type)
    whichLen0 <- lengths(used) == 0
    message("No any 'job' detected: \n\t", paste0(alls$.dir[whichLen0], collapse = "\n\t "))
    used <- lst_clear0(used)
    used <- as_df.lst(used)
    used <- dplyr::mutate(used, value = 1,
      order = paste0("Order ", gs(type, "##.*$", "")),
      rank = as.integer(as.Date(type)),
      name = gs(name, "^job_|^new_|(?<=^pb).*", "", perl = T),
      name = gs(name, "^pb$", "Publication"),
      type = gs(type, ".*##(.*)$", "\\1")
    )
    data <- as_tibble(used)
    p.ind.pop <- ggplot(data, aes(x = reorder(order, rank, decreasing = T), y = value)) +
      geom_col(aes(fill = name), position = "stack", width = .6) +
      geom_point(
        data = dplyr::distinct(data, order, type, rank),
        aes(reorder(order, rank, decreasing = T), y = -3, shape = type, color = type),
        size = 5) +
      coord_flip() +
      labs(y = "Order (Date + N)", x = "Frequence", fill = "Name",
        color = "Type", shape = "Type") +
      scale_fill_manual(values = color_set(T)) +
      scale_color_manual(values = rstyle("pal", seed = 6)) +
      theme_minimal()
    p.alls <- new_pie(data$name, fun_text = ggrepel::geom_label_repel)
    p.alls <- wrap(p.alls, 8, 6)
    namel(p.ind.pop, p.alls)
  }
  ###################################
  ###################################
  fun.jobde <- function() {
    omit <- circleGrob(seq(.2, .8, , 3), .5, .07, gp = gpar(fill = "grey20"),
      vp = viewport(, , u(1.5, line), u(1.5, line)))
    ids <- n(n, 5)
    graph <- random_graph(ids, layout = "kk")
    theme <- theme_void() + theme(legend.position = "none")
    layer0 <- ggraph::ggraph(graph)
    scale_size <- scale_size(range = c(1, 3))
    layer5 <- layer0 + 
      ggraph::geom_edge_fan(aes(edge_width = width), color = "lightblue") +
      ggraph::geom_node_point(aes(size = size, fill = name), shape = 21) +
      scale_size + ggsci::scale_fill_npg() + theme
    layer5_ <- ggather(as_grob(layer5))
    ## weight
    weight.chi <- list("..." = 1, tables = 3, plots = 3, params = 3, `...` = 1)
    ## draw
    type <- c("class", "slot", "sub.slot", "function", "custom")
    pal <- MCnebula2:::.as_dic(color_set(), type, na.rm = T)
    grobs.chi <- lst_grecti(names(weight.chi), pal, "custom", grecti2)
    grobs.chi$... <- omit
    ## signal grid layout
    n <- 3
    seq <- seq(.2, .8, length.out = n)
    grid <- segmentsGrob(c(rep(0, n), seq),
      c(seq, rep(1, n)),
      c(rep(1, n), seq),
      c(seq, rep(0, n)), gp = gpar(lty = "dashed"))
    grid <- ggather(clipGrob(height = .8), rectGrob(), grid, vp = viewport(, , .8, .8))
    grobs.chi$tables %<>% into(grid)
    ## signal viewport
    n <- 6
    seq <- seq(135, 0, length.out = n)
    seq2 <- seq(45, 0, length.out = n)
    vps <- lapply(1:length(seq),
      function(n, fx = .5, fy = .5, x = .5) {
        ang <- seq[n]
        ang2 <- seq2[n]
        vp <- viewport(cospi(ang / 180) * fx + x, sinpi(ang / 180) * fy + x,
          u(2, line), u(2, line), angle = ang2,
          just = c("centre", "bottom"))
        ggather(rectGrob(gp = gpar(fill = "lightyellow", col = pal[["sub.slot"]])),
          gtext(paste0("n", rev(1:length(seq))[n]),
            gpar(fontface = "plain"), form = F), vp = vp)
      })
    vps <- do.call(ggather, c(vps, list(vp = viewport(, .1, .5, .5))))
    grobs.chi$params %<>% into(vps)
    ## signal ggset
    ggsets <- frame_col(c(n = 1, x = 1, mn = 1),
      list(n = gtext("n", list(cex = 2.5), form = F),
        x = gtext("×", list(cex = 2, font = c(plain = 1))),
        mn = into(glayer(3), layer5_)))
    ggsets <- ggather(ggsets, vp = viewport(, , .7, .7))
    grobs.chi$plots %<>% into(ggsets)
    ## child... gather
    chi <- frame_row(weight.chi, grobs.chi)
    chi <- ggather(chi, vp = viewport(, , .9, .95))
    job <- grecti2("Job")
    job <- into(job, chi)
    route <- as_network(
      list("Job:step1",
        "step1:step2",
        "step2:step3",
        "step3:step4",
        "step4:..."
        ), "tree"
    )
    p.route <- as_grob(flowChart(route, 1.1, 1))
    job_de <- xf(job = 1, p.route = .4)
    p.jobde <- wrap(job_de, 5, 4)
    p.jobde
  }
  message("Plot use frequency")
  p.use_freq <- fun.use_freq(data)
  if (get_jobs_used) {
    jobs_used <- p.use_freq$p.ind.pop$data
    return(jobs_used)
  }
  message("Plot correlation of jobs")
  p.co_job <- fun.co_job()
  message("Plot job demonstration")
  p.jobde <- fun.jobde()
  p.fr1 <- frame_col(c(p.jobde = 1, p.co_job = 1.5),
    list(p.jobde = ggather(p.jobde@data, vp = viewport(, , .9, .8)),
      p.co_job = as_grob(p.co_job@data)))
  p.fr1 <- wrap(p.fr1, 10, 5)
  p.fr2 <- frame_col(c(p.alls = 1, p.ind.pop = 1),
    list(p.alls = ggather(p.use_freq$p.alls@data, vp = viewport(, , .9, .9)),
      p.ind.pop = as_grob(p.use_freq$p.ind.pop)))
  p.fr2 <- wrap(p.fr2, 14, 7)
  lst <- c(namel(p.co_job, p.jobde), p.use_freq, namel(p.fr1, p.fr2))
  jobs_used <- p.use_freq$p.ind.pop$data
  attr(lst, "jobs_used") <- jobs_used
  ##########################################
  p.tags <- new_pie(unlist(data$tags.cn))
  lst$p.tags <- wrap(p.tags, 10, 8)
  lst
}

# ==========================================================================
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.local_db <- setClass("local_db", 
  contains = c(),
  representation = representation(
    db = "df",
    idcol = "character",
    db_file = "character",
    query = "ANY"),
  prototype = NULL)

new_db <- function(db_file, idcol) {
  path <- dirname(db_file)
  if (!dir.exists(path)) {
    dir.create(path)
  }
  .local_db(db_file = db_file, idcol = idcol)
}

setMethod("not", signature = c(x = "local_db"),
  function(x, ids, force = F){
    ids <- unique(ids)
    if (nrow(x@db) & !force) {
      stop("The `x` is not an empty data object.")
    }
    if (file.exists(x@db_file)) {
      x@db <- readRDS(x@db_file)
      if (!is(x@db, "tbl_df")) {
        x@db <- dplyr::as_tibble(x@db)
      }
      x@query <- ids[ !ids %in% x@db[[ x@idcol ]] ]
    } else {
      x@query <- ids
    }
    return(x)
  })

setMethod("upd", signature = c(x = "local_db"),
  function(x, data, ids = x@query){
    idcol <- x@idcol
    if (!is(data, "df")) {
      stop('is(data, "df") == F')
    }
    if (!is.null(ids)) {
      if (nrow(data)) {
        notIn <- ids[ !ids %in% data[[ idcol ]] ]
        if (length(notIn)) {
          args <- c(list(data), nl(idcol, list(notIn)))
          data <- do.call(tibble::add_row, args)
        }
      } else {
        data <- do.call(tibble::tibble, nl(idcol, list(ids)))
      }
    }
    if (nrow(data)) {
      data <- dplyr::filter(data, !(!!rlang::sym(idcol) %in% !!x@db[[ idcol ]]))
    }
    if (nrow(data)) {
      x@db <- dplyr::bind_rows(x@db, tibble::as_tibble(data))
      saveRDS(x@db, x@db_file)
    }
    return(x)
  })

setMethod("vis", signature = c(x = "df"),
  function(x, n = 100, width = 200) {
    print(x, n = n, width = width)
  })

setMethod("res", signature = c(x = "local_db"),
  function(x, what){
    dplyr::filter(as_tibble(x@db), !!rlang::sym(x@idcol) %in% !!what)
  })

.get_savesInternalRds <- function(file) {
  if (file == "workflow.rds") {
    rds <- ".injobs.rds"
  } else {
    rds <- paste0(".injobs.", get_realname(file), ".rds")
  }
}

saves <- function(file = "workflow.rdata", saveInjobs = T, ...) {
  injobs <- getOption("internal_job")
  if (!is.null(injobs) && saveInjobs) {
    rds <- .get_savesInternalRds(file)
    saveRDS(injobs, rds)
  }
  if (file.exists(file)) {
    isThat <- usethis::ui_yeah("The `file` exists, overwrite that?")
    if (!isThat) {
      return(NULL)
    }
  }
  save.image(file, ...)
}

loads <- function(file = "workflow.rdata", ...) {
  if (file.exists(file)) {
    load(file, envir = .GlobalEnv)
  }
  rds <- .get_savesInternalRds(file)
  if (file.exists(rds)) {
    injobs <- readRDS(rds)
    options(internal_job = injobs)
  }
  writeAllCompletion()
  # if (isNamespaceLoaded("nvimcom")) {
    # pkgload::unload("nvimcom")
    # require("nvimcom")
  # }
}

remotejob <- function(wd, remote = "remote") {
  .job(params = list(set_remote = T, wd = wd, remote = remote))
}

set_prefix <- function(wd, db, op) {
  options(wd_prefix = wd, db_prefix = db, op_prefix = op)
}

.prefix <- function(path, name = c("wd", "db", "op"), warning = T)
{
  name <- match.arg(name)
  pr <- getOption(That <- paste0(name, "_prefix"), NULL)
  if (is.null(pr)) {
    stop("`", That, "` not set. Use `set_prefix` or `options` to set that.")
  }
  if (missing(path)) {
    res <- pr
  } else {
    res <- file.path(pr, path)
  }
  if (!file.exists(res) & warning) {
    warning("The file not exists in local")
  }
  res
}

stop_debug <- function(x) {
  Terror <<- x
  stop("Stop for debugging. Use `Terror` to get data.")
}

setMethod("show", signature = c(object = "character"),
  function(object){
    if (length(object) > 200) {
      message(crayon::silver("Length: ", length(object)))
      object <- head(object, n = 200)
    }
    if (any(tooMany <- (nchar(object) > 1000))) {
      message(crayon::silver("Chunk strings too long."))
      object[ tooMany ] <- paste0(substr(object, 1, 1000), crayon::silver("..."))
    }
    print(object)
  })

genes <- function(x) {
  if (any(grpl(x, "\\.[0-9]+$"))) {
    message("Trunc ther version info after Symbol.")
    x <- gs(x, "\\.[0-9]+$", "")
  }
  x
}
