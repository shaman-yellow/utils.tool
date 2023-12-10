# ==========================================================================
# for wrap various of analysis workflow contains many steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @exportClass job
.job <- setClass("job", 
  contains = c(),
  representation = representation(
    "VIRTUAL",
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
    sig = "character"
    ),
  prototype = prototype(step = 0L,
    params = list(wd = ".", remote = "remote")
    ))

setValidity("job", 
  function(object){
    class <- class(object)
    options(method_name = class)
    if (is(object, "job")) T else F
  })

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
  writeLines(Hmisc::capitalize(gs(lab(x), " ", "-")))
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
        lab(obj) <- gs(gs(paste(sig, group, body, suffix), "[ ]+", " "), "^[ \\-]|[ \\-]$", "")
        return(obj)
      })
  } else {
    lab(x) <- gs(gs(paste(sig, group, body, suffix), "[ ]+", " "), "^[ \\-]|[ \\-]$", "")
  }
  return(x)
}

setMethod("sig", signature = c(x = "job"),
  function(x){
    x@sig
  })
setReplaceMethod("sig", signature = c(x = "job"),
  function(x, value){
    initialize(x, sig = value)
  })

#' @exportMethod pg
setMethod("pg", signature = c(x = "job"),
  function(x, recode = getOption("pg_remote_recode", pg_remote_recode())){
    if (is.remote(x)) {
      recode <- recode[[ x@pg ]]
      if (is.null(recode))
        stop("No program running command obtained.")
      recode
    } else {
      x@pg
    }
  })

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
    x <- standardGeneric("set_remote")
    set_remote.default(x, tmpdir = "/data/hlc/tmp",
      map_local = paste0(gs(class(x), "^job_", ""), "_local"),
      remote = "remote"
    )
  })

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

setGeneric("step1",
  function(x, ...) {
    if (identical(sig(x), character(0))) {
      name <- substitute(x)
      sig <- toupper(gs(gs(gs(name, "^[a-zA-Z0-9]+", ""), "\\.", " "), "^[ ]*", ""))
      sig(x) <- sig
    }
    x$seed <- sample(1:100, 1)
    x <- checkAddStep(x, 1L)
    x <- standardGeneric("step1")
    stepPostModify(x, 1)
  })

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
  x <- convertPlots(x, n)
  if (length(x@plots) >= n) {
    if (!is.null(x@plots[[ n ]]))
      names(x@plots)[ n ] <- paste0("step", n)
  }
  if (length(x@tables) >= n) {
    if (!is.null(x@tables[[ n ]]))
      names(x@tables)[ n ] <- paste0("step", n)
  }
  x
}

convertPlots <- function(x, n) {
  if (length(x@plots) >= n)
    x@plots[[ n ]] <- rapply2_asWrap(x@plots[[ n ]], "gg.obj")
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
    names <- names(x)
    x <- lapply(1:length(x),
      function(n) {
        name <- names[[ n ]]
        obj <- x[[ n ]]
        if (!is.null(obj))
          attr(obj, ".name") <- name
        obj
      })
    names(x) <- names
    lapply(x, rapply2_asWrap, class = class)
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

checkAddStep <- function(x, n) {
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
  message(crayon::green("Running", STEP), ":")
  x@step <- n
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
  if (isNamespaceLoaded("future")) {
    e(pkgload::unload("future"))
  }
  options(future.globals.maxSize = Inf)
  e(future::plan(future::multicore, workers = workers))
  x <- fun(x)
  e(future::plan(future::sequential))
  x
}

show_nonstandardGenericFunction <- selectMethod(
  "show", "nonstandardGenericFunction"
)

show_standardGeneric <- selectMethod(
  "show", "standardGeneric"
)

.show_method <- function(object, default, filter_by_options){
  if (!hasMethods(object, package = environmentName(topenv()))) {
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

setGeneric("fill", 
  function(x, ...) standardGeneric("fill"))

setGeneric("upd", 
  function(x, ...) standardGeneric("upd"))

setGeneric("not", 
  function(x, ...) standardGeneric("not"))

setMethod("not", signature = c(x = "job"),
  function(x){
    validObject(x)
    return(x)
  })

setMethod("upd", signature = c(x = "ANY"),
  function(x){
    new <- new(class(x))
    new@object <- x@object
    new@params <- x@params
    new@plots <- x@plots
    new@tables <- x@tables
    new@step <- x@step
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

remoteRun <- function(..., path, run_after_cd = NULL,
  postfix = NULL, remote = "remote", tmpdir = "/data/hlc/tmp", x)
{
  expr <- paste0(unlist(list(...)), collapse = "")
  if (missing(x)) {
    x <- get("x", parent.frame(1))
  }
  if (missing(remote)) {
    remote <- x@params$remote
  }
  if (missing(postfix)) {
    postfix <- x@params$postfix
  }
  if (!is.null(postfix)) {
    expr <- postfix(expr)
  }
  if (missing(tmpdir)) {
    tmpdir <- x@params$tmpdir
  }
  if (!is.null(tmpdir))
    expr <- c(paste0("export TMPDIR=", tmpdir), expr)
  if (missing(path)) {
    path <- x@params$wd
  }
  if (!is.null(run_after_cd)) {
    expr <- c(run_after_cd, expr)
  }
  if (!is.null(path)) {
    expr <- c(paste0("cd ", path), expr)
  }
  script <- tempfile("remote_Script_", fileext = ".sh")
  writeLines(expr, script)
  writeLines(crayon::yellow(paste0("The script file for remote is: ", script)))
  system(paste0("ssh ", remote, " < ", script))
}

setGeneric("is.remote", 
  function(x, ...) standardGeneric("is.remote"))

setMethod("is.remote", signature = c(x = "job"),
  function(x){
    if (!is.null(x@params$set_remote)) T
    else F
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

pg_remote_recode <- function() {
  list(
    qiime = "~/miniconda3/bin/conda run -n qiime2 qiime",
    fastp = "conda run -n base fastp",
    bcftools = "conda run -n base bcftools",
    elprep = "conda run -n base elprep",
    biobakery_workflows = "conda run -n biobakery biobakery_workflows"
  )
}

setMethod("map", signature = c(x = "df"),
  function(x, ref, y, y.ref, y.get, rename = T){
    x[[ ref ]] <- y[[ y.get ]][match(x[[ ref ]], y[[ y.ref ]])]
    if (rename)
      colnames(x)[ colnames(x) == ref ] <- y.get
    x
  })

setMethod("map", signature = c(x = "list"),
  function(x, y, y.ref, y.get){
    lapply(x,
      function(x) {
        y[[ y.get ]][match(x, y[[ y.ref ]])]
      })
  })

setMethod("gname", signature = c(x = "character"),
  function(x){
    gs(x, "\\.[0-9]", "")
  })

move_rds <- function(from, to, exclude = ".items.rds") {
  lapply(from,
    function(path) {
      dirname <- get_filename(gs(path, "/$", ""))
      target <- paste0(to, "/", dirname)
      dir.create(target, F)
      files <- list.files(path, pattern = ".rds$|.RData$|.Rdata$", full.names = T, all.files = T)
      thefilenames <- get_filename(files)
      files <- files[ !thefilenames %in% exclude ]
      lapply(files,
        function(file) {
          system(paste0("cp ", file, " -t ", target))
          file.remove(file)
        })
    })
}

activate_celldancer <- function(env_pattern = "cellDancer", conda = "~/miniconda3/bin/conda") {
  activate_base(env_pattern, conda = conda)
}

activate_base <- function(env_pattern = "base", env_path = "~/miniconda3/envs/", conda = "~/miniconda3/bin/conda")
{
  meta <- dplyr::filter(e(reticulate::conda_list()), grepl(env_pattern, name))
  conda_env <- meta$name[1]
  python <- meta$python[1]
  e(base::Sys.setenv(RETICULATE_PYTHON = python))
  ## e(reticulate::py_config())
  e(reticulate::use_condaenv(conda_env, conda, required = TRUE))
  e(reticulate::import("platform"))
}

.view_obj <- function(x, howto = "vsplit") {
  nvimcom:::nvim_viewobj(x, howto = howto)
}

view_obj_for_vim <- function(x, y) {
  if (isGeneric(x)) {
    obj <- selectMethod(x, class(get(y, envir = .GlobalEnv)))
    assign(name <- ".tmpfun", obj, envir = .GlobalEnv)
    .view_obj(name)
  } else {
    .view_obj(x)
  }
}
