# ==========================================================================
# for wrap various of analysis workflow contains many steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
    info = "ANY"
    ),
  prototype = prototype(step = 0L))

.marker_list <- setClass("marker_list", 
  contains = c(),
  representation = representation(level = "numeric", marker = "list"),
  prototype = NULL)

setMethod("show", 
  signature = c(object = "marker_list"),
  function(object){
    .show(object)
  })

setMethod("show", 
  signature = c(object = "job"),
  function(object){
    message("A workflow of '", class(object), "' in step (done): ", object@step)
    message("Object size: ", obj.size(object))
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


setMethod("object", signature = c(x = "job"),
  function(x){
    x@object
  })
setReplaceMethod("object", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })

setMethod("plots", signature = c(x = "job"),
  function(x){
    x@plots[[ x@step ]]
  })
setMethod("plots", signature = c(x = "job", step = "numeric"),
  function(x, step){
    x@plots[[ step ]]
  })
setReplaceMethod("plots", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })

setMethod("tables", signature = c(x = "job"),
  function(x){
    x@tables[[ x@step ]]
  })
setMethod("tables", signature = c(x = "job", step = "double"),
  function(x, step){
    x@tables[[ step ]]
  })
setMethod("tables", signature = c(x = "job", step = "double"),
  function(x, step){
    x@tables[[ step ]]
  })
setReplaceMethod("tables", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })


setMethod("params", signature = c(x = "job"),
  function(x){
    x@params
  })
setReplaceMethod("params", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
  })


setMethod("others", signature = c(x = "job"),
  function(x){
    x@others[[ x@step ]]
  })
setReplaceMethod("others", signature = c(x = "job"),
  function(x, value){
    initialize(x, object = value)
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

e <- function(expr, text = NULL) {
  expr <- substitute(expr)
  names <- stringr::str_extract(deparse(expr), "[a-zA-Z0-9_.]*:::?[a-zA-Z0-9_.]*")
  names <- names[ !is.na(names) ]
  if (!length(names))
    name <- "FUNCTION"
  else name <- names[1]
  if (is.null(text))
    cli::cli_alert_info(name)
  else
    cli::cli_alert_info(paste0(name, " (", text, ")" ))
  suppressMessages(eval(expr, envir = parent.frame(1)))
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
      str <- strsplit(str, "")[[1]]
      sig <- 0L
      pair <- 0L
      char <- c()
      for (i in str) {
        if (i == "(") {
          pair <- 1L
          sig <- sig + 1L
        } else if (i == ")") {
          sig <- sig - 1L
        }
        if (sig == 0L & pair)
          break
        else if (pair) {
          if (i == "(" & sig == 1L)
            next
          char <- c(char, i)
        }
      }
      args <- paste(char, collapse = "")
      message(crayon::green(gsub("#", ", ", methods@names[[n]])), ":")
      message(strwrap(args, indent = 4))
    })
  cli::cli_h1("Methods parameters")
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

setGeneric("map", 
  function(x, ref, ...) standardGeneric("map"))

setGeneric("getsub", 
  function(x, ...) standardGeneric("getsub"))

setGeneric("active", 
  function(x, ...) standardGeneric("active"))


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
