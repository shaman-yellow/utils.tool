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

setMethod("show", 
  signature = c(object = "job"),
  function(object){
    message("A workflow of '", class(object), "' in step (done): ", object@step)
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
  function(x, ...) standardGeneric("step0"))
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

setMethod("step0", signature = c(x = "job"),
  function(x){
    message("This step do nothing but description.")
    message(crayon::silver("Use `step0()` to show all available workflow.\n"))
  })

lapply(1:12,
  function(n) {
    STEP <- paste0("step", n)
    setGeneric(STEP, function(x, ...) x)
    setMethod(STEP, signature = c(x = "job"),
      function(x, ...){
        if (!existsMethod(STEP, class(x)))
          stop("existsMethod(", STEP, ", class(x)) == F")
        if (x@step >= n)
          stop("x@step >= ", n, ". This step has been executed.")
        if (x@step < n - 1)
          stop("x@step <", n - 1, ". The previous step has not been executed.")
        if (STEP == "step1") {
          if (!is.null(x@info)) {
            message(crayon::yellow("Workflow information:"))
            step_message(x@info)
          }
        }
        message(crayon::green("Running", STEP), ":")
        x@step <- n
        x
      })
  })

step_message <- function(..., show_end = "Description end") {
  str <- paste0(unlist(list(...)), collapse = "")
  str <- gs(str, "red\\{\\{(.*?)\\}\\}", "\033[31m\\1\033[39m")
  str <- gs(str, "grey\\{\\{(.*?)\\}\\}", "\033[90m\\1\033[39m")
  str <- gs(str, "yellow\\{\\{(.*?)\\}\\}", "\033[33m\\1\033[39m")
  str <- gs(str, "blue\\{\\{(.*?)\\}\\}", "\033[34m\\1\033[39m")
  textSh(str, pre_trunc = F)
  if (!is.null(show_end))
    cli::cli_h1(show_end)
}

