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
    step = "integer"
    ),
  prototype = prototype(step = -1L))

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
setMethod("plots", signature = c(x = "job", step = "double"),
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

setGeneric("step0", 
  function(x, ...) standardGeneric("step0"))
setMethod("step0", signature = c(x = "missing"),
  function(x){
    message("This step do nothing but description.\n")
  })

setGeneric("step1", 
  function(x, ...) standardGeneric("step1"))
setMethod("step1", signature = c(x = "job"),
  function(x){
    if (x@step >= 1)
      stop("x@step >= 1. This step has been executed.")
    x@step <- 1L
    x
  })

setGeneric("step2", 
  function(x, ...) standardGeneric("step2"))
setMethod("step2", signature = c(x = "job"),
  function(x){
    if (x@step >= 2)
      stop("x@step >= 2. This step has been executed.")
    x@step <- 2L
    x
  })

setGeneric("step3", 
  function(x, ...) standardGeneric("step3"))
setMethod("step3", signature = c(x = "job"),
  function(x){
    if (x@step >= 3)
      stop("x@step >= 3. This step has been executed.")
    x@step <- 3L
    x
  })

setGeneric("step4", 
  function(x, ...) standardGeneric("step4"))
setMethod("step4", signature = c(x = "job"),
  function(x){
    if (x@step >= 4)
      stop("x@step >= 4. This step has been executed.")
    x@step <- 4L
    x
  })

setGeneric("step5", 
  function(x, ...) standardGeneric("step5"))
setMethod("step5", signature = c(x = "job"),
  function(x){
    if (x@step >= 5)
      stop("x@step >= 5. This step has been executed.")
    x@step <- 5L
    x
  })

setGeneric("step6", 
  function(x, ...) standardGeneric("step6"))
setMethod("step6", signature = c(x = "job"),
  function(x){
    if (x@step >= 6)
      stop("x@step >= 6. This step has been executed.")
    x@step <- 6L
    x
  })

setGeneric("step7", 
  function(x, ...) standardGeneric("step7"))
setMethod("step7", signature = c(x = "job"),
  function(x){
    if (x@step >= 7)
      stop("x@step >= 7. This step has been executed.")
    x@step <- 7L
    x
  })

setGeneric("step8", 
  function(x, ...) standardGeneric("step8"))
setMethod("step8", signature = c(x = "job"),
  function(x){
    if (x@step >= 8)
      stop("x@step >= 8. This step has been executed.")
    x@step <- 8L
    x
  })

setGeneric("step9", 
  function(x, ...) standardGeneric("step9"))
setMethod("step9", signature = c(x = "job"),
  function(x){
    if (x@step >= 9)
      stop("x@step >= 9. This step has been executed.")
    x@step <- 9L
    x
  })

step_message <- function(...) {
  str <- paste0(unlist(list(...)), collapse = "")
  textSh(str, pre_trunc = F)
}
