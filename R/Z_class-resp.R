# ==========================================================================
# setGeneric
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMissing <- function(generic, ..., .SIG = "missing"){
    args <- list(...)
    sig <- getGeneric(generic)@signature
    res <- vapply(sig, FUN.VALUE = "character",
      function(name){
        if (is.null(args[[ name ]]))
          .SIG
        else
          args[[ name ]]
      })
    names(res) <- sig
    return(res)
  }

setGeneric("new_resp", 
  function(..., yaml) standardGeneric("new_resp"))

setGeneric("resp_args", 
  function(x) standardGeneric("resp_args"))
setGeneric("resp_args<-", 
  function(x, value) standardGeneric("resp_args<-"))
setGeneric("chrunk_args", 
  function(x) standardGeneric("chrunk_args"))
setGeneric("chrunk_args<-", 
  function(x, value) standardGeneric("chrunk_args<-"))
setGeneric("%<in%", 
  function(x, y) standardGeneric("%<in%"))

# ==========================================================================
# class defination
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.yaml_default <- function(style = c("resp_default")) {
  style <- match.arg(style)
  readLines(paste0(.expath, "/", style, ".yml"))
}

.resp <- 
  setClass("resp", 
    contains = c("layerSet"),
    representation = 
      representation(
        yaml = "character",
        resp_args = "list",
        chrunk_args = "list"),
      prototype = prototype(yaml = "")
  )

# ==========================================================================
# validity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setValidity("resp", 
  function(object){
    recepts <- c("section", "heading", "code_block", "character")
    tip <- paste0("'layer' in 'resp' must either be: ",
      paste0("'", recepts, "'", collapse = ", "))
    validate_class_in_list(layers(object), recepts, tip)
  })

# ==========================================================================
# method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("resp_args", 
  signature = c(x = "resp"),
  function(x){
    x@resp_args
  })

setReplaceMethod("resp_args", 
  signature = c(x = "resp"),
  function(x, value){
    initialize(x, resp_args = value)
  })

setMethod("chrunk_args", 
  signature = c(x = "resp"),
  function(x){
    x@chrunk_args
  })

setReplaceMethod("chrunk_args", 
  signature = c(x = "resp"),
  function(x, value){
    initialize(x, chrunk_args = value)
  })

#' @exportMethod new_resp
#' @aliases new_resp
#' @description \code{new_resp}: Create a [resp-class] object.
#' @param ... An arbitrary number of [heading-class],
#' [section-class] or [code_block-class] in sequence.
#' Specially, \code{NULL} can be passed herein, but would be ignored.
#' @param yaml character. Passed to .Rmd for setting format of documentation.
#' @rdname resp-class
setMethod("new_resp", 
  signature = c(yaml = "character"),
  function(..., yaml){
    layers <- list(...)
    layers <- layers[!vapply(layers, is.null, logical(1))]
    .resp(yaml = yaml, layers = layers)
  })

#' @exportMethod new_resp
#' @description \code{new_resp()}: get the default parameters for the method \code{new_resp}.
#' @description \code{new_resp(x, ...)}: use the default parameters whatever 'missing'
#' while performing the method \code{new_resp}.
#' @rdname resp-class
setMethod("new_resp", 
  signature = setMissing("new_resp",
    yaml = "missing"),
  function(...){
    args <- list(yaml = .yaml_default())
    if (missing(...))
      return(args)
    else
      reCallMethod("new_resp", args, ...)
  })

#' @exportMethod call_command
#' @aliases call_command
#' @description \code{call_command}: Format 'resp' object as character, which can be output
#' by \code{writeLines()} function as '.Rmd' file and than use \code{rmarkdown::render} output
#' as pdf, html, or other format files.
#' @family call_commands
#' @seealso [writeLines()], [rmarkdown::render()]...
#' @rdname resp-class
setMethod("call_command", 
  signature = c(x = "resp"),
  function(x){
    yaml <- c("---", yaml(x), "---")
    layers <- unlist(lapply(layers(x), call_command))
    c(yaml, "", layers)
  })

setMethod("call_command", 
  signature = c(x = "character"),
  function(x){ 
    if (tail(x, n = 1) != "")
      c(x, "")
    else
      x
  })


# ==========================================================================
# functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# ==========================================================================
# operator
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setMethod("%<in%", 
  signature = c(x = "resp"),
  function(x, y){
  })


