# ==========================================================================
# utilites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases utilites
#'
#' @title utilites for programming
#'
#' @description This is a combination of tools that are not always used.
#'
#' @name utilites
NULL
#> NULL

reCallMethod <- 
  function(funName, args, ...){
    arg.order <- unname(getGeneric(funName)@signature)
    args.missing <- !arg.order %in% names(args)
    if (any(args.missing)) {
      args.missing <- arg.order[args.missing]
      args.missing <- sapply(args.missing, simplify = F,
                             function(x) structure(0L, class = "missing"))
      args <- c(args, args.missing)
    }
    args <- lapply(arg.order, function(i) args[[i]])
    sig <- get_signature(args)
    method <- selectMethod(funName, sig)
    last_fun <- sys.function(sys.parent())
    n <- 0
    while (identical(last_fun, method@.Data, ignore.environment = T)) {
      if (n == 0) {
        mlist <- getMethodsForDispatch(getGeneric(funName))
      }
      n <- n + 1
      rm(list = paste0(method@defined, collapse = "#"), envir = mlist)
      method <- selectMethod(funName, sig, mlist = mlist)
    }
    expr <- paste0("method@.Data(",
                   paste0(paste0(arg.order, " = args[[",
                                 1:length(arg.order), "]]"),
                          collapse = ", "),
                   ", ...)")
    eval(parse(text = expr))
  }

setMissing <- 
  function(generic, ..., .SIG = "missing"){
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

.fresh_param <- 
  function(default, args){
    if (missing(args))
      args <- as.list(parent.frame())
    args <- args[ !vapply(args, is.name, T) ]
    sapply(unique(c(names(default), names(args))),
           simplify = F,
           function(name){
             if (any(name == names(args)))
               args[[ name ]]
             else
               default[[ name ]]
           })
  }

.fresh_param2 <- 
  function(default, args){
    if (missing(args))
      return(default)
    if (length(args) == 0)
      return(default)
    .fresh_param(default, args)
  }

.fresh_param2f <- 
  function(default, args, class = "gpar"){
    structure(.fresh_param2(default, args), class = class)
  }

.check_columns <- 
  function(obj, lst, tip){
    if (!is.data.frame(obj))
      stop(paste0("'", tip, "' must be a 'data.frame'."))
    lapply(lst, function(col){
             if (is.null(obj[[ col ]]))
               stop(paste0("'", tip, "' must contains a column of '", col, "'."))
           })
  }

.message_info_viewport <- 
  function(info = "info"){
    .message_info(info, "current.viewport:",
                  paste0("\n\t", paste0(grid::current.viewport())))
  }

.message_info <- 
  function(main, sub, arg = NULL, sig = "##"){
    message(sig, " ", main, ": ", sub, " ", arg)
  }

.suggest_bio_package <- 
  function(pkg){
    if (!requireNamespace(pkg, quietly = T))
      stop("package '", pkg, "' not installed. use folloing to install:\n",
           '\nif (!require("BiocManager", quietly = TRUE))',
           '\n\tinstall.packages("BiocManager")',
           '\nBiocManager::install("', pkg, '")\n\n')
  }

#' @export tmp_pdf
#' @aliases tmp_pdf
#' @description \code{tmp_pdf}: ...
#' @rdname utilites
tmp_pdf <- function() {
  paste0(tempdir(), "/tmp_pdf.pdf")
}

#' @export op
#' @aliases op
#' @description \code{op}: ...
#' @rdname utilites
op <- function(file) {
  system(paste0("xdg-open ", file))
}

#' @export .cairosvg_to_grob
#' @aliases .cairosvg_to_grob
#' @description \code{.cairosvg_to_grob}: Convert cairo svg to 'grob'.
#' @rdname utilites
.cairosvg_to_grob <- 
  function(path){
    grImport2::grobify(grImport2::readPicture(path))
  }

#' @export .as_dic
#' @aliases .as_dic
#' @description \code{.as_dic}: ...
#' @rdname utilites
.as_dic <- 
  function(vec, names, default,
           fill = T, as.list = T, na.rm = F){
    if (is.null(names(vec)))
      names(vec) <- names[1:length(vec)]
    if (fill) {
      if (any(!names %in% names(vec))) {
        ex.names <- names[!names %in% names(vec)]
        ex <- rep(default, length(ex.names))
        names(ex) <- ex.names
        vec <- c(vec, ex)
      }
    }
    if (as.list) {
      if (!is.list(vec))
        vec <- as.list(vec)
    }
    if (na.rm) {
      vec <- vec[!is.na(names(vec))]
    }
    vec
  }

#' @export fill_list
#' @aliases fill_list
#' @description \code{fill_list}: ...
#' @rdname utilites
fill_list <- function(names, vec, default = vec[1]) {
  .as_dic(vec, names, default, fill = T, as.list = F, na.rm = F)
}

#' @export n
#' @aliases n
#' @description \code{n}: ...
#' @rdname utilites
n <- function(name, n){
  name <- as.character(substitute(name))
  paste0(name, 1:n)
}

#' @export namel
#' @aliases namel
#' @description \code{namel}: ...
#' @rdname utilites
namel <- function(...){
  call <- substitute(list(...))
  lst <- list(...)
  if (is.null(names(call))) {
    names <- lapply(2:length(call), function(n) as.character(call)[n])
  } else {
    names <- vapply(2:length(call), FUN.VALUE = character(1),
                    function(n) {
                      if (names(call)[n] == "")
                        as.character(call)[n]
                      else
                        names(call)[n]
                    })
  }
  names(lst) <- names
  lst
}

#' @export repSuffix
#' @aliases repSuffix
#' @description \code{repSuffix}: ...
#' @rdname utilites
repSuffix <- 
  function(chs, anno = ".rep."){
    gsub(paste0(anno, 1, "$"), "",
         vapply(1:length(chs), FUN.VALUE = character(1),
                function(n){
                  paste0(chs[n], anno, length(chs[1:n][ chs[1:n] == chs[n] ]))
                }))
  }

#' @export %>%
#' @aliases %>%
#' @description \code{%>%}: ...
#' @rdname utilites
`%>%` <- magrittr::`%>%`

#' @export %<>%
#' @aliases %<>%
#' @description \code{%<>%}: ...
#' @rdname utilites
`%<>%` <- magrittr::`%<>%`

#' @export .expath
#' @aliases .expath
#' @description \code{.expath}: ...
#' @rdname utilites
.expath <- function() {
  .expath <- system.file("extdata", ".", package = gsub("^.*:", "", environmentName(topenv())))
  assign('.expath', .expath, envir = topenv())
}

.onLoad <- function(libname, pkgname) {
  .expath()
  .expathsvg()
  .check_external_svg()
}

#' @export agroup
#' @aliases agroup
#' @description \code{agroup}: ...
#' @rdname utilites
agroup <- function(group, value, FUN.VALUE = character(1)) {
  ug <- unique(group)
  if (length(ug) > length(value))
    stop( "the length of 'value' not enough to assign" )
  dic <- .as_dic(value, ug, fill = F, na.rm = F)
  vapply(group, function(g) dic[[g]], FUN.VALUE)
}

#' @export write_tsv
#' @aliases write_tsv
#' @description \code{write_tsv}: ...
#' @rdname utilites
write_tsv <-
  function(x, filename, col.names = T, row.names = F){
    write.table(x, file = filename, sep = "\t",
                col.names = col.names, row.names = row.names, quote = F)
  }

#' @export read_tsv
#' @aliases read_tsv
#' @description \code{read_tsv}: ...
#' @rdname utilites
read_tsv <- function(path){
  file <- data.table::fread(input = path, sep = "\t",
                            header = T, quote = "", check.names = F)
  return(file)
}

#' @export mapply_rename_col
#' @aliases mapply_rename_col
#' @description \code{mapply_rename_col}: ...
#' @rdname utilites
mapply_rename_col <- 
  function(
           mutate_set,
           replace_set,
           names,
           fixed = F
           ){
    envir <- environment()
    mapply(mutate_set, replace_set,
           MoreArgs = list(envir = envir, fixed = fixed),
           FUN = function(mutate, replace, envir,
                          fixed = F, names = get("names", envir = envir)){
             names <- gsub(mutate, replace, names, perl = ifelse(fixed, F, T), fixed = fixed)
             assign("names", names, envir = envir)
           })
    return(names)
  }

#' @export turn_vector
#' @aliases turn_vector
#' @description \code{turn_vector}: ...
#' @rdname utilites
turn_vector <- function(vec) {
  names <- names(vec)
  names(vec) <- unname(vec)
  vec[] <- names
  vec
}

#' @export group_switch
#' @aliases group_switch
#' @description \code{group_switch}: ...
#' @rdname utilites
group_switch <- function(data, meta.lst, by) {
  if (!is.character(data[[ by ]]))
    stop( "is.character(data[[ by ]]) == F" )
  meta <- unlist(meta.lst)
  names(meta) <- rep(names(meta.lst), lengths(meta.lst))
  meta <- as.list(turn_vector(meta))
  data <- data[data[[by]] %in% names(meta), ]
  group <- data.frame(order = 1:length(data[[ by ]]), col = data[[ by ]])
  group <- split(group, ~ col)
  group <- lapply(names(group),
                  function(name){
                    data <- group[[ name ]]
                    data$col <- meta[[ name ]]
                    return(data)
                  })
  group <- data.table::rbindlist(group)
  group <- group[order(group$order), ]$col
  split(data, group)
}

#' @export .find_and_sort_strings
#' @aliases .find_and_sort_strings
#' @description \code{.find_and_sort_strings}: ...
#' @rdname utilites
.find_and_sort_strings <- 
  function(strings, patterns){
    lapply(patterns,
           function(pattern){
             strings[grepl(pattern, strings, perl = T)]
           })
  }

#' @export maps
#' @aliases maps
#' @description \code{maps}: ...
#' @rdname utilites
maps <- function(data, value, from, to) {
  if (!is.list(value))
    value <- list(value)
  lapply(value,
         function(value) {
           data <- data[data[[from]] %in% value, ]
           vec <- data[[ to ]]
           names(vec) <- data[[ from ]]
           vec
         })
}

#' @export order_list
#' @aliases order_list
#' @description \code{order_list}: ...
#' @rdname utilites
order_list <- 
  function(
           list
           ){
    lt <- list()
    length(lt) <- length(list)
    names(lt) <- sort(names(list))
    for(i in names(lt)){
      lt[[i]] <- list[[i]]
    }
    return(lt)
  }

#' @export molconvert_structure
#' @aliases molconvert_structure
#' @description \code{molconvert_structure}: ...
#' @rdname utilites
## use 'molconvert' ...
## https://chemaxon.com/marvin
molconvert_structure <-
  function(smile, path){
    system(paste0("molconvert mol \"", smile, "\" -o ", path))
    src <- paste(readLines(path), collapse = "\n")
    ChemmineOB::convertToImage("MOL", "SVG", source = src, toFile = path)
    rsvg::rsvg_svg(path, path)
  }

#' @export obj.size
#' @aliases obj.size
#' @description \code{obj.size}: ...
#' @rdname utilites
obj.size <- function(x, ...) {
  format(object.size(x), units = "MB", ...)
}
