# ==========================================================================
# utilites
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

binary_search <- function(fun = function(x) x <= 5, low = 1, high = 12)
{
  n <- 0L
  while (low < high) {
    n <- n + 1L
    mayThat <- floor((low + high) / 2)
    has <- fun(mayThat)
    if (has) {
      low <- mayThat + 1
    } else {
      high <- mayThat
    }
  }
  message("Loop times: ", n)
  return(low - 1)
}

new_from_package <- function(Class, package, ...) {
  ClassDef <- getClass("EList", where = asNamespace(package))
  value <- .Call(methods:::C_new_object, ClassDef)
  initialize(value, ...)
}

setFakeClasses <- function(classes) {
  lapply(classes,
    function(name) {
      if (inherits(try(getClass(name, where = .GlobalEnv), silent = TRUE), "try-error")) {
        suppressWarnings(setClass(name, where = topenv()))
      }
    })
}

split_text_by_width <- function(text, width) {
  if (is.null(text) || !length(text)) {
    return(text)
  }
  parts <- unlist(strsplit(text, "(?<=\\p{Han})", perl = TRUE))
  # parts <- unlist(strsplit(parts, "(?=\\p{Han})", perl = TRUE))
  parts <- unlist(strsplit(parts, "(?<=\\s)", perl = TRUE))
  parts <- unlist(strsplit(parts, "(?=\\s)", perl = TRUE))
  current_width <- 0
  current_line <- ""
  lines <- c()
  for (part in parts) {
    part <- trimws(part)
    if (part == "") next
    part_width <- sum(nchar(stringr::str_replace_all(part, "[\u4e00-\u9fa5]", "")) + 
                      2 * stringr::str_count(part, "[\u4e00-\u9fa5]"))
    if (current_width + part_width > width) {
      lines <- c(lines, current_line)
      current_line <- part
      current_width <- part_width
    } else {
      if (current_line != "") {
        current_line <- paste0(current_line, part)
      } else {
        current_line <- part
      }
      current_width <- current_width + part_width
    }
  }
  if (nchar(current_line) > 0) {
    lines <- c(lines, current_line)
  }
  return(lines)
}

#' @aliases utilites
#'
#' @title utilites for programming
#'
#' @description This is a combination of tools that are not always used.
#'
#' @name utilites
NULL
#> NULL

reCallMethod <- function(funName, args, ...){
    arg.order <- unname(getGeneric(funName)@signature)
    args.missing <- !arg.order %in% names(args)
    if (any(args.missing)) {
      args.missing <- arg.order[args.missing]
      args.missing <- sapply(args.missing, simplify = FALSE,
                             function(x) structure(0L, class = "missing"))
      args <- c(args, args.missing)
    }
    args <- lapply(arg.order, function(i) args[[i]])
    sig <- get_signature(args)
    method <- selectMethod(funName, sig)
    last_fun <- sys.function(sys.parent())
    n <- 0
    while (identical(last_fun, method@.Data, ignore.environment = TRUE)) {
      if (n == 0) {
        mlist <- getMethodsForDispatch(getGeneric(funName))
      }
      n <- n + 1
      rm(list = paste0(method@defined, collapse = "#"), envir = mlist)
      method <- selectMethod(funName, sig, mlist = mlist)
    }
    expr <- paste0("method@.Data(",
                   paste0(paste0(arg.order, " = args[[",
                                 seq_along(arg.order), "]]"),
                          collapse = ", "),
                   ", ...)")
    eval(parse(text = expr))
  }

get_signature <- function(args){
  vapply(args, function(arg) class(arg)[1], FUN.VALUE = "ch")
}

.fresh_param <- function(default, args){
    if (missing(args))
      args <- as.list(parent.frame())
    args <- args[ !vapply(args, is.name, TRUE) ]
    sapply(unique(c(names(default), names(args))),
           simplify = FALSE,
           function(name){
             if (any(name == names(args)))
               args[[ name ]]
             else
               default[[ name ]]
           })
  }

.fresh_param2 <- function(default, args){
    if (missing(args))
      return(default)
    if (length(args) == 0)
      return(default)
    .fresh_param(default, args)
  }

.fresh_param2f <- function(default, args, class = "gpar"){
    structure(.fresh_param2(default, args), class = class)
  }

.check_columns <- function(obj, lst, tip = "`obj`"){
    if (!is.data.frame(obj))
      stop(paste0("'", tip, "' must be a 'data.frame'."))
    lapply(lst, function(col){
             if (is.null(obj[[ col ]]))
               stop(paste0("'", tip, "' must contains a column of '", col, "'."))
           })
  }

.message_info_viewport <- function(info = "info"){
    .message_info(info, "current.viewport:",
                  paste0("\n\t", paste0(grid::current.viewport())))
  }

.message_info <- function(main, sub, arg = NULL, sig = "##"){
    message(sig, " ", main, ": ", sub, " ", arg)
  }

.suggest_bio_package <- function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE))
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
  file.path(tempdir(), "tmp_pdf.pdf")
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

#' @export fill_list
#' @aliases fill_list
#' @description \code{fill_list}: ...
#' @rdname utilites
fill_list <- function(names, vec, default = vec[1]) {
  .as_dic(vec, names, default, fill = TRUE, as.list = FALSE, na.rm = FALSE)
}

#' @export n
#' @aliases n
#' @description \code{n}: ...
#' @rdname utilites
n <- function(name, n){
  if (n == 0) {
    return(NULL)
  }
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
  names(lst) <- getNamesOfCall(call)
  lst
}

getNamesOfCall <- function(call) {
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
  return(names)
}

#' @export repSuffix
#' @aliases repSuffix
#' @description \code{repSuffix}: ...
#' @rdname utilites
repSuffix <- 
  function(chs, anno = ".rep."){
    gsub(paste0(anno, 1, "$"), "",
         vapply(seq_along(chs), FUN.VALUE = character(1),
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
  if (requireNamespace("rsvg", quietly = TRUE)) {
    .check_external_svg()
  }
}

#' @export agroup
#' @aliases agroup
#' @description \code{agroup}: ...
#' @rdname utilites
agroup <- function(group, value, FUN.VALUE = character(1)) {
  ug <- unique(group)
  if (length(ug) > length(value))
    stop( "the length of 'value' not enough to assign" )
  dic <- .as_dic(value, ug, fill = FALSE, na.rm = FALSE)
  vapply(group, function(g) dic[[g]], FUN.VALUE)
}

#' @export write_tsv
#' @aliases write_tsv
#' @description \code{write_tsv}: ...
#' @rdname utilites
write_tsv <-
  function(x, filename, col.names = TRUE, row.names = FALSE){
    write.table(x, file = filename, sep = "\t",
                col.names = col.names, row.names = row.names, quote = FALSE)
  }

#' @export read_tsv
#' @aliases read_tsv
#' @description \code{read_tsv}: ...
#' @rdname utilites
read_tsv <- function(path){
  file <- data.table::fread(input = path, sep = "\t",
                            header = TRUE, quote = "", check.names = FALSE)
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
           fixed = FALSE
           ){
    envir <- environment()
    mapply(mutate_set, replace_set,
           MoreArgs = list(envir = envir, fixed = fixed),
           FUN = function(mutate, replace, envir,
                          fixed = FALSE, names = get("names", envir = envir)){
             names <- gsub(mutate, replace, names, perl = ifelse(fixed, FALSE, TRUE), fixed = fixed)
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
    stop( "is.character(data[[ by ]]) == FALSE" )
  meta <- unlist(meta.lst)
  names(meta) <- rep(names(meta.lst), lengths(meta.lst))
  meta <- as.list(turn_vector(meta))
  data <- data[data[[by]] %in% names(meta), ]
  group <- data.frame(order = seq_along(data[[ by ]]), col = data[[ by ]])
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

group_strings <- function(strings, patterns, target = NA){
    if (is.null(names(patterns)))
      stop("`patterns` must be characters with names.")
    lst <- .find_and_sort_strings(strings, patterns)
    lst <- lapply(names(lst), function(name){
                    data.frame(target = lst[[name]], group = name)
           })
    df <- do.call(rbind, lst)
    if (!is.na(target)) {
      colnames(df)[1] <- target
    }
    tibble::as_tibble(df)
  }

#' @export .find_and_sort_strings
#' @aliases .find_and_sort_strings
#' @description \code{.find_and_sort_strings}: ...
#' @rdname utilites
.find_and_sort_strings <- function(strings, patterns, ...)
{
  sapply(patterns, simplify = FALSE,
    function(pattern){
      strings[grepl(pattern, strings, ...)]
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
order_list <- function(list) {
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
molconvert_structure <- function(smile, path)
{
  system(paste0("molconvert mol \"", smile, "\" -o ", path))
  src <- paste(readLines(path), collapse = "\n")
  ChemmineOB::convertToImage("MOL", "SVG", source = src, toFile = path)
  rsvg::rsvg_svg(path, path)
}

clearMatch <- function(strs)
{
  strs <- unlist(strs)
  strs <- strs[ !vapply(strs, is.na, logical(1)) ]
  strs
}

eg <- data.frame(x = sample(1:10, 10), y = sample(1:10, 10), z = sample(1:10, 10))
eg2 <- data.frame(a = sample(1:10, 10), b = sample(1:10, 10), c = sample(1:10, 10))

textSh <- function(..., sep = "", exdent = 4, ending = "\n",
           pre_collapse = FALSE, collapse = "\n",
           pre_trunc = TRUE, trunc_width = 200,
           pre_wrap = TRUE, wrap_width = 60){
    text <- list(...)
    if (pre_collapse) {
      text <- vapply(text, paste, "ch", collapse = collapse)
    }
    text <- paste(text, sep = sep)
    if (pre_trunc) {
      text <- .text_fold(text, trunc_width)
    }
    if (pre_wrap) {
      text <- paste0(strwrap(text, width = wrap_width), collapse = "\n")
    }
    exdent <- paste0(rep(" ", exdent), collapse = "")
    writeLines(gsub("(?<=\n)|(?<=^)", exdent, text, perl = TRUE))
    if (!is.null(ending))
      cat(ending)
  }

.text_fold <- function(text, width = 200, ellipsis = crayon::silver("...(fold)")){
  stringr::str_trunc(text, width = width, ellipsis = ellipsis)
}

# ==========================================================================
# get external grob
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @export .expathsvg
.expathsvg <- function() {
  .expathsvg <- system.file("extdata", "svg", package = gsub("^.*:", "",
      environmentName(topenv())))
  assign('.expathsvg', .expathsvg, envir = topenv())
}

prefix <- c()

#' @export .check_external_svg
.check_external_svg <- function(){
  files <- list.files(.expathsvg, "\\.svg$", full.names = TRUE)
  log.path <- file.path(.expathsvg, "log")
  if (file.exists(log.path)) {
    log <- readLines(log.path)
    log <- log[vapply(log, file.exists, TRUE)]
  } else {
    log <- c()
  }
  if (length(files) > 0) {
    new.log <-
      lapply(files,
        function(file) {
          if (!file %in% log) {
            rsvg::rsvg_svg(file, file)
            file
          }
        })
    new.log <- unlist(new.log)
    log <- c(log, new.log)
  }
  if (!is.null(log))
    writeLines(log, log.path)
}

#' @export ex_grob
ex_grob <- function(name, fun = .cairosvg_to_grob){
  file <- file.path(.expathsvg, paste0(name, ".svg"))
  if (file.exists(file)) {
    fun(file)
  } else {
    stop("file.exsits(file) == FALSE")
  }
}

#' @export ex_pic
ex_pic <- function(name){
  ex_grob(name, fun = grImport2::readPicture)
}


.show <- function(object){
    cat(class(object), "\n")
    slots_mapply(object, function(names, slots){
              cat(names, ":\n", sep = "")
              cat(str(slots))
              cat("\n\n")
           })
  }

slots_mapply <- function(x, fun, ...){
    slots <- attributes(x)
    slots <- slots[-length(slots)]
    res <- mapply(fun, slot = slots, name = names(slots), ...)
    return(res)
  }

spacing <- space <- function(env = .GlobalEnv, all.names = TRUE) {
  names <- ls(envir = env, all.names = all.names)
  info <- vapply(names, FUN.VALUE = character(1),
    function(name) {
      obj.size(get(name, envir = env))
    })
  info <- tibble::tibble(name = names(info), size = unname(info))
  info <- dplyr::mutate(info, value = as.double(stringr::str_extract(size, "^[0-9.]+")))
  info <- dplyr::arrange(info, dplyr::desc(value))
  info
}

save_small <- function(cutoff = 50, file = "small.rdata") {
  space <- space()
  space <- dplyr::filter(space, .data$value < !!cutoff)
  save(list = space$name, file = file, envir = .GlobalEnv)
}

# ==========================================================================
# colors
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

wgcna_colors <- function() {
  c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink",
    "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue",
    "lightcyan", "grey60", "lightgreen", "lightyellow", "royalblue", "darkred",
    "darkgreen", "darkturquoise", "darkgrey", "orange", "darkorange", "white",
    "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet",
    "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen", "skyblue3",
    "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1",
    "ivory", "floralwhite", "darkorange2", "brown4", "bisque4", "darkslateblue",
    "plum2", "thistle2", "thistle1", "salmon4", "palevioletred3", "navajowhite2",
    "maroon", "lightpink4", "lavenderblush3", "honeydew1", "darkseagreen4",
    "coral1", "antiquewhite4", "coral2", "mediumorchid", "skyblue2", "yellow4",
    "skyblue1", "plum", "orangered3", "mediumpurple2", "lightsteelblue",
    "lightcoral", "indianred4", "firebrick4", "darkolivegreen4", "brown2",
    "blue2", "darkviolet", "plum3", "thistle3", "thistle", "salmon2",
    "palevioletred2", "navajowhite1", "magenta4", "lightpink3", "lavenderblush2",
    "honeydew", "darkseagreen3", "coral", "antiquewhite2", "coral3",
    "mediumpurple4", "skyblue4", "yellow3", "sienna4", "pink4", "orangered1",
    "mediumpurple1", "lightslateblue", "lightblue4", "indianred3", "firebrick3",
    "darkolivegreen2", "blueviolet", "blue4", "deeppink", "plum4", "thistle4",
    "tan4", "salmon1", "palevioletred1", "navajowhite", "magenta3", "lightpink2",
    "lavenderblush1", "green4", "darkseagreen2", "chocolate4", "antiquewhite1",
    "coral4", "mistyrose", "slateblue", "yellow2", "sienna2", "pink3",
    "orangered", "mediumpurple", "lightskyblue4", "lightblue3", "indianred2",
    "firebrick2", "darkolivegreen1", "blue3", "brown1", "deeppink1",
    "powderblue", "tomato", "tan3", "royalblue3", "palevioletred", "moccasin",
    "magenta2", "lightpink1", "lavenderblush", "green3", "darkseagreen1",
    "chocolate3", "aliceblue", "cornflowerblue", "navajowhite3", "slateblue1",
    "whitesmoke", "sienna1", "pink2", "orange4", "mediumorchid4",
    "lightskyblue3", "lightblue2", "indianred1", "firebrick", "darkgoldenrod4",
    "blue1", "brown3", "deeppink2", "purple2", "tomato2", "tan2", "royalblue2",
    "paleturquoise4", "mistyrose4", "magenta1", "lightpink", "lavender",
    "green2", "darkseagreen", "chocolate2", "antiquewhite", "cornsilk",
    "navajowhite4", "slateblue2", "wheat3", "sienna", "pink1", "orange3",
    "mediumorchid3", "lightskyblue2", "lightblue1", "indianred", "dodgerblue4",
    "darkgoldenrod3", "blanchedalmond", "burlywood", "deepskyblue", "red1",
    "tomato4", "tan1", "rosybrown4", "paleturquoise3", "mistyrose3", "linen",
    "lightgoldenrodyellow", "khaki4", "green1", "darksalmon", "chocolate1",
    "antiquewhite3", "cornsilk2", "oldlace", "slateblue3", "wheat1", "seashell4",
    "peru", "orange2", "mediumorchid2", "lightskyblue1", "lightblue", "hotpink4",
    "dodgerblue3", "darkgoldenrod1", "bisque3", "burlywood1", "deepskyblue4",
    "red4", "turquoise2", "steelblue4", "rosybrown3", "paleturquoise1",
    "mistyrose2", "limegreen", "lightgoldenrod4", "khaki3", "goldenrod4",
    "darkorchid4", "chocolate", "aquamarine", "cyan1", "orange1", "slateblue4",
    "violetred4", "seashell3", "peachpuff4", "olivedrab4", "mediumorchid1",
    "lightskyblue", "lemonchiffon4", "hotpink3", "dodgerblue1", "darkgoldenrod",
    "bisque2", "burlywood2", "dodgerblue2", "rosybrown2", "turquoise4",
    "steelblue3", "rosybrown1", "palegreen4", "mistyrose1", "lightyellow4",
    "lightgoldenrod3", "khaki2", "goldenrod3", "darkorchid3", "chartreuse4",
    "aquamarine1", "cyan4", "orangered2", "snow", "violetred2", "seashell2",
    "peachpuff3", "olivedrab3", "mediumblue", "lightseagreen", "lemonchiffon3",
    "hotpink2", "dodgerblue", "darkblue", "bisque1", "burlywood3", "firebrick1",
    "royalblue1", "violetred1", "steelblue1", "rosybrown", "palegreen3",
    "mintcream", "lightyellow3", "lightgoldenrod2", "khaki1", "goldenrod2",
    "darkorchid2", "chartreuse3", "aquamarine2", "darkcyan", "orchid", "snow2",
    "violetred", "seashell1", "peachpuff2", "olivedrab2", "mediumaquamarine",
    "lightsalmon4", "lemonchiffon2", "hotpink1", "deepskyblue3", "cyan3",
    "bisque", "burlywood4", "forestgreen", "royalblue4", "violetred3",
    "springgreen3", "red3", "palegreen1", "mediumvioletred", "lightyellow2",
    "lightgoldenrod1", "khaki", "goldenrod1", "darkorchid1", "chartreuse2",
    "aquamarine3", "darkgoldenrod2", "orchid1", "snow4", "turquoise3",
    "seashell", "peachpuff1", "olivedrab1", "maroon4", "lightsalmon3",
    "lemonchiffon1", "hotpink", "deepskyblue2", "cyan2", "beige", "cadetblue",
    "gainsboro", "salmon3", "wheat", "springgreen2", "red2", "palegreen",
    "mediumturquoise", "lightyellow1", "lightgoldenrod", "ivory4", "goldenrod",
    "darkorchid", "chartreuse1", "aquamarine4", "darkkhaki", "orchid3",
    "springgreen1", "turquoise1", "seagreen4", "peachpuff", "olivedrab",
    "maroon3", "lightsalmon2", "lemonchiffon", "honeydew4", "deepskyblue1",
    "cornsilk4", "azure4", "cadetblue1", "ghostwhite", "sandybrown", "wheat2",
    "springgreen", "purple4", "palegoldenrod", "mediumspringgreen",
    "lightsteelblue4", "lightcyan4", "ivory3", "gold3", "darkorange4",
    "chartreuse", "azure", "darkolivegreen3", "palegreen2", "springgreen4",
    "tomato3", "seagreen3", "papayawhip", "navyblue", "maroon2", "lightsalmon1",
    "lawngreen", "honeydew3", "deeppink4", "cornsilk3", "azure3", "cadetblue2",
    "gold", "seagreen", "wheat4", "snow3", "purple3", "orchid4",
    "mediumslateblue", "lightsteelblue3", "lightcyan3", "ivory2", "gold2",
    "darkorange3", "cadetblue4", "azure1", "darkorange1", "paleturquoise2",
    "steelblue2", "tomato1", "seagreen2", "palevioletred4", "navy", "maroon1",
    "lightsalmon", "lavenderblush4", "honeydew2", "deeppink3", "cornsilk1",
    "azure2", "cadetblue3", "gold4", "seagreen1", "yellow1", "snow1", "purple1",
    "orchid2", "mediumseagreen", "lightsteelblue2", "lightcyan2", "ivory1",
    "gold1")
}

color_set <- function(more = FALSE) {
  shiny <- c("#9EDAE5FF", "#DBDB8DFF", "#F7B6D2FF", "#C49C94FF", "#C5B0D5FF", "#FF9896FF",
    "#98DF8AFF", "#FFBB78FF", "#AEC7E8FF", "#17BECFFF", "#BCBD22FF", "#7F7F7FFF",
    "#E377C2FF", "#8C564BFF", "#9467BDFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF",
    "#1F77B4FF", "#FED439FF", "#709AE1FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF",
    "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF", "#370335FF", "#075149FF",
    "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF", "#FF0000FF", "#FF9900FF",
    "#FFCC00FF", "#00FF00FF", "#6699FFFF", "#CC33FFFF")
  if (more) {
    cols <- c(shiny, wgcna_colors())
    cols <- t(col2rgb(cols))
    cols <- cols[ !duplicated(cols), ]
    apply(cols, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  } else {
    shiny
  }
}

color_gradient <- function() {
  c("#67001FFF", "#B2182BFF", "#D6604DFF",
    "#F4A582FF", "#FDDBC7FF", "#F7F7F7FF",
    "#D1E5F0FF", "#92C5DEFF", "#4393C3FF",
    "#2166ACFF", "#053061FF")
}

fun_color <- function(from = -1, to = 1, sample = TRUE, category = c("div", "seq"),
  values = NULL)
{
  if (!is.null(values)) {
    from <- -max(ceiling(abs(range(values))))
    to <- -from
  }
  if (sample) {
    category <- match.arg(category)
    pals <- dplyr::filter(RColorBrewer::brewer.pal.info, category == !!category)
    which <- sample(seq_len(nrow(pals)), 1)
    pal <- rownames(pals)[ which ]
    max <- pals$maxcolors[ which ]
  } else {
    pal <- "RdBu"
    max <- 11L
  }
  colors <- rev(RColorBrewer::brewer.pal(max, pal))
  circlize::colorRamp2(seq(from, to, length.out = max), colors)
}

draw_smile <- function(smile, file, pdf = TRUE) {
  file <- .smiles_to_cairosvg(smile, file, pdf)
  .file_fig(file)
}

#' @importFrom ChemmineOB convertToImage
#' @importFrom rsvg rsvg_svg
.smiles_to_cairosvg <- function(smile, path, pdf = FALSE){
  ChemmineOB::convertToImage("SMI", "SVG", source = smile, toFile = path)
  if (pdf) {
    rsvg::rsvg_pdf(path, file <- gs(path, "\\.svg$", ".pdf"))
    return(file)
  } else {
    rsvg::rsvg_svg(path, path)
    return(path)
  }
}


.get_legend <- function(p){
  obj <- cowplot::get_plot_component(p, "guide-box", TRUE)
  obj[ vapply(obj, function(x) is(x, "gtable"), logical(1)) ][[ 1 ]]
  # p <- ggplot2:::ggplot_build.ggplot(p)$plot
  # theme <- ggplot2:::plot_theme(p)
  # position <- theme$legend.position
  # ggplot2:::build_guides(p$scales, p$layers, p$mapping,
  # position, theme, p$guides, p$labels)
}

.pattern.cn <- function() {
  "[\u4e00-\u9fa5]"
}

replaceFunInPackage <- function(name, fun, package) {
  nameSpace <- asNamespace(package)
  unlockBinding(name, nameSpace)
  assign(name, fun, envir = nameSpace)
  lockBinding(name, nameSpace)
}

get_grid <- function(len, ncols = 2, nrows = NULL) {
  if (is.null(nrows)) {
    nrows <- len %/% ncols
  } else if (is.null(ncols)) {
    ncols <- len %/% nrows
  }
  i <- seq(len)
  col <- ((i - 1) %% ncols) + 1
  row <- (((i) - 1) %/% ncols) + 1
  namel(row, col)
}

reset_object_name <- function(old, new, env = .GlobalEnv) {
  obj <- get(old, envir = env)
  assign(new, obj, envir = env)
  rm(list = old, envir = env)
  message(glue::glue("Rename object name: {old} -> {new}"))
}

bind <- function(..., co = ", ", quote = FALSE) {
  if (quote) {
    paste0('"', ..., '"', collapse = co)
  } else {
    paste0(..., collapse = co)
  }
}

collate_common_function <- function(package, freq = 10) {
  codes <- lapply(list.files(package, "\\.R$", full.names = TRUE), readLines)
  codes <- stringr::str_extract_all(unlist(codes), "(?<!=::)[a-zA-Z0-9_.]*(?=\\()")
  codes <- table(unlist(codes))
  codes <- names(codes[ codes > freq ])
  codes <- codes[ nchar(codes) >= 3 ]
  return(codes)
}
