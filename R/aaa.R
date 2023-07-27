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

get_signature <- 
  function(args){
    vapply(args, function(arg) class(arg)[1], FUN.VALUE = "ch")
  }

.fresh_param <- function(default, args){
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
  message("Export `.expath` for inst files usage.")
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

clearMatch <- function(strs)
{
  strs <- unlist(strs)
  strs <- strs[ !vapply(strs, is.na, logical(1)) ]
  strs
}

eg <- data.frame(x = sample(1:10, 10), y = sample(1:10, 10), z = sample(1:10, 10))
eg2 <- data.frame(a = sample(1:10, 10), b = sample(1:10, 10), c = sample(1:10, 10))

textSh <- function(..., sep = "", exdent = 4, ending = "\n",
           pre_collapse = F, collapse = "\n",
           pre_trunc = T, trunc_width = 200,
           pre_wrap = T, wrap_width = 60){
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
    writeLines(gsub("(?<=\n)|(?<=^)", exdent, text, perl = T))
    if (!is.null(ending))
      cat(ending)
  }

.text_fold <- function(text, width = 200, ellipsis = crayon::silver("...(fold)")){
  stringr::str_trunc(text, width = width, ellipsis = ellipsis)
}

gett <- function(obj){
    system(paste("echo", obj, "| xsel -b -i"))
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

color_set <- function() {
  c("#9EDAE5FF", "#DBDB8DFF", "#F7B6D2FF", "#C49C94FF", "#C5B0D5FF", "#FF9896FF",
    "#98DF8AFF", "#FFBB78FF", "#AEC7E8FF", "#17BECFFF", "#BCBD22FF", "#7F7F7FFF",
    "#E377C2FF", "#8C564BFF", "#9467BDFF", "#D62728FF", "#2CA02CFF", "#FF7F0EFF",
    "#1F77B4FF", "#FED439FF", "#709AE1FF", "#D2AF81FF", "#FD7446FF", "#D5E4A2FF",
    "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF", "#370335FF", "#075149FF",
    "#C80813FF", "#91331FFF", "#1A9993FF", "#FD8CC1FF", "#FF0000FF", "#FF9900FF",
    "#FFCC00FF", "#00FF00FF", "#6699FFFF", "#CC33FFFF")
}

