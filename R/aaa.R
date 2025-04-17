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

add_filename_suffix <- function(file, suffix, sep = "_") {
  if (is.null(file)) {
    stop('is.null(file), no `file` specified.')
  }
  name <- tools::file_path_sans_ext(file)
  ext <- tools::file_ext(file)
  if (ext == "gz") {
    ext <- paste0(tools::file_ext(name), ".", ext)
    name <- tools::file_path_sans_ext(name)
  }
  paste0(name, sep, suffix, ".", ext)
}

new_from_package <- function(Class, package, ...) {
  ClassDef <- getClass(Class, where = asNamespace(package))
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
  if (is.null(text) || !length(text) || (length(text) == 1 && !nchar(text))) {
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
.cairosvg_to_grob <- function(path){
    grImport2::pictureGrob(
      grImport2::readPicture(path),
      width = unit(1, "npc"),
      height = unit(1, "npc")
    )
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
    message(glue::glue("The file is: {file}"))
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
  c("Turquoise", "Blue", "Brown", "Yellow", "Green", "Red", "Black", "Pink",
    "Magenta", "Purple", "Greenyellow", "Tan", "Salmon", "Cyan", "Midnightblue",
    "Lightcyan", "Grey60", "Lightgreen", "Lightyellow", "Royalblue", "Darkred",
    "Darkgreen", "Darkturquoise", "Darkgrey", "Orange", "Darkorange", "White",
    "Skyblue", "Saddlebrown", "Steelblue", "Paleturquoise", "Violet",
    "Darkolivegreen", "Darkmagenta", "Sienna3", "Yellowgreen", "Skyblue3",
    "Plum1", "Orangered4", "Mediumpurple3", "Lightsteelblue1", "Lightcyan1",
    "Ivory", "Floralwhite", "Darkorange2", "Brown4", "Bisque4", "Darkslateblue",
    "Plum2", "Thistle2", "Thistle1", "Salmon4", "Palevioletred3", "Navajowhite2",
    "Maroon", "Lightpink4", "Lavenderblush3", "Honeydew1", "Darkseagreen4",
    "Coral1", "Antiquewhite4", "Coral2", "Mediumorchid", "Skyblue2", "Yellow4",
    "Skyblue1", "Plum", "Orangered3", "Mediumpurple2", "Lightsteelblue",
    "Lightcoral", "Indianred4", "Firebrick4", "Darkolivegreen4", "Brown2",
    "Blue2", "Darkviolet", "Plum3", "Thistle3", "Thistle", "Salmon2",
    "Palevioletred2", "Navajowhite1", "Magenta4", "Lightpink3", "Lavenderblush2",
    "Honeydew", "Darkseagreen3", "Coral", "Antiquewhite2", "Coral3",
    "Mediumpurple4", "Skyblue4", "Yellow3", "Sienna4", "Pink4", "Orangered1",
    "Mediumpurple1", "Lightslateblue", "Lightblue4", "Indianred3", "Firebrick3",
    "Darkolivegreen2", "Blueviolet", "Blue4", "Deeppink", "Plum4", "Thistle4",
    "Tan4", "Salmon1", "Palevioletred1", "Navajowhite", "Magenta3", "Lightpink2",
    "Lavenderblush1", "Green4", "Darkseagreen2", "Chocolate4", "Antiquewhite1",
    "Coral4", "Mistyrose", "Slateblue", "Yellow2", "Sienna2", "Pink3",
    "Orangered", "Mediumpurple", "Lightskyblue4", "Lightblue3", "Indianred2",
    "Firebrick2", "Darkolivegreen1", "Blue3", "Brown1", "Deeppink1",
    "Powderblue", "Tomato", "Tan3", "Royalblue3", "Palevioletred", "Moccasin",
    "Magenta2", "Lightpink1", "Lavenderblush", "Green3", "Darkseagreen1",
    "Chocolate3", "Aliceblue", "Cornflowerblue", "Navajowhite3", "Slateblue1",
    "Whitesmoke", "Sienna1", "Pink2", "Orange4", "Mediumorchid4",
    "Lightskyblue3", "Lightblue2", "Indianred1", "Firebrick", "Darkgoldenrod4",
    "Blue1", "Brown3", "Deeppink2", "Purple2", "Tomato2", "Tan2", "Royalblue2",
    "Paleturquoise4", "Mistyrose4", "Magenta1", "Lightpink", "Lavender",
    "Green2", "Darkseagreen", "Chocolate2", "Antiquewhite", "Cornsilk",
    "Navajowhite4", "Slateblue2", "Wheat3", "Sienna", "Pink1", "Orange3",
    "Mediumorchid3", "Lightskyblue2", "Lightblue1", "Indianred", "Dodgerblue4",
    "Darkgoldenrod3", "Blanchedalmond", "Burlywood", "Deepskyblue", "Red1",
    "Tomato4", "Tan1", "Rosybrown4", "Paleturquoise3", "Mistyrose3", "Linen",
    "Lightgoldenrodyellow", "Khaki4", "Green1", "Darksalmon", "Chocolate1",
    "Antiquewhite3", "Cornsilk2", "Oldlace", "Slateblue3", "Wheat1", "Seashell4",
    "Peru", "Orange2", "Mediumorchid2", "Lightskyblue1", "Lightblue", "Hotpink4",
    "Dodgerblue3", "Darkgoldenrod1", "Bisque3", "Burlywood1", "Deepskyblue4",
    "Red4", "Turquoise2", "Steelblue4", "Rosybrown3", "Paleturquoise1",
    "Mistyrose2", "Limegreen", "Lightgoldenrod4", "Khaki3", "Goldenrod4",
    "Darkorchid4", "Chocolate", "Aquamarine", "Cyan1", "Orange1", "Slateblue4",
    "Violetred4", "Seashell3", "Peachpuff4", "Olivedrab4", "Mediumorchid1",
    "Lightskyblue", "Lemonchiffon4", "Hotpink3", "Dodgerblue1", "Darkgoldenrod",
    "Bisque2", "Burlywood2", "Dodgerblue2", "Rosybrown2", "Turquoise4",
    "Steelblue3", "Rosybrown1", "Palegreen4", "Mistyrose1", "Lightyellow4",
    "Lightgoldenrod3", "Khaki2", "Goldenrod3", "Darkorchid3", "Chartreuse4",
    "Aquamarine1", "Cyan4", "Orangered2", "Snow", "Violetred2", "Seashell2",
    "Peachpuff3", "Olivedrab3", "Mediumblue", "Lightseagreen", "Lemonchiffon3",
    "Hotpink2", "Dodgerblue", "Darkblue", "Bisque1", "Burlywood3", "Firebrick1",
    "Royalblue1", "Violetred1", "Steelblue1", "Rosybrown", "Palegreen3",
    "Mintcream", "Lightyellow3", "Lightgoldenrod2", "Khaki1", "Goldenrod2",
    "Darkorchid2", "Chartreuse3", "Aquamarine2", "Darkcyan", "Orchid", "Snow2",
    "Violetred", "Seashell1", "Peachpuff2", "Olivedrab2", "Mediumaquamarine",
    "Lightsalmon4", "Lemonchiffon2", "Hotpink1", "Deepskyblue3", "Cyan3",
    "Bisque", "Burlywood4", "Forestgreen", "Royalblue4", "Violetred3",
    "Springgreen3", "Red3", "Palegreen1", "Mediumvioletred", "Lightyellow2",
    "Lightgoldenrod1", "Khaki", "Goldenrod1", "Darkorchid1", "Chartreuse2",
    "Aquamarine3", "Darkgoldenrod2", "Orchid1", "Snow4", "Turquoise3",
    "Seashell", "Peachpuff1", "Olivedrab1", "Maroon4", "Lightsalmon3",
    "Lemonchiffon1", "Hotpink", "Deepskyblue2", "Cyan2", "Beige", "Cadetblue",
    "Gainsboro", "Salmon3", "Wheat", "Springgreen2", "Red2", "Palegreen",
    "Mediumturquoise", "Lightyellow1", "Lightgoldenrod", "Ivory4", "Goldenrod",
    "Darkorchid", "Chartreuse1", "Aquamarine4", "Darkkhaki", "Orchid3",
    "Springgreen1", "Turquoise1", "Seagreen4", "Peachpuff", "Olivedrab",
    "Maroon3", "Lightsalmon2", "Lemonchiffon", "Honeydew4", "Deepskyblue1",
    "Cornsilk4", "Azure4", "Cadetblue1", "Ghostwhite", "Sandybrown", "Wheat2",
    "Springgreen", "Purple4", "Palegoldenrod", "Mediumspringgreen",
    "Lightsteelblue4", "Lightcyan4", "Ivory3", "Gold3", "Darkorange4",
    "Chartreuse", "Azure", "Darkolivegreen3", "Palegreen2", "Springgreen4",
    "Tomato3", "Seagreen3", "Papayawhip", "Navyblue", "Maroon2", "Lightsalmon1",
    "Lawngreen", "Honeydew3", "Deeppink4", "Cornsilk3", "Azure3", "Cadetblue2",
    "Gold", "Seagreen", "Wheat4", "Snow3", "Purple3", "Orchid4",
    "Mediumslateblue", "Lightsteelblue3", "Lightcyan3", "Ivory2", "Gold2",
    "Darkorange3", "Cadetblue4", "Azure1", "Darkorange1", "Paleturquoise2",
    "Steelblue2", "Tomato1", "Seagreen2", "Palevioletred4", "Navy", "Maroon1",
    "Lightsalmon", "Lavenderblush4", "Honeydew2", "Deeppink3", "Cornsilk1",
    "Azure2", "Cadetblue3", "Gold4", "Seagreen1", "Yellow1", "Snow1", "Purple1",
    "Orchid2", "Mediumseagreen", "Lightsteelblue2", "Lightcyan2", "Ivory1",
    "Gold1")
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
  values = NULL, rev = TRUE)
{
  category <- match.arg(category)
  if (!is.null(values)) {
    from <- -max(ceiling(abs(range(values))))
    to <- -from
    if (category == "seq") {
      from <- min(values)
    }
  }
  if (sample) {
    pals <- dplyr::filter(RColorBrewer::brewer.pal.info, category == !!category)
    which <- sample(seq_len(nrow(pals)), 1)
    pal <- rownames(pals)[ which ]
    max <- pals$maxcolors[ which ]
  } else {
    pal <- "RdBu"
    max <- 11L
  }
  colors <- RColorBrewer::brewer.pal(max, pal)
  if (rev) {
    colors <- rev(colors)
  }
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

.forest_lines <- function(data, header = TRUE) {
  if (header) {
    end <- nrow(data) + 2
  } else {
    end <- nrow(data) + 1
  }
  lst <- rep(list(gpar(lty = 1, lwd = 2, col = "black")), 3)
  names(lst) <- c(1:2, end)
  return(lst)
}

rbind_list <- function(lst, .id = NULL, keep_missing = TRUE) {
  if (keep_missing) {
    data <- lapply(lst,
      function(x) {
        if (!nrow(x)) {
          dplyr::add_row(x)
        } else x
      })
  } else {
    data <- lst
  }
  data <- do.call(
    dplyr::bind_rows, c(data, .id = .id)
  )
  data
}

ors <- function(lst) {
  .operator_for_list(lst, `|`)
}

ands <- function(lst) {
  .operator_for_list(lst, `&`)
}

.operator_for_list <- function(lst, fun) {
  if (!length(lst)) {
    stop('!length(lst).')
  }
  res <- lst[[1]]
  for (i in lst[-1]) {
    res <- fun(res, i)
  }
  res
}

.show_layers <- function(layers) {
  cat(crayon::silver("layers of", length(layers), "\n"))
  mapply(names(layers), seq_along(layers),
    FUN = function(name, seq){
      cat(crayon::silver("  +++ layer", seq, ": ", name, "+++\n"))
      print(layers[[ seq ]])
      cat("\n")
    })
  cat("\n")
}

calculate_layout <- function(n, max_ratio = 3, prefer_rows = FALSE) {
  if (n < 1) {
    stop('n < 1.')
  }
  if (max_ratio < 1) {
    stop('max_ratio < 1.')
  }

  # Generate candidate layout schemes
  generate_candidates <- function(n) {
    candidates <- list()
    max_dim <- ceiling(sqrt(n) * max_ratio)
    
    for (r in 1:max_dim) {
      c <- ceiling(n / r)
      ratio <- max(r, c) / min(r, c)
      
      if (ratio <= max_ratio && r * c >= n) {
        candidates[[length(candidates)+1]] <- c(r, c, r*c - n, abs(r - c))
      }
    }
    do.call(rbind, candidates)
  }

  # Scoring criteria: Blank quantity>Column difference>Aspect ratio
  evaluate_candidates <- function(candidates) {
    candidates[order(
      candidates[,3],
      candidates[,4],
      candidates[,2]/candidates[,1]
    ), ]
  }

  # Execute the calculation process
  candidates <- generate_candidates(n)
  if (nrow(candidates) == 0) {
    stop('nrow(candidates) == 0.')
  }
  
  sorted <- evaluate_candidates(candidates)
  best <- sorted[1, 1:2]
  
  # Processing direction preference
  if (!prefer_rows && best[1] > best[2]) best <- rev(best)
  
  setNames(as.integer(best), c("rows", "cols"))
}

calculate_layout_size <- function(layout, base_size = 5, 
  decay_rate = 0.7, min_size = 1)
{
  rows <- layout[[1]]
  cols <- layout[[2]]
  layout_complexity <- sqrt(rows + cols)

  # Calculation formula for dynamic attenuation coefficient
  size <- base_size * (decay_rate ^ (layout_complexity - 1))

  # Apply size limitations and maintain consistent width and height
  final_size <- max(min_size, min(base_size, size))

  list(
    dimension = c(rows, cols),
    size = rep(round(final_size, 1), 2),
    area_ratio = (final_size^2) / (base_size^2)
  )
}

formal_name <- function(x) {
  gs(make.names(x), "[.]+", "_")
}

try_get_url <- function(url, maxTry = 10, pattern = "0A000126:SSL routines") {
  n <- 0L
  res <- NULL
  while (is.null(res) || (n <= maxTry && inherits(res, "try-error") && grpl(res, pattern)))
  {
    n <- n + 1L
    res <- try(RCurl::getURL(url), TRUE)
    if (!is.null(res)) {
      message(glue::glue("\nTry again ({n}): {url}\n"))
    }
  }
  return(res)
}

pngDpi <- function(filename, width = 7, height = 7, dpi = 300) {
  png(filename, width = width, height = height, units = "in", res = dpi)
}
