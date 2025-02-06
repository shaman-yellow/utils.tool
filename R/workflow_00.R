# ==========================================================================
# for wrap various of analysis workflow contains many steps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setClass("virtual_job", "VIRTUAL")

setClassUnion("function_or_NULL", c("function", "NULL"))

.expect_col <- setClass("expect_col",
  representation = representation(
    name = "character", pattern_find = "character",
    fun_check = "function_or_NULL",
    pattern_recode = "character", fun_mutate = "function_or_NULL"),
  prototype = prototype(
    name = character(1), pattern_find = character(1),
    fun_check = NULL, pattern_recode = character(0), fun_mutate = NULL
    ))

gcol <- general_col <- function(name, pattern, check = NULL, 
  recode = character(0), mutate = NULL)
{
  .expect_col(
    name = name, pattern_find = pattern, fun_check = check, 
    pattern_recode = recode, fun_mutate = mutate
  )
}

gcols <- general_cols <- function(..., global = NULL, 
  uniqueness = character(0), db_file = character(0))
{
  .expect_cols(
    list(...), global = global,
    uniqueness = uniqueness, db_file = db_file
  )
}

setMethod("show", signature = c(object = "expect_col"),
  function(object){
    content <- c(
      paste(crayon::yellow("Name:"), object@name), 
      paste(crayon::silver("Pattern:"), bind(object@pattern_find))
    )
    if (length(object@pattern_recode)) {
      extra <- paste(
        crayon::silver("Recode:"), 
        bind(names(object@pattern_recode), " = ", unname(object@pattern_recode))
      )
      content <- append(content, extra)
    }
    extra <- paste(
      crayon::silver("Mutate:"), !is.null(object@fun_mutate)
    )
    content <- append(content, extra)
    writeLines(content)
  })

.expect_cols <- setClass("expect_cols",
  contains = c("list"),
  representation = representation(
    global = "function_or_NULL", uniqueness = "character",
    db_file = "character"
  ),
  prototype = prototype(global = NULL))

setMethod("show", signature = c(object = "expect_cols"),
  function(object){
    if (length(object@uniqueness)) {
      writeLines(
        c(paste(crayon::blue("Uniqueness:"), object@uniqueness), 
          paste(crayon::blue("DB file:"), object@db_file),
          crayon::silver("++++++"))
      )
    }
    invisible(lapply(object@.Data, show))
  })

setValidity("expect_cols",
  function(object){
    all(vapply(object@.Data, is, logical(1), "expect_col"))
  })

.meth <- setClass("meth",
  contains = c("environment"),
  representation = representation(init = "logical"),
  prototype = prototype(TRUE, init = TRUE))

setMethod("$", signature = c(x = "meth"),
  function(x, name){
    x@.xData[[ name ]]
  })

setMethod("$<-", signature = c(x = "meth"),
  function(x, name, value){
    if (!is(value, "character")) {
      stop('!is(value, "character"), must be character.')
    }
    x@.xData[[ name ]] <- value
    if (!any(x@.xData$.order == name)) {
      x@.xData$.order <- append(x@.xData$.order, name)
    }
    return(x)
  })

setMethod("[[", signature = c(x = "meth"),
  function(x, i, ...){
    if (is(i, "numeric")) {
      x@.xData[[ x@.xData$.order[i] ]]
    } else {
      x@.xData[[ i ]]
    }
  })

setMethod("[[<-", signature = c(x = "meth"),
  function(x, i, ..., value){
    if (!is(i, "character")) {
      stop('!is(i, "character"), `i` must be "character".')
    }
    if (is(value, "NULL")) {
      x@.xData$.order <- x@.xData$.order[ x@.xData$.order != i ]
      x@.xData[[ i ]] <- NULL
      return(x)
    }
    if (!is(value, "character")) {
      stop('!is(value, "character"), must be character.')
    }
    if (!any(x@.xData$.order == i)) {
      x@.xData$.order <- append(x@.xData$.order, i)
    }
    x@.xData[[ i ]] <- value
    return(x)
  })

setMethod("initialize", "meth",
  function(.Object, ...) {
    .Object <- callNextMethod()
    if (is.null(.Object@.xData)) {
      .Object@.xData <- new.env()
      .Object@.xData$.order <- character(0)
    }
    .Object
  })

setMethod("length", signature = c(x = "meth"),
  function(x){
    length(x@.xData$.order)
  })

.snap <- setClass("snap", contains = c("meth"))

setClassUnion("meth_or_NULL", c("meth", "NULL"))

#' @exportClass job
.job <- setClass("job",
  contains = c("virtual_job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY",
    step = "integer",
    info = "ANY",
    pg = "character",
    cite = "character",
    method = "character",
    meth = "meth_or_NULL",
    snap = "meth_or_NULL",
    sig = "character",
    tag = "character",
    analysis = "character"
    ),
  prototype = prototype(step = 0L,
    tag = character(1),
    analysis = character(1),
    meth = NULL,
    snap = NULL,
    params = list(wd = ".", remote = "remote"),
    others = list()
    ))

setMethod("initialize", "job",
  function(.Object, ...) {
    .Object <- callNextMethod()
    if (is.null(.Object@meth)) {
      .Object@meth <- .meth()
    }
    if (is.null(.Object@snap)) {
      .Object@snap <- .snap()
    }
    .Object
  })


setClassUnion("numeric_or_character", c("numeric", "character"))

.feature <- setClass("feature",
  contains = c("vector"),
  representation = representation(
    object = "ANY", sig = "character",
    snap = "character", nature = "character", type = "character"),
  prototype = NULL)

.feature_list <- setClass("feature_list", contains = c("feature"))
.feature_char <- setClass("feature_char", contains = c("feature"))

setValidity("feature_list",
  function(object){
    is(object@.Data, "list")
  })

setMethod("[[", signature = c(x = "feature"),
  function(x, i, ...) {
    x@.Data <- x@.Data[[ i ]]
    if (is(x, "feature_list") && is(x@.Data, "character")) {
      new <- .feature_char()
      for (i in slotNames(new)) {
        slot(new, i) <- slot(x, i)
      }
      x <- new
    }
    names(x) <- names(x@.Data)
    return(x)
  })

setMethod("[", signature = c(x = "feature"),
  function(x, i, ...){
    x@.Data <- x@.Data[ i ]
    names(x) <- names(x)[ i ]
    return(x)
  })


setValidity("feature_char",
  function(object){
    is(object@.Data, "character")
  })

setGeneric("feature",
  function(x, ...) standardGeneric("feature"))

setGeneric("feature<-",
  function(x, value) standardGeneric("feature<-"))

setMethod("feature", signature = c(x = "job"),
  function(x, ...){
    feas <- x$.feature
    if (!is(feas, "feature")) {
      feas <- as_feature(feas, x, ...)
    }
    feas
  })

setReplaceMethod("feature", signature = c(x = "job"),
  function(x, value){
    x$.feature <- value
    return(x)
  })

setReplaceMethod("feature", signature = c(x = "wrap"),
  function(x, value){
    x$.feature <- value
    x
  })

setMethod("feature", signature = c(x = "wrap"),
  function(x){
    x$.feature
  })

setGeneric("as_feature",
  function(x, ref, ...) standardGeneric("as_feature"))

setMethod("as_feature", signature = c(x = "ANY", ref = "job"),
  function(x, ref, nature = c("genes", "compounds", "feature"),
    type = "disease", analysis = NULL)
  {
    if (is.null(sig(ref))) {
      stop('is.null(sig(ref)), no "sig" in that "job".')
    }
    nature <- match.arg(nature)
    nature <- switch(nature, "genes" = "基因集", "compounds" = "化合物", "feature" = "特征集")
    type <- dplyr::recode(type, "disease" = "疾病", .default = type)
    if (is(x, "character")) {
      x <- .feature_char(x, type = type, nature = nature)
    } else if (is(x, "list")) {
      x <- .feature_list(x, type = type, nature = nature)
    } else {
      stop("Now, only 'character' or 'list' types of 'feature' support.")
      x <- .feature(x, type = type, nature = nature)
    }
    if (missing(analysis)) {
      analysis <- ref@analysis
    } else {
      ## upper the generic environment.
      analysis <- glue::glue(analysis, .envir = parent.frame(2))
    }
    sig <- sig(ref)
    if (length(sig)) {
      sig <- glue::glue("[Section: {sig}]")
    } else {
      sig <- ""
    }
    snap(x) <- glue::glue("来自于{analysis}{sig}")
    return(x)
  })

setGeneric("snap",
  function(x, ref, ...) standardGeneric("snap"))

setGeneric("snap<-",
  function(x, value, ...) standardGeneric("snap<-"))

setMethod("snap", signature = c(x = "ANY", ref = "missing"),
  function(x){
    attr(x, ".SNAP")
  })

setReplaceMethod("snap", signature = c(x = "ANY"),
  function(x, value, ...){
    attr(x, ".SNAP") <- glue::glue(value, ..., .envir = parent.frame(2))
    return(x)
  })

setMethod("snap", signature = c(x = "job", ref = "missing"),
  function(x){
    x@snap
  })

setReplaceMethod("snap", signature = c(x = "job"),
  function(x, value){
    x@snap <- value
    return(x)
  })

setMethod("snap", signature = c(x = "feature", ref = "missing"),
  function(x){
    snap(x, FALSE)
  })

setMethod("snap", signature = c(x = "feature", ref = "logical"),
  function(x, ref = FALSE, limit = 50, num = 3){
    if (ref) {
      if (is(x, "feature_list")) {
        stop('is(x, "feature_list"), only "feature_char" support.')
      }
      if (length(x) > limit) {
        stop('length(x) > limit, too many features to show.')
      }
      glue::glue("{bind(x)} ({x@snap}) ")
    } else {
      if (is(x, "feature_list")) {
        str <- names(x)
        if (is.null(str)) {
          str <- ""
        }
      } else if (is(x, "feature_char")) {
        str <- x@.Data
      }
      if (length(str) >= num) {
        n <- length(str)
        if (length(str) == num) {
          add <- FALSE
        } else {
          add <- TRUE
        }
        str <- bind(head(str, n = num))
        if (nchar(str) > 25) {
          str <- stringr::str_trunc(str, limit)
        }
        if (add) {
          str <- paste0(str, glue::glue(", ...[n = {n}]"))
        }
      } else {
        str <- bind(str)
      }
      sep <- if (nchar(str)) ", " else ""
      glue::glue("{x@nature} ({str}{sep}{x@snap}) ")
    }
  })

setReplaceMethod("snap", signature = c(x = "feature"),
  function(x, value){
    x@snap <- value
    return(x)
  })

setMethod("snap", signature = c(x = "job", ref = "numeric_or_character"),
  function(x, ref){
    fun <- function(ref) {
      if (!is.null(snap <- snap(x)[[ paste0("step", ref) ]])) {
        snap
      } else {
        ""
      }
    }
    if (length(ref) > 1) {
      res <- unlist(lapply(ref, fun))
      paste0(res[ res != "" ], collapse = "")
    } else {
      fun(ref)
    }
  })

trace_filter <- function(data, ..., quosures = NULL) {
  cli::cli_alert_info("`trace_filter`: auto-recognizing can not do very well.")
  if (is.null(quosures)) {
    quosures <- rlang::enquos(...)
  }
  before <- data
  data <- dplyr::filter(data, !!!(quosures))
  cols <- colnames(data)
  labels <- lapply(quosures, rlang::as_label)
  snaps <- vapply(labels,
    function(label) {
      alls <- stringr::str_extract_all(label, "[a-zA-Z0-9_.]+|`.*?`")[[1]]
      alls <- gs(alls, "`", "")
      col <- alls[ alls %in% cols ]
      if (length(col) > 1 || !length(col)) {
        stop(glue::glue("`trace_filter` can not found the target column: {col} in {label}."))
      }
      if (is.numeric(data[[ col ]])) {
        symbol <- strx(label, ">=|<=|==|>|<")
        if (!length(symbol)) {
          stop("`trace_filter`: numeric only support: ", ">=|<=|==|>|<")
        }
        thred <- strx(label, "[0-9.]+")
        thred <- thred[ !is.na(thred) ]
        if (length(thred)) {
          thred <- as.numeric(thred)
        }
        if (symbol == ">=" || symbol == ">") {
          cutoff <- min(data[[ col ]])
        } else if (symbol == "<=" || symbol == "<") {
          cutoff <- max(data[[ col ]])
        } else if (symbol == "==") {
          cutoff <- unique(data[[ col ]])
          if (length(cutoff) > 1) {
            stop("`trace_filter`: expected '==', but not unique value.")
          }
        } else {
          stop(glue::glue("`trace_filter`: can not recognize the symbol: {symbol}"))
        }
        if (is.numeric(thred)) {
          if (eval(parse(text = paste0("cutoff ", symbol, " thred")))) {
            cutoff <- thred
          }
        }
        des <- switch(symbol,
          ">=" = "大于或等于",
          "<=" = "小于或等于",
          "==" = "等于",
          ">" = "大于", "<" = "小于"
        )
        glue::glue("筛选 {col} {des} {cutoff}")
      } else if (is.character(data[[ col ]]) || is.factor(data[[ col ]])) {
        if (grpl(label, "gr[e]?pl")) {
          if (grpl(label, "!gr[e]?pl")) {
            neg <- "不"
          } else neg <- ""
          match <- strx(label, '".*"|\'.*\'')
          if (is.na(match)) {
            stop(
              'is.na(match), detected "grepl" or "grpl", but not matched the pattern.'
            )
          }
          glue::glue("匹配 {col} 中{neg}包含{match}的描述")
        } else {
          types <- unique(data[[ col ]])
          if (length(types) > 20) {
            stop("`trace_filter`: too many '{col}' ({length(types)}) to display. ")
          }
          glue::glue("筛选 {col} 为 {bind(types, quote = TRUE)}")
        }
      } else {
        stop(
          "`trace_filter`: detected: ",
          col, "(", class(data[[ col ]]), "), Only 'character' or 'numeric' support."
        )
      }
    }, character(1))
  if (nrow(data) == nrow(before)) {
    message("trace_filter: no filter out.")
    snap <- paste0(bind(snaps), glue::glue("，未过滤掉任何数据 (n = {nrow(data)})。"))
  } else {
    snap <- paste0(bind(snaps), glue::glue("，最终得到 {nrow(data)} 例数据。"))
  }
  message(snap)
  snap(data, .open = "{{{", .close = "}}}") <- snap
  return(data)
}

setGeneric("try_snap",
  function(x, ...) standardGeneric("try_snap"))

setMethod("try_snap", signature = c(x = "data.frame"),
  function(x, main, sub){
    items <- split(x[[ sub ]], x[[ main ]])
    len <- vapply(items, function(x) length(unique(x)), integer(1))
    bind(paste0(names(items), " (n=", len, ") "))
  })

setGeneric("try_summary",
  function(x, ...) standardGeneric("try_summary"))

setMethod("try_summary", signature = c(x = "list"),
  function(x, merge = TRUE){
    lens <- lengths(x)
    alls <- list(max = names(lens)[ lens == max(lens) ],
      min = names(lens)[ lens == min(lens) ],
      mean = round(mean(lens))
    )
    alls$max <- glue::glue("数目最多的为 {bind(alls$max)} (n={max(lens)})")
    alls$min <- glue::glue("数目最少的为 {bind(alls$min)} (n={min(lens)})")
    alls$mean <- glue::glue("平均数目 {alls$mean}")
    if (merge) {
      alls <- paste0(alls$max, "；", alls$min, "；", alls$mean, "。")
    }
    return(alls)
  })

setClassUnion("character_or_factor", c("character", "factor"))

setMethod("try_snap", signature = c(x = "character_or_factor"),
  function(x){
    if (is(x, "factor")) {
      x <- as.character(x)
    }
    items <- split(x, x)
    len <- vapply(items, length, integer(1))
    bind(paste0(names(items), " (n=", len, ") "))
  })

setMethod("try_snap", signature = c(x = "list"),
  function(x, col){
    len <- vapply(x,
      function(x) {
        if (is(x, "data.frame")) {
          length(unique(x[[ col ]]))
        } else {
          length(unique(x))
        }
      }, integer(1))
    bind(paste0(names(x), " (n=", len, ") "))
  })

setGeneric("init",
  function(x, ...) standardGeneric("init"))

setGeneric("init<-",
  function(x, value) standardGeneric("init<-"))

setMethod("init", signature = c(x = "meth"),
  function(x){
    x@init
  })

setMethod("init", signature = c(x = "NULL"),
  function(x){
    TRUE
  })

setReplaceMethod("init", signature = c(x = "meth"),
  function(x, value){
    x@init <- value
    return(x)
  })

setGeneric("meth",
  function(x, ref, ...) standardGeneric("meth"))
setMethod("meth", signature = c(x = "ANY", ref = "missing"),
  function(x){
    x@meth
  })

setMethod("meth", signature = c(x = "job", ref = "logical"),
  function(x, ref){
    x <- meth(x)
    if (ref) {
      lapply(x@.xData$.order, function(name) x@.xData[[ name ]])
    } else {
      vapply(x@.xData$.order, function(name) x@.xData[[ name ]], character(1))
    }
  })

setMethod("meth", signature = c(x = "job", ref = "numeric_or_character"),
  function(x, ref){
    fun <- function(ref) {
      if (!is.null(meth <- meth(x)[[ paste0("step", ref) ]])) {
        meth
      } else {
        ""
      }
    }
    if (length(ref) > 1) {
      res <- unlist(lapply(ref, fun))
      paste0(res[ res != "" ], collapse = "")
    } else {
      fun(ref)
    }
  })


setGeneric("meth<-",
  function(x, value) standardGeneric("meth<-"))
setReplaceMethod("meth", signature = c(x = "ANY"),
  function(x, value){
    x@meth <- value
    return(x)
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

get_meth <- function(x, setHeader = TRUE) {
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
    c("", header, "", meth(x, TRUE))
  } else {
    c("", meth(x, TRUE), "")
  }
}

setMethod("show", signature = c(object = "meth"),
  function(object){
    lapply(object@.xData$.order,
      function(name) {
        writeLines(
          paste(crayon::silver(paste0(name, ":")), object@.xData[[ name ]])
        )
      })
    writeLines("\n")
  })

# setValidity("job",
#   function(object){
#     class <- class(object)
#     options(method_name = class)
#     if (is(object, "job")) TRUE else FALSE
#   })

setClass("hclust")
## cutree, the tree higher, the sort number lower

.marker_list <- setClass("marker_list",
  contains = c(),
  representation = representation(level = "numeric", marker = "list"),
  prototype = NULL)

setMethod("show", signature = c(object = "marker_list"),
  function(object){
    .show(object)
  })

#' @exportMethod show
setMethod("show", signature = c(object = "job"),
  function(object){
    message("A workflow of '", class(object), "' in step (done): ", object@step)
    message("Object size: ", obj.size(object))
    if (!is.null(object@params$set_remote))
      message("Remote: ", object@params$remote)
  })

setGeneric("hunt",
  function(x, ref, ...) standardGeneric("hunt"))
setGeneric("upd",
  function(x, ...) standardGeneric("upd"))

.db_expect <- setClass(
  "db_expect", contains = c("list"), 
  representation = representation(db_file = "character")
)

setMethod("show", signature = c(object = "db_expect"),
  function(object){
    content <- c(
      paste(crayon::yellow("DB file:"), object@db_file), 
      paste(crayon::silver("Number:"), length(object)),
      paste(crayon::silver("Duplicated:", length(which(duplicated(object@.Data)))))
    )
    writeLines(content)
  })

new_db_expect <- function(db_file, create_dir = TRUE) {
  if (!length(db_file) || length(db_file) > 1) {
    stop('!length(db_file) || length(db_file) > 1, illegal file name.')
  }
  if (file.exists(db_file)) {
    message(glue::glue("The file '{db_file}' already exists."))
  } else {
    if (create_dir) {
      dir <- dirname(db_file)
      dir.create(dir, FALSE, recursive = TRUE)
    }
  }
  .db_expect(db_file = db_file)
}

.hash_expect <- setClass("hash_expect",
  contains = c("character"),
  representation = representation(hash = "character"),
  prototype = NULL)

setMethod("show", signature = c(object = "hash_expect"),
  function(object){
    writeLines(c(crayon::silver(object@hash), bind(object@.Data)))
  })

setGeneric("hash_expect",
  function(x, ...) standardGeneric("hash_expect"))
setGeneric("hash_expect<-",
  function(x, value) standardGeneric("hash_expect<-"))
setMethod("hash_expect", signature = c(x = "ANY"),
  function(x){
    attr(x, "__HASH_EXPECT__")
  })
setReplaceMethod("hash_expect", signature = c(x = "ANY"),
  function(x, value){
    if (!is(value, "hash_expect")) {
      stop('!is(value, "hash_expect"), value should be a "hash_expect" object.')
    }
    attr(x, "__HASH_EXPECT__") <- value
    return(x)
  })

setMethod("hunt", signature = c(x = "expect_cols", ref = "ANY"),
  function(x, ref){
    hunt(suppressMessages(new_db_expect(x@db_file)), ref)
  })

setMethod("hunt", signature = c(x = "db_expect", ref = "hash_expect"),
  function(x, ref){
    if (!length(x)) {
      if (length(x@db_file) && file.exists(x@db_file)) {
        x@.Data <- readRDS(x@db_file)
      }
    }
    x[[ ref@hash ]]
  })

setMethod("hunt", signature = c(x = "db_expect", ref = "ANY"),
  function(x, ref){
    hunt(x, hash_expect(ref))
  })

setMethod("hunt", signature = c(x = "db_expect", ref = "NULL"),
  function(x, ref){
    stop("`ref` should not be NULL, or without 'hash_expect'.")
  })

setMethod("upd", signature = c(x = "expect_cols"),
  function(x, data, ...){
    upd(suppressMessages(new_db_expect(x@db_file)), data, ...)
  })

setMethod("upd", signature = c(x = "db_expect"),
  function(x, data, save = TRUE){
    ref <- hash_expect(data)
    if (!is(ref, "hash_expect")) {
      stop('!is(hash_expect(data), "hash_expect"), should be "hash_expect" object.')
    }
    if (!length(x@.Data)) {
      if (!length(x@db_file)) {
        stop('!length(x@db_file), no "db_file"?')
      }
      if (file.exists(x@db_file)) {
        x@.Data <- readRDS(x@db_file)
      }
    }
    item <- sapply(
      ref@.Data, function(name) data[[ name ]], simplify = FALSE
    )
    if (!is.null(x[[ ref@hash ]])) {
      if (identical(x[[ ref@hash ]], item)) {
        message("Already exists, and identical, need not update.")
        return(x)
      } else {
        message("Already exists, update now.")
      }
    } else {
      message("Add new record.")
    }
    x[[ ref@hash ]] <- item
    if (save) {
      lst <- x@.Data
      names(lst) <- names(x)
      saveRDS(lst, x@db_file)
    }
    return(x)
  })

setGeneric("expect",
  function(x, ref, ...) standardGeneric("expect"))

set_hash_expect <- function(x, ref, id = NULL) {
  if (length(ref@uniqueness)) {
    if (is.null(x[[ref@uniqueness]])) {
      stop('is.null(x[[ref@uniqueness]]), element should not be NULL.')
    }
    if (is(x, "df")) {
      len <- nrow(x)
    } else {
      len <- length(x)
    }
    if (!len) {
      stop('!len, should has at least one depth.')
    }
    hash <- digest::digest(c(id, x[[ref@uniqueness]]), "md5")
    hash_expect(x) <- .hash_expect(
      vapply(ref, function(obj) obj@name, character(1)),
      hash = hash
    )  
  }
  return(x)
}

setMethod("expect", signature = c(x = "data.frame", ref = "expect_cols"),
  function(x, ref, force = FALSE, id = NULL, sleep = NULL,
    silent_select = TRUE, keep = TRUE)
  {
    if (!length(ref)) {
      stop('!length(ref), not found any "expect_col".')
    }
    if (!is.null(id)) {
      cli::cli_h2(paste("START:", id))
    }
    if (length(ref@uniqueness)) {
      x <- set_hash_expect(x, ref, id)
      if (length(ref@db_file)) {
        records <- hunt(ref, x)
        if (!is.null(records) && !force) {
          message(
            crayon::yellow("Recodes exists, use that.\n"),
            paste0(
              vapply(records, showStrings, character(1)), collapse = "\n"
            )
          )
          x <- dplyr::mutate(x, !!!records)
          if (!is.null(ref@global)) {
            x <- ref@global(x)
          }
          if (!is.null(id)) {
            cli::cli_h3("END")
          }
          return(x)
        }
      }
    }
    showChoices <- function(colnames) {
      choices <- vapply(colnames, 
        function(name) {
          paste0(name, ": ", bind("\'", head(x[[name]], n = 20), "\'"))
        }, character(1))
      ## Add unknown.
      choices <- c(choices, "Unknown")
      stringr::str_trunc(choices, 160L)
    }
    lapply(ref,
      function(i) {
        message(crayon::silver("+++ "), crayon::blue(i@name), crayon::silver(" +++"))
        which <- integer(0L)
        for (pat in i@pattern_find) {
          if (length(which <- grp(colnames(x), pat))) {
            if (!is.null(i@fun_check)) {
              isThats <- vapply(which,
                function(n) {
                  if (is.null(i@fun_mutate)) {
                    i@fun_check(x[[n]]) 
                  } else {
                    res <- try(i@fun_mutate(x[[n]]), TRUE)
                    if (inherits(res, "try-error")) {
                      return(FALSE)
                    } else {
                      i@fun_check(res)
                    }
                  }
                }, logical(1))
              if (all(!isThats)) {
                which <- integer(0)
                next
              } else {
                which <- which[ isThats ]
                break
              }
            } else break
          }
        }
        re_check <- FALSE
        if (!length(which)) {
          which <- menu(
            showChoices(colnames(x)), title = glue::glue(
              "Can not match '{i@name}' by pattern, please specify that."
            )
          )
          re_check <- TRUE
        } else if (length(which) > 1) {
          whiches <- which
          which <- menu(
            showChoices(colnames(x)[which]), title = glue::glue(
              "Match too many '{i@name}' with pattern, please specify one."
            )
          )
          which <- which(colnames(x) == colnames(x)[whiches][which])
          re_check <- TRUE
        }
        if (re_check && silent_select) {
          re_check <- FALSE
        }
        if (which <= length(x)) {
          if (re_check && !is.null(i@fun_check)) {
            if (!i@fun_check(x[[ which ]])) {
              if (!usethis::ui_yeah('!i@fun_check(x[[ which ]]), continue?')) {
                stop("Can not pass function checking.")
              }
            }
          }
          x[[ i@name ]] <- x[[ which ]]
          if (length(i@pattern_recode)) {
            if (is.null(names(i@pattern_recode))) {
              stop('is.null(names(i@pattern_recode)), should has names.')
            }
            meta <- group_strings(
              x[[ i@name ]], i@pattern_recode
            )
            matches <- match(x[[ i@name ]], meta[[ "target" ]])
            x[[ i@name ]] <- meta[[ "group" ]][ matches ]
            x[[ i@name ]] <- ifelse(
              is.na(x[[ i@name ]]), "Others", x[[ i@name ]]
            )
          }
          if (!is.null(i@fun_mutate)) {
            x[[ i@name ]] <- i@fun_mutate(x[[ i@name ]])
          }
          show <- stringr::str_wrap(
            bind(c(bind(paste0("\'", head(x[[i@name]], n = 10), "\'")),
                if (length(x[[i@name]]) > 10) "..." else NULL)),
            indent = 4, exdent = 4
          )
          note <- crayon::silver(paste0("(From column: ", colnames(x)[which], ")"))
          message(glue::glue("Matched {crayon::yellow(i@name)} {note}:\n{show}"))
          if (!keep && colnames(x)[which] != i@name) {
            x <- x[, -which]
          }
        } else {
          message(crayon::red("You selected 'Unknown', so the results will all be 'Unknown'."))
          x[[ i@name ]] <- "Unknown"
        }
        if (!is.null(sleep) && is.numeric(sleep)) {
          Sys.sleep(sleep)
        }
        x <<- x
      })
    if (!is.null(ref@global)) {
      x <- ref@global(x)
    }
    if (length(ref@uniqueness) && length(ref@db_file)) {
      upd(ref, x)
    }
    if (!is.null(id)) {
      cli::cli_h3("END")
    }
    return(x)
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

setGeneric("less",
  function(x, ...) standardGeneric("less"))

setGeneric("label",
  function(x, ...) standardGeneric("label"))
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
setGeneric("transmute",
  function(x, ref, ...) standardGeneric("transmute"))

setMethod("less", signature = c(x = "numeric_or_character"),
  function(x, n = 3, ..., quote = NULL){
    if (!is.null(snap(x))) {
      return(snap(x))
    }
    if (length(x) > n) {
      out <- head(x, n = n)
    } else {
      out <- x
    }
    if (is.null(quote) && any(grpl(out, "[^a-zA-Z0-9]"))) {
      out <- paste0("\"", out, "\"")
    }
    if (length(x) > n) {
      out <- c(out, glue::glue("...(n = {length(x)})"))
    }
    bind(out, ...)
  })

setMethod("transmute", signature = c(x = "feature", ref = "missing"),
  function(x){
    generator <- .test_job_generator(sys.function(2))
    transmute(x, generator())
  })

.test_job_generator <- function(method) {
  if (!isS4(method)) {
    stop('!isS4(method), was not in S4 generic call.')
  }
  convertTo <- sub("^[^_]*(_.*)", ".job\\1", method@generic)
  if (!is(generator <- get_fun(convertTo), "classGeneratorFunction")) {
    stop(
      '!is(generator <- get_fun(convertTo), "classGeneratorFunction"), the matched `',
      convertTo, '` is not a "classGeneratorFunction".'
    )
  }
  return(generator)
}

setMethod("transmute", signature = c(x = "feature", ref = "job"),
  function(x, ref){
    glue::glue("对{snap(x)}进行{ref@analysis}。")
  })

setMethod("transmute", signature = c(x = "feature", ref = "character"),
  function(x, ref){
    glue::glue("对{snap(x)}进行{ref}。")
  })

.character_ref <- setClass("character_ref",
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

setMethod("ref", signature = c(x = "character_ref"),
  function(x, add_internal_job = TRUE, ...){
    ref <- strx(x, "[A-Z][A-Za-z]{10,}[0-9]{4}")
    message("Extracted reference: ", ref)
    if (add_internal_job) {
      job <- .job_publish(cite = paste0("[@", ref, "]"))
      .add_internal_job(job)
    }
    return(x)
  })

setMethod("label", signature = c(x = "ANY"),
  function(x, ...){
    if (is.null(lab(x))) {
      stop('is.null(lab(x)), no "lab" found, use `lab` to set it.')
    }
    x <- as_chunk_label(lab(x))
    ref(x, ...)
  })

setMethod("ref", signature = c(x = "character"),
  function(x, legend = TRUE, try = FALSE, ...){
    codes_list <- knitr::knit_code$get()
    if (!length(codes_list)) {
      message("Not in rendering circumstance.")
      return("")
    }
    codes <- codes_list[[x]]
    if (is.null(codes)) {
      if (try) {
        return("")
      }
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
    if (legend) {
      placeHolder <- glue::glue("<!-- autor_legend:{x} -->")
    } else {
      placeHolder <- ""
    }
    if (is(obj, "df")) {
      glue::glue("Tab. \\@ref(tab:{x}) {placeHolder}")
    } else if (is(obj, "can_not_be_draw") || is(obj, "can_be_draw") || is(obj, "fig")) {
      glue::glue("Fig. \\@ref(fig:{x}) {placeHolder}")
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
  Hmisc::capitalize(gs(gs(x, "(?<![a-zA-Z])ids[ \\-]?", "", perl = TRUE), " |_|\\.", "-"))
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
    x@sig <- value
    return(x)
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
  function(x, is.remote = FALSE,
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

touch_dir <- function(dir_output, check = TRUE) {
  if (check && dir.exists(dir_output)) {
    if (sureThat("dir.exists(dir_output), removed?")) {
      unlink(dir_output, TRUE)
    }
  }
  dir.create(dir_output, FALSE)
  return(dir_output)
}

pg_local_recode <- function() {
  conda <- getOption("conda", "~/miniconda3")
  lst <- list(
    fusion = "~/fusion_twas",
    ldscPython = "{conda}/bin/conda run -n ldsc python",
    ldsc = "~/ldsc",
    annovar = "~/disk_sda1/annovar",
    vep = "~/ensembl-vep/vep",
    vep_cache = "~/disk_sda1/.vep",
    vina = "vina",
    python = "{conda}/bin/python3",
    conda = "{conda}/bin/conda",
    conda_env = "{conda}/envs",
    qiime = "{conda}/bin/conda run -n qiime2 qiime",
    musitePython = "{conda}/bin/conda run -n musite python3",
    musitePTM = "~/MusiteDeep_web/MusiteDeep/predict_multi_batch.py",
    musitePTM2S = "~/MusiteDeep_web/PTM2S/ptm2Structure.py",
    hobPython = "{conda}/bin/conda run -n hobpre python",
    hobPredict = "~/HOB/HOB_predict.py",
    hobModel = "~/HOB/model",
    hobExtra = "~/HOB/pca_hob.m",
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
    pymol = "/usr/bin/python3 -m pymol",
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
    x@object <- value
    return(x)
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
    x@object <- value
    return(x)
  })

setGeneric("Legend",
  function(x, ...) standardGeneric("Legend"))

setGeneric("Legend<-",
  function(x, value) standardGeneric("Legend<-"))

setMethod("Legend", signature = c(x = "ANY"),
  function(x){
    attr(x, ".LEGEND")
  })

setReplaceMethod("Legend", signature = c(x = "ANY"),
  function(x, value){
    attr(x, ".LEGEND") <- value
    return(x)
  })

setLegend <- function(x, from, ..., check_list = TRUE, env = parent.frame(1))
{
  if (check_list && is(x, "list") && length(from) == length(x)) {
    x <- mapply(x, from, FUN = 
      function(obj, leg) {
        Legend(obj) <- glue::glue(leg, .envir = env, ...)
        return(obj)
      }, SIMPLIFY = FALSE)
    return(x)
  } else {
    if (is.null(env)) {
      Legend(x) <- from
    } else {
      Legend(x) <- glue::glue(from, .envir = env, ...)
    }
    return(x)
  }
}

plotsAdd <- function(x, ...) {
  jobSlotAdd(x, "plots", ...)
}

tablesAdd <- function(x, ...) {
  jobSlotAdd(x, "tables", ...)
}

jobSlotAdd <- function(x, name, ..., reset, step) {
  calls <- getNamesOfCall(substitute(list(...)))
  objs <- list(...)
  names(objs) <- calls
  labs <- formatNames(calls)
  if (missing(reset)) {
    symbol <- paste0(".", name, "_initial")
    isInitial <- x[[ symbol ]]
    if (is.null(isInitial) || isInitial) {
      x[[ symbol ]] <- FALSE
      reset <- TRUE
    } else {
      reset <- FALSE
    }  
  }
  if (missing(step)) {
    step <- x@step
  }
  db <- if (reset) {
    list()
  } else if (length(slot(x, name)) >= step) {
    slot(x, name)[[ step ]]
  } else {
    list()
  }
  db <- db[ !names(db) %in% names(objs) ]
  slot(x, name)[[ step ]] <- c(db, .set_lab_batch(objs, labs, sig(x)))
  return(x)
}

setGeneric("methodAdd",
  function(x, from, ...) standardGeneric("methodAdd"))

setMethod("methodAdd", signature = c(x = "job", from = "character"),
  function(x, from, add = FALSE, env = parent.frame(2), step = NULL) {
    if (missing(step)) {
      if (!is.null(x$.map_step) && x$.map_step) {
        step <- "m"
      } else {
        step <- x@step
      }
    }
    if (missing(add)) {
      if (init(meth(x))) {
        add <- FALSE
        init(meth(x)) <- FALSE
      } else {
        add <- TRUE
      }
    } else {
      if (!add) {
        init(meth(x)) <- FALSE
      }
    }
    if (add) {
      former <- meth(x)[[ paste0("step", x@step) ]]  
    } else {
      former <- ""
    }
    if (is.null(env)) {
      meth(x)[[ paste0("step", x@step) ]] <- paste0(former, from)
    } else {
      meth(x)[[ paste0("step", x@step) ]] <- paste0(
        former, glue::glue(from, .envir = env)
      )
    }
    return(x)
  })

setMethod("methodAdd", signature = c(x = "job", from = "job"),
  function(x, from, ...){
    methodAdd(x, meth(from), ...)
  })

setMethod("methodAdd", signature = c(x = "job", from = "meth"),
  function(x, from, ...){
    for (name in from@.xData$.order) {
      x <- methodAdd(x, from[[ name ]], ...)
    }
    return(x)
  })

setGeneric("snapAdd",
  function(x, from, ...) standardGeneric("snapAdd"))

setMethod("snapAdd", signature = c(x = "job", from = "character"),
  function(x, from, add = FALSE, env = parent.frame(2), step = NULL) {
    if (missing(step)) {
      if (!is.null(x$.map_step) && x$.map_step) {
        step <- "m"
      } else {
        step <- x@step
      }
    }
    if (missing(add)) {
      if (init(snap(x))) {
        add <- FALSE
        init(snap(x)) <- FALSE
      } else {
        add <- TRUE
      }
    } else {
      if (!add) {
        init(snap(x)) <- FALSE
      }
    }
    if (add) {
      former <- snap(x)[[ paste0("step", step) ]]
    } else {
      former <- ""
    }
    if (is.null(env)) {
      snap(x)[[ paste0("step", step) ]] <- paste0(former, from)
    } else {
      snap(x)[[ paste0("step", step) ]] <- paste0(former, glue::glue(from, .envir = env))
    }
    return(x)
  })

setMethod("snapAdd", signature = c(x = "job", from = "job"),
  function(x, from, ...){
    snapAdd(x, snap(from), ...)
  })

setMethod("snapAdd", signature = c(x = "job", from = "meth"),
  function(x, from, ...){
    for (name in from@.xData$.order) {
      x <- snapAdd(x, from[[ name ]], ...)
    }
    return(x)
  })

# static_call_on_exit <- function(static, dynamic, funName, args = "",
#   expr = "{dynamic} <- {funName}({dynamic}, {static}, {args})",
#   add = TRUE, after = TRUE,
#   env = parent.frame(1), saveStatic = ".__TEMP__")
# {
#   expr <- glue::glue("on.exit({expr}, add = {add}, after = {after})")
#   eval(parse(text = glue::glue(expr)), envir = env)
# }

resolve_feature_snapAdd_onExit <- function(x, feature, funName = "snapAdd",
  saveStatic = ".__SNAP__", env = parent.frame(1),
  fun_generator = sys.function(2))
{
  if (is(feature, "feature")) {
    generator <- .test_job_generator(fun_generator)
    snap <- transmute(feature, generator())
    snapAdd_onExit(x, snap, funName, saveStatic, env)
    return(feature@.Data)  
  } else {
    feature
  }
}

resolve_feature <- function(feature, unlist = TRUE, 
  use.names = FALSE, recursive = TRUE)
{
  if (is(feature, "feature")) {
    snap <- snap(feature)
    feature <- feature@.Data
    if (unlist) {
      feature <- unlist(feature, use.names = use.names, recursive = recursive)
    }
    snap(feature) <- snap
  } 
  return(feature)
}

snapAdd_onExit <- function(x, static, funName = "snapAdd",
  saveStatic = ".__SNAP__", env = parent.frame(1))
{
  methodAdd_onExit(x, static, funName, saveStatic, env)
}

methodAdd_onExit <- function(x, static, funName = "methodAdd",
  saveStatic = ".__METH__", env = parent.frame(1))
{
  if (!is(x, "character") || !is(static, "character")) {
    stop('!is(x, "character") || !is(static, "character"), should all be character.')
  }
  saveEnv <- try(get(saveStatic, envir = env), TRUE)
  if (inherits(saveEnv, "try-error")) {
    saveEnv <- .meth()
    assign(saveStatic, saveEnv, envir = env)
    hasSetUp <- FALSE
  } else {
    hasSetUp <- TRUE
  }
  name <- paste0("value", length(saveEnv) + 1L)
  saveEnv[[ name ]] <- glue::glue(static, .envir = env)
  if (!hasSetUp) {
    expr <- glue::glue(
      "{funName}(get({x}, env = env), saveEnv, env = NULL)"
    )
    withr::defer(eval(parse(text = expr)), envir = env, "last")
  }
}

.set_lab_batch <- function(objs, labs, sig) {
  mapply(objs, labs, SIMPLIFY = FALSE,
    FUN = function(obj, lab) {
      if (is.null(lab(obj))) {
        if (length(obj) >= 1 && is(obj, "list")) {
          if (is.null(names(obj))) {
            names <- paste0("n", seq_along(obj))
          } else {
            names <- names(obj)
          }
          lab <- paste0(lab, "_", names)
        }
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
    x@object <- value
    return(x)
  })


#' @exportMethod params
setMethod("params", signature = c(x = "job"),
  function(x){
    x@params
  })
setReplaceMethod("params", signature = c(x = "job"),
  function(x, value){
    x@object <- value
    return(x)
  })


#' @exportMethod others
setMethod("others", signature = c(x = "job"),
  function(x){
    x@others[[ x@step ]]
  })
setReplaceMethod("others", signature = c(x = "job"),
  function(x, value){
    x@object <- value
    return(x)
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
    x$set_remote <- TRUE
  if (is.null(x$remote))
    x$remote <- remote
  if (is.null(x$map_local))
    x$map_local <- map_local
  if (is.null(x$tmpdir))
    x$tmpdir <- tmpdir
  if (is.null(x$wait)) {
    x$wait <- 10L
  }
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

setGroupGeneric("asjob_series",
  function(x, ...) {
    stop("...")
  }
)

setGroupGeneric("step_series",
  function(x, ...) {
    stop("...")
  })

setGeneric("step0", group = list("step_series"),
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
    mess <- paste0(paste0(seq_along(signatures), ": ", signatures), collapse = "\n")
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

setGeneric("step1", group = list("step_series"),
  function(x, ...) {
    # if (identical(sig(x), character(0))) {
    legal <- TRUE
    if (identical(parent.frame(1), .GlobalEnv)) {
      name <- rlang::as_label(substitute(x))
      if (grepl("^[A-Za-z][A-Za-z0-9_.]*$", name)) {
        sig <- toupper(gs(gs(gs(name, "^[a-zA-Z0-9]+", ""), "\\.", " "), "^[ ]*", ""))
        attr(sig, "name") <- name
        sig(x) <- sig 
      } else {
        legal <- FALSE
      }
    } else {
      legal <- FALSE
    }
    # }
    if (is.null(seed <- getOption("step_seed"))) {
      x$seed <- 987456L
    } else {
      x$seed <- seed
    }
    x <- checkAddStep(x, 1L)
    x$.append_heading <- TRUE
    x <- standardGeneric("step1")
    x <- stepPostModify(x, 1)
    if (legal && x$.append_heading) {
      job_append_heading(x)
    }
    x$.append_heading <- NULL
    x
  })

job_append_heading <- function (x, mutate = TRUE, heading = NULL) {
  if (getOption("job_appending", FALSE)) {
    if (is(x, "job")) {
      if (missing(heading)) {
        heading <- x@analysis
      }
      if (length(x@sig)) {
        heading <- paste0(heading, " (", x@sig, ")")
      }
      .append_heading(heading, mutate)
    }
  }
}

getParentClasses <- function(class) {
  names(attributes(getClassDef(class))$contains)
}

findMaxStepMethod <- function(class, max = 12L, inherited = TRUE) {
  allMax <- vapply(
    c(getParentClasses(class), class),
    function(name) {
      for (i in seq_len(max)) {
        sigs <- findMethodSignatures(
          f = paste0("step", i), inherited = TRUE
        )
        if (!any(sigs == name)) {
          return(i - 1L)
        }
      }
    }, integer(1)
  )
  return(max(allMax))
}

job_append_method <- function(x) {
  if (getOption("job_appending", FALSE)) {
    if (is(x, "job")) {
      if (requireNamespace("nvimcom")) {
        oname <- substitute(x)
        max <- findMaxStepMethod(class(x))
        args <- glue::glue("CheckAndAppendMethod('{oname}', {max})")
        .C("nvimcom_msg_to_nvim", args, PACKAGE = "nvimcom")
      }
    }
  }
}

.append_heading <- function(x, mutate = FALSE) {
  args <- paste0("CheckAndAppendHeading('", x, "'", if (mutate) ", 1" else ", 0", ")")
  if (requireNamespace("nvimcom")) {
    if (nchar(x)) {
      .C("nvimcom_msg_to_nvim", args, PACKAGE = "nvimcom")
    }
  }
}

.append_class <- function(expr, short = TRUE, envir = .GlobalEnv) {
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

setGeneric("step2", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 2L)
    x <- standardGeneric("step2")
    stepPostModify(x, 2)
  })

setGeneric("step3", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 3L)
    x <- standardGeneric("step3")
    stepPostModify(x, 3)
  })

setGeneric("step4", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 4L)
    x <- standardGeneric("step4")
    stepPostModify(x, 4)
  })

setGeneric("step5", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 5L)
    x <- standardGeneric("step5")
    stepPostModify(x, 5)
  })

setGeneric("step6", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 6L)
    x <- standardGeneric("step6")
    stepPostModify(x, 6)
  })

setGeneric("step7", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 7L)
    x <- standardGeneric("step7")
    stepPostModify(x, 7)
  })

setGeneric("step8", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 8L)
    x <- standardGeneric("step8")
    stepPostModify(x, 8)
  })

setGeneric("step9", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 9L)
    x <- standardGeneric("step9")
    stepPostModify(x, 9)
  })

setGeneric("step10", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 10L)
    x <- standardGeneric("step10")
    stepPostModify(x, 10)
  })

setGeneric("step11", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 11L)
    x <- standardGeneric("step11")
    stepPostModify(x, 11)
  })

setGeneric("step12", group = list("step_series"),
  function(x, ...) {
    x <- checkAddStep(x, 12L)
    x <- standardGeneric("step12")
    stepPostModify(x, 12)
  })

stepPostModify <- function(x, n = NULL) {
  cli::cli_h1("Job finished & Start post modify")
  xname <- attr(x@sig, "name")
  news <- c()
  if (!is.null(n)) {
    if (getOption("auto_convert_plots", FALSE)) {
      x <- convertPlots(x, n)
    }
    if (length(x@plots) >= n) {
      if (!is.null(x@plots[[ n ]])) {
        names(x@plots)[ n ] <- paste0("step", n)
      }
      if (!is.null(xname)) {
        new_plots <- paste0(xname, "@plots$step", x@step, "$", names(x@plots[[ x@step ]]))
        cli::cli_alert_info("Plot news:")
        lapply(new_plots, cli::cli_code)
        news <- c(news, new_plots)
      }
    }
    if (length(x@tables) >= n) {
      if (!is.null(x@tables[[ n ]]))
        names(x@tables)[ n ] <- paste0("step", n)
      if (!is.null(xname)) {
        new_tables <- paste0(xname, "@tables$step", x@step, "$", names(x@tables[[ x@step ]]))
        cli::cli_alert_info("Table news:")
        lapply(new_tables, cli::cli_code)
        news <- c(news, new_tables)
      }
    }
  }
  if (!is.null(xname)) {
    isNewParams <- !names(x@params) %in% x@others$.oldParams
    if (any(isNewParams)) {
      cli::cli_alert_info("Params news:")
      new_params <- paste0(xname, "@params$", names(x@params)[ isNewParams ])
      lapply(new_params, cli::cli_code)
      news <- c(news, new_params)
    }
    if (!is.null(xname) && length(xname)) {
      assign(xname, x)
      writeJobSlotsAutoCompletion(xname, environment())
    }
    meth <- meth(x, x@step)
    if (any(nchar(meth))) {
      if (length(meth) > 1) {
        meth <- paste0(meth, collapse = "")
      }
      message("METH:\n", .strwrapCh(meth))
    }
    snap <- snap(x, x@step)
    if (any(nchar(snap))) {
      if (length(snap) > 1) {
        snap <- paste0(snap, collapse = "")
      }
      message("SNAP:\n", .strwrapCh(snap))
    }
  }
  x
}

.strwrapCh <- function(x) {
  paste0(strwrap(split_text_by_width(x, 70), , 4), collapse = "\n")
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
        if (inherits(try(validObject(obj), TRUE), "try-error")) {
          message(glue::glue("Object invalid, update that: {x}"))
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
  dir.create(dir, FALSE, recursive = TRUE)
  obj <- try(get(x, envir = envir), TRUE)
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
  gc <- FALSE
  if (length(x@plots) >= n) {
    if (object.size(x@plots[[ n ]]) > 1e8) {
      gc <- TRUE
    }
    x@plots[[ n ]] <- rapply2_asWrap(x@plots[[ n ]], "gg.obj")
  }
  if (gc) {
    gc()
  }
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
  } else if (identical(class(x), "list") && length(x)) {
    attrs <- attributes(x)
    names <- names(x)
    x <- lapply(seq_along(x),
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

checkAddStep <- function(x, n, clear_tables_plots = TRUE) {
  STEP <- paste0("step", n)
  if (getOption("step_check", TRUE)) {
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
  init(snap(x)) <- TRUE
  init(meth(x)) <- TRUE
  x$.tables_initial <- TRUE
  x$.plots_initial <- TRUE
  x
}

step_message <- function(..., show_end = "Job Status") {
  str <- paste0(unlist(list(...)), collapse = "")
  str <- gs(str, "red\\{\\{(.*?)\\}\\}", "\033[31m\\1\033[39m")
  str <- gs(str, "grey\\{\\{(.*?)\\}\\}", "\033[90m\\1\033[39m")
  str <- gs(str, "yellow\\{\\{(.*?)\\}\\}", "\033[33m\\1\033[39m")
  str <- gs(str, "blue\\{\\{(.*?)\\}\\}", "\033[34m\\1\033[39m")
  textSh(str, pre_trunc = FALSE)
  if (!is.null(show_end))
    cli::cli_h1(show_end)
}

.argEnv <- new.env()

E <- function(expr, text = NULL, internal = FALSE) {
  e(expr, text, internal, n = 1)
}

e <- function(expr, text = NULL, internal = TRUE, n = NULL) {
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

setMethod("show", signature = c(object = "standardGeneric"),
  function(object){
    .show_method(object, show_standardGeneric, FALSE)
  })

setMethod("show", signature = c(object = "nonstandardGenericFunction"),
  function(object) {
    .show_method(object, show_nonstandardGenericFunction, TRUE)
  }
)

setMethod("print", signature = c(x = "nonstandardGenericFunction"),
  function(x){
    show(x)
  })

setGeneric("clear",
  function(x, ...) standardGeneric("clear"))

setMethod("clear", signature = c(x = "job"),
  function(x, save = TRUE, suffix = NULL, name = substitute(x, parent.frame(1))){
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
  function(x, ref, ...) {
    if (is(x, "job")) {
      init(snap(x)) <- TRUE
      init(meth(x)) <- TRUE
      x$.map_step <- TRUE
    }
    x <- standardGeneric("map")
    x$.map_step <- FALSE
    if (is(x, "job")) {
      if (identical(parent.frame(1), .GlobalEnv)) {
        if (!is.null(x$.map_heading)) {
          job_append_heading(x, heading = x$.map_heading)
        }
      }
      x$.map_heading <- NULL
      stepPostModify(x)
    } else {
      x
    }
  })

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

setGeneric("upload",
  function(x, ...) standardGeneric("upload"))

setGeneric("pull",
  function(x, ...) standardGeneric("pull"))

setMethod("intersect", signature = c(x = "ANY", y = "ANY"),
  function(x, y){
    base::intersect(x, y)
  })

setGeneric("not",
  function(x, ...) standardGeneric("not"))

setMethod("not", signature = c(x = "job"),
  function(x){
    options(method_name = class(x))
    invisible(x)
  })

setGeneric("merge",
  function(x, y, ...) standardGeneric("merge"))

setMethod("upd", signature = c(x = "job"),
  function(x){
    new <- new(class(x))
    new@object <- x@object
    new@params <- x@params
    if (inherits(try(x@snap@.xData, TRUE), "try-error")) {
      new@snap <- .snap()
      if (!is.null(x@params$.snap) && is(x@params$.snap, "list")) {
        new@snap@.xData <- as.environment(x@params$.snap)
        new@snap@.xData$.order <- names(x@params$.snap)
        new@params$.snap <- NULL  
      }
    }
    new@plots <- x@plots
    new@tables <- x@tables
    new@step <- x@step
    if (inherits(try(x@meth@.xData, TRUE), "try-error")) {
      new@meth <- .meth()
      lst <- as.list(unlist(x@meth))
      new@meth@.xData <- as.environment(lst)
      new@meth@.xData$.order <- names(lst)
    }
    if (!inherits(try(x@sig, TRUE), "try-error")) {
      new@sig <- x@sig
    }
    new
  })

setMethod("fill", signature = c(x = "df"),
  function(x, ref, fill, lst){
    for (i in seq_along(lst)) {
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
      if (x$set_remote) TRUE else FALSE
    } else FALSE
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
  function(x, ref, y, y.ref, y.get, rename = TRUE, col = NULL)
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
    if (any(grpl(x, " |/"))) {
      message(
        crayon::red("Deal with gene names, multiple gene names covered in One name, use first.")
      )
      x <- strx(x, "[^ ^/]+")
    }
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

treeObj <- function(obj, stops = c("gg", "wrap", "data.frame"), maxdepth = 10L, envir = parent.frame(1))
{
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
            paste0(oname, "[[", seq_along(obj), "]]"),
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

view_obj_for_vim <- function(x, y, view = TRUE) {
  if (isGeneric(x)) {
    if (nchar(y)) {
      classY <- class(eval(parse(text = y), envir = .GlobalEnv))
      if (length(classY) > 1) {
        classY <- classY[1]
      }
      message("Class of target: ", classY)
      obj <- try(selectMethod(x, classY), silent = TRUE)
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
          alls <- alls[alls[, 1] %in% getParentClasses(classY), ]
        }
        if (is.matrix(alls)) {
          useWhich <- menuThat(apply(alls, 1, paste0, collapse = ", "), "Show which?")
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
    if (as.double(strx(obj.size(eval(parse(text = x))), "[0-9.]+")) > .5) {
      stop(
        'as.double(strx(obj.size(eval(parse(text = x))), "[0-9.]+")) > .5, too large object to show.'
      )
    }
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
  x <- fun(plot_workflow_summary(later, TRUE))
  y <- fun(plot_workflow_summary(former, TRUE))
  x[ !x %in% y ]
}

plot_workflow_summary <- function(data, get_jobs_used = FALSE) {
  fun.co_job <- function() {
    alls <- available_signatures("step1")
    name <- gs(alls, "^job_", "")
    nodes <- data.frame(name = name)
    asjobs <- paste0("asjob_", name)
    asjobs <- sapply(asjobs, simplify = FALSE,
      function(x) {
        res <- try(get_fun(x), TRUE)
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
    graph <- fast_layout(edges, "linear", nodes = nodes, circular = TRUE)
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
    names(used) <- paste0(alls$receive_date, "_", seq_along(used), "##", alls$.type)
    whichLen0 <- lengths(used) == 0
    message("No any 'job' detected: \n\t", paste0(alls$.dir[whichLen0], collapse = "\n\t "))
    used <- lst_clear0(used)
    used <- as_df.lst(used)
    used <- dplyr::mutate(used, value = 1,
      order = paste0("Order ", gs(type, "##.*$", "")),
      rank = as.integer(as.Date(type)),
      name = gs(name, "^job_|^new_|(?<=^pb).*", "", perl = TRUE),
      name = gs(name, "^pb$", "Publication"),
      type = gs(type, ".*##(.*)$", "\\1")
    )
    data <- as_tibble(used)
    p.ind.pop <- ggplot(data, aes(x = reorder(order, rank, decreasing = TRUE), y = value)) +
      geom_col(aes(fill = name), position = "stack", width = .6) +
      geom_point(
        data = dplyr::distinct(data, order, type, rank),
        aes(reorder(order, rank, decreasing = TRUE), y = -3, shape = type, color = type),
        size = 5) +
      coord_flip() +
      labs(y = "Order (Date + N)", x = "Frequence", fill = "Name",
        color = "Type", shape = "Type") +
      scale_fill_manual(values = color_set(TRUE)) +
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
    pal <- MCnebula2:::.as_dic(color_set(), type, na.rm = TRUE)
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
    vps <- lapply(seq_along(seq),
      function(n, fx = .5, fy = .5, x = .5) {
        ang <- seq[n]
        ang2 <- seq2[n]
        vp <- viewport(cospi(ang / 180) * fx + x, sinpi(ang / 180) * fy + x,
          u(2, line), u(2, line), angle = ang2,
          just = c("centre", "bottom"))
        ggather(rectGrob(gp = gpar(fill = "lightyellow", col = pal[["sub.slot"]])),
          gtext(paste0("n", rev(seq_along(seq))[n]),
            gpar(fontface = "plain"), form = FALSE), vp = vp)
      })
    vps <- do.call(ggather, c(vps, list(vp = viewport(, .1, .5, .5))))
    grobs.chi$params %<>% into(vps)
    ## signal ggset
    ggsets <- frame_col(c(n = 1, x = 1, mn = 1),
      list(n = gtext("n", list(cex = 2.5), form = FALSE),
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
  function(x, ids, force = FALSE){
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
      stop('is(data, "df") == FALSE')
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

saves <- function(file = "workflow.rdata", saveInjobs = TRUE, ...) {
  injobs <- getOption("internal_job")
  if (!is.null(injobs) && saveInjobs) {
    rds <- .get_savesInternalRds(file)
    saveRDS(injobs, rds)
  }
  if (file.exists(file)) {
    isThat <- sureThat("The `file` exists, overwrite that?")
    if (!isThat) {
      return(NULL)
    }
  }
  save.image(file, ...)
}

menuThat <- function(qs, x, ...) {
  if (getOption("job_appending", FALSE) && requireNamespace("nvimcom")) {
    choices <- bind(
      paste0("\"", c(x, "exit", qs), "\""), co = ", "
    )
    .C("nvimcom_msg_to_nvim", glue::glue('RSelect({choices})'), PACKAGE = "nvimcom")
  }
  menu(qs, x, ...)
}

sureThat <- function(x, yes = c("Yes", "Sure", "Yup",
    "Yeah", "Agree", "Absolutely"), no = c("No way", "Not now", 
    "Negative", "No", "Nope", "Never"), n_yes = 1, n_no = 2, 
  shuffle = TRUE, .envir = parent.frame())
{
  x <- glue::glue_collapse(x, "\n")
  x <- glue::glue(x, .envir = .envir)
  if (!rlang::is_interactive()) {
    usethis::ui_stop(c("User input required, but session is not interactive.",
        "Query: {x}"))
  }
  n_yes <- min(n_yes, length(yes))
  n_no <- min(n_no, length(no))
  if (shuffle) {
    qs <- c(sample(yes, n_yes), sample(no, n_no))
    qs <- sample(qs)
  } else {
    qs <- c(yes[n_yes], no[n_no])
  }
  if (getOption("job_appending", FALSE) && requireNamespace("nvimcom")) {
    choices <- bind(
      paste0("\"", c(x, "exit", qs), "\""), co = ", "
    )
    .C("nvimcom_msg_to_nvim", glue::glue('RSelect({choices})'), PACKAGE = "nvimcom")
  }
  rlang::inform(x)
  out <- utils::menu(qs)
  out != 0L && qs[[out]] %in% yes
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

.auto_loads <- function() {
  if (!length(ls(envir = .GlobalEnv))) {
    if (length(save <- list.files(".", "workflow.*rdata"))) {
      if (sureThat("Workflow file exists, load that?", shuffle = FALSE)) {
        if (length(save) > 1) {
          save <- save[menuThat(save, "Which file you want to load?")]
        }
        loads(save)
      }
    }
  }
}

remotejob <- function(wd, remote = "remote") {
  .job(params = list(set_remote = TRUE, wd = wd, remote = remote))
}

set_prefix <- function(wd, db = wd, op = wd) {
  options(wd_prefix = wd, db_prefix = db, op_prefix = op)
}

.prefix <- function(path, name = c("wd", "db", "op"), warning = TRUE)
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
