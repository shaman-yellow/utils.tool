# ==========================================================================
# workflow of anno
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.assist_anno <- setClass("assist_anno",
  contains = c("virtual_job", "list"),
  representation = representation(
    object = "ANY",
    others = "ANY",
    info = "character",
    sig = "character",
    plots = "list",
    tables = "list",
    step = "integer"),
  prototype = prototype(
    step = 0L,
    pg = "anno",
    info = c("Tutorial: ...")
    ))

assist_anno <- function(file_md = "index.Rmd") {
  x <- .assist_anno()
  x$data <- get_md_titles(file_md)
  return(x)
}

setMethod("step1", signature = c(x = "assist_anno"),
  function(x, get = NULL, target = "# 附：分析流程 {#workflow}", result = "# 分析结果 {#results}")
  {
    step_message("Obtain document title.", show_end = "Assistant annotation")
    x$target <- target
    x$dataTarget <- matchMdContent(x$data, target, T)
    x$dataResult <- matchMdContent(x$data, result, T)
    x$query <- c(gptMessages()$heading, "",
      "实际呈现的标题：", "", x$dataResult$lines[-1], "",
      "提示性标题：", "", x$dataTarget$lines)
    if (is.null(get)) {
      get <- !getOption("assist_anno_hasGet", F)
    }
    if (get) {
      gett(x$query)
      message(crayon::yellow("Query has sent to clipboard."))
      options(assist_anno_hasGet = T)
    }
    message(crayon::red("Please manually copy results from ChatGpt."))
    return(x)
  })

setMethod("step2", signature = c(x = "assist_anno"),
  function(x, anno = get_clipboard(), target = T, result = T)
  {
    step_message("Embedded AI annotation text.")
    if (!requireNamespace("nvimcom")) {
      stop("nvimcom required.")
    }
    if (!all(x$dataTarget$lines %in% anno)) {
      stop(crayon::red("No paste into clipboard?"))
    }
    x$anno <- get_md_titles(lines = anno)
    fun_append <- function(meta, anno) {
      lapply(nrow(meta):1,
        function(n) {
          target <- meta$lines[n]
          anno <- matchMdContent(anno, target, mode = "d")
          if (is.null(anno)) {
            message("\nNo results for matching:\n", target)
            return()
          }
          anno <- anno[-1, ]
          if (any(nchar(anno$lines) >= 1)) {
            line <- meta$number[n]
            content <- anno$lines
            if (content[1] != "") {
              content <- c("", content)
            }
            if (tail(content, 1) != "") {
              content <- c(content, "")
            }
            writeLines(paste0(line, ": ", paste0(content, collapse = "\n")))
            .append_content(line, content)
            Sys.sleep(.3)
          }
        })
    }
    .C("nvimcom_msg_to_nvim", "ToggleWinheight()", PACKAGE = "nvimcom")
    .C("nvimcom_msg_to_nvim", "ShowNotifyWindow(['Appending contents.'])", PACKAGE = "nvimcom")
    if (target) {
      meta <- x$dataTarget
      anno <- matchMdContent(x$anno, x$target)
      fun_append(meta, anno)
    }
    if (result) {
      meta <- x$dataResult
      anno <- matchMdContent(x$anno, "# 总结内容")
      fun_append(meta, anno)
    }
    return(x)
  })

.append_content <- function(line, contents) {
  contents <- paste0(shQuote(contents), collapse = ", ")
  contents <- paste0("[", contents, "]")
  arg <- paste0("append(", line, ", ", contents, ")")
  .C("nvimcom_msg_to_nvim", arg, PACKAGE = "nvimcom")
}

matchMdContent <- function(data, target, only.title = F, mode = c("all", "description")) {
  mode <- match.arg(mode)
  nTarget <- which(data$lines == target)
  if (length(nTarget)) {
    if (length(nTarget) > 1) {
      stop("Too many match: ", paste0(nTarget, collapse = ", "))
    }
  } else {
    return(NULL)
  }
  if (mode == "all") {
    level <- nchar(strx(target, "^[#]+"))
  } else {
    level <- 100L
  }
  end <- 0L
  for (i in (nTarget + 1L):nrow(data)) {
    if (data$level[i] <= level && data$level[i] != -1L) {
      end <- i - 1L
      break
    }
  }
  if (!end) {
    end <- nrow(data)
  }
  dataTarget <- dplyr::slice(data, nTarget:end)
  if (only.title) {
    dplyr::filter(dataTarget, isTitle)
  } else {
    dataTarget
  }
}

get_clipboard <- function() {
  if (.Platform$OS.type != "unix") {
    stop("The value got TRUE: `.Platform$OS.type != 'unix'`")
  }
  system("xsel -b -o", intern = TRUE)
}

gett <- function(obj) {
  if (.Platform$OS.type != "unix") {
    stop("The value got TRUE: `.Platform$OS.type != 'unix'`")
  }
  if (length(obj) > 1) {
    obj <- paste0(obj, collapse = "\n")
  }
  system(paste("echo", shQuote(obj), "| xsel -b -i"))
}

gptMessages <- function() {
  list(heading = readLines(file.path(.expath, "gpt", "analysis_workflow_query_heading.txt")))
}

get_md_titles <- function(file = "index.Rmd", lines = NULL) {
  if (is.null(lines)) {
    lines <- readLines(file)
  }
  data <- tibble::tibble(lines = lines, number = seq_along(lines))
  inchunk <- F
  data$isInChunk <- vapply(lines, FUN.VALUE = integer(1), USE.NAMES = F,
    function(x) {
      if (grpl(x, "^```")) {
        inchunk <<- !inchunk
        return(-1L)
      } else {
        if (inchunk) {
          return(1L)
        } else {
          return(0L)
        }
      }
    })
  data <- dplyr::mutate(data,
    level = ifelse(isInChunk == 0L, nchar(strx(lines, "^[#]+")), -1L),
    level = ifelse(is.na(level), -1L, level),
    isTitle = level >= 1)
  data
}
