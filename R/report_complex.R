# ==========================================================================
# complex
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


order_publish.complex <- function(file = "index.Rmd", output = "output.Rmd", title = "",
  MyInstitution = getOption("MyInstitution"), ...)
{
  if (MyInstitution == "complex") {
    fun <- get_fun("write_complexDocx")
  }
  browseURL(fun(file, output, title), "wps")
}

write_complexDocx <- function(report, savename, title,
  yml = file.path(.expath, "complex.yml"), ...)
{
  write_biocStyle(report, savename, title, yml, ...)
}

formatName.complex <- function(file, info) {
  if (missing(info)) {
    info <- get("info", envir = parent.frame(1))
  }
  ext <- tools::file_ext(file)
  nfile <- paste0(info$id, "-", info$client, "-",
    info$type, "-", info$title, "-", format.Date(info$finish, "%Y.%m.%d"), ".", ext)
  file.copy(file, nfile, TRUE)
  gett(dirname(normalizePath(nfile)))
  tools::file_path_sans_ext(nfile)
}


set_cover.complex <- function(info, title = info$title, author = "黄礼闯", date = Sys.Date(),
  coverpage = NULL, institution = "铂赛")
{
  env <- environment()
  if (is.null(coverpage)) {
    if (info$type == "思路设计") {
      coverpage <- .prefix("cover_page_idea.pdf")
    } else {
      coverpage <- .prefix("cover_page_analysis.pdf")
    }
  }
  content <- NULL
  if (info$type == "思路设计") {
    if (knitr::is_latex_output()) {
      stop('knitr::is_latex_output()')
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.complex_idea(env)
    }
  } else {
    if (knitr::is_latex_output()) {
      stop('knitr::is_latex_output()')
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.complex_analysis(env)
    }
  }
  content <- vapply(content, glue::glue, character(1),
    .open = "[[[", .close = "]]]", .envir = env)
  writeLines(content)
}


get_docx_cover.complex_idea <- function(env) {
  items <- c(
    "  研究方向：", "[[[title]]]",
    "  委托人：", "[[[info$client]]]",
    "  受托人：", "..."
  )
  .get_docx_coverStyle.complex("生物医药合作项目开发", items, 6, 3, env = env)
}

get_docx_cover.complex_analysis <- function(env) {
  items <- c(
    "  项目标题：", "[[[title]]]",
    "  单    号：", "[[[info$id]]]",
    "  分析人员：", "[[[author]]]",
    "  分析类型：", "[[[info$type]]]",
    "  委 托 人：", "[[[info$client]]]",
    "  受 托 人：", "..."
  )
  .get_docx_coverStyle.complex("生信分析报告", items, 4, 3, env = env)
}

.get_docx_coverStyle.complex <- function(header, items, header.before = 4, header.after = 3,
  env, lineNChar = 45)
{
  fp_header <- officer::fp_par("center", line_spacing = 2)
  fp_items <- officer::fp_par(line_spacing = 2.5)
  fp_tlarge <- officer::fp_text(font.size = 26, font.family = "Microsoft YaHei", bold = TRUE)
  fp_tsmall <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = TRUE)
  fp_tsmall_und <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = TRUE, underlined = TRUE)
  blank <- officer::fpar()
  par_title <- list(officer::fpar(officer::ftext(header, fp_tlarge), fp_p = fp_header))
  N <- 0L
  NEnd <- length(items) / 2
  par_items <- mapply(seq(1, length(items), 2), seq(2, length(items), 2), SIMPLIFY = FALSE,
    FUN = function(n1, n2) {
      N <<- N + 1L
      query <- officer::ftext(items[n1], fp_tsmall)
      answer <- glue::glue(items[n2], .open = "[[[", .close = "]]]", .envir = env)
      answer <- split_text_by_width(answer, lineNChar - 5)
      if (length(answer) >= 2) {
        if (nchar(tail(answer, n = 1)) <= 2) {
          answer[ length(answer) - 1 ] %<>% paste0(., tail(answer, n = 1))
          answer <- answer[ -length(answer) ]
        }  
      }
      answer <- stringr::str_pad(answer, lineNChar, "both")
      if (length(answer) >= 2) {
        lineBreak <- officer::run_linebreak()
        anlist <- list()
        length(anlist) <- length(answer) * 3 - 2
        for (i in seq_along(answer)) {
          if (i == length(answer)) {
            postfix <- if (N == NEnd) "." else ";"
          } else {
            postfix <- ""
          }
          anlist[[ 3 * i - 2 ]] <- officer::ftext(
            paste0(answer[i], postfix), fp_tsmall_und
          )
          if (i != length(answer)) {
            anlist[[ 3 * i - 1 ]] <- lineBreak
            anlist[[ 3 * i ]] <- officer::ftext(
              stringr::str_pad("", nchar(items[1]) * 2 - 2), fp_tsmall
            )
          }
        }
        answer <- anlist
      } else {
        answer <- paste0(answer, if (N == NEnd) "." else ";")
        answer <- list(officer::ftext(answer, fp_tsmall_und))
      }
      args <- c(list(query), answer, list(fp_p = fp_items))
      do.call(officer::fpar, args)
    })
  content <- c(
    rep(list(blank), header.before), par_title,
    rep(list(blank), header.after), par_items
  )
  content <- unlist(lapply(content, assis_docx_par))
  c("", content, "", "\\newpage", "")
}

