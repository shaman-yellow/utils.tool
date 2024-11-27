
summary_week.bosai <- function(
  time = Sys.Date(),
  orders = get_orders(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  week = lubridate::week(time),
  path = .prefix(),
  templ_dir = .prefix("summary/"),
  rm = F)
{
  dir <- file.path(path, paste0("summary_", week))
  targets <- c(ass = "summary.xlsx")
  if (rm) {
    unlink(list.files(dir, full.names = T, all.files = T, recursive = T), T, T)
    dir.create(dir)
    file.copy(file.path(templ_dir, targets), dir)
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
    file.copy(file.path(templ_dir, targets), dir)
  }
  targets[] <- file.path(dir, targets)
  wb <- openxlsx2::wb_load(targets[[ "ass" ]])
  ## prepare orders (this week or next week)
  orders <- dplyr::filter(orders,
    is.na(lubridate::week(finish)) | (lubridate::week(finish) == !!week)
  )
  ## modify the assess table
  data_ass <- dplyr::mutate(orders, expect = finish)
  data_ass <- dplyr::select(data_ass,
    id, client, type, title, start, end, expect, finish, note
  )
  fun <- function(wb, data) {
    openxlsx2::wb_add_data(wb, 1, data, col_names = F,
      dims = do.call(openxlsx2::wb_dims, pos.data_ass), na.strings = "")
  }
  pos.data_ass <- list(3, 1)
  wb <- fun(wb, dplyr::filter(data_ass, !is.na(finish)))
  pos.data_ass <- list(15, 1)
  wb <- fun(wb, dplyr::filter(data_ass, is.na(finish)))
  openxlsx2::wb_save(wb, targets[[ "ass" ]])
  ## check
  browseURL(normalizePath(targets[[ "ass" ]]))
}

td <- function(character_date) {
  as.Date(character_date, tryFormats = "%Y%m%d")
}

order_publish.bosai <- function(file = "index.Rmd", output = "output.Rmd", title = "",
  funs = list(write_thesisDocxEn))
{
  n <- 0L
  lapply(funs,
    function(fun) {
      n <<- n + 1L
      res <- fun(file, output, title)
      if (n == length(funs)) {
        browseURL(res)
      }
    })
  invisible()
}

set_cover.bosai <- function(info, title = info$title, author = "黄礼闯", date = Sys.Date(),
  coverpage = NULL, institution = "铂赛")
{
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
      content <- get_tex_cover.bosai_idea()
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.bosai_idea()
    }
  } else {
    if (knitr::is_latex_output()) {
      content <- get_tex_cover.bosai_analysis()
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.bosai_analysis()
    }
  }
  env <- environment()
  content <- vapply(content, glue::glue, character(1),
    .open = "[[[", .close = "]]]", .envir = env)
  writeLines(content)
}

get_docx_cover.bosai_idea <- function() {
  fp_pcenter <- officer::fp_par("center", line_spacing = 2)
  fp_tlarge <- officer::fp_text(font.size = 26, font.family = "SimHei", bold = T)
  fp_tsmall <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = T)
  blank <- officer::fpar()
  content <- c(
    rep(list(blank), 9),
    list(officer::fpar(officer::ftext("生物医药合作项目开发", fp_tlarge), fp_p = fp_pcenter)),
    rep(list(blank), 5),
    list(officer::fpar(officer::ftext("研究方向：[[[title]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("委托人：[[[info$client]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("受托人：杭州铂赛生物科技有限公司", fp_tsmall), fp_p = fp_pcenter))
  )
  content <- unlist(lapply(content, assis_docx_officer))
  c("", content, "", "\\newpage", "")
}

get_docx_cover.bosai_analysis <- function() {
  fp_pcenter <- officer::fp_par("center", line_spacing = 2)
  fp_tlarge <- officer::fp_text(font.size = 26, font.family = "SimHei", bold = T)
  fp_tsmall <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = T)
  blank <- officer::fpar()
  content <- c(
    rep(list(blank), 5),
    list(officer::fpar(officer::ftext("生信分析报告", fp_tlarge), fp_p = fp_pcenter)),
    rep(list(blank), 5),
    list(officer::fpar(officer::ftext("项目标题：[[[title]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("单号：[[[info$id]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("分析人员：[[[author]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("分析类型：[[[info$type]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("委托人：[[[info$client]]]", fp_tsmall), fp_p = fp_pcenter)),
    list(officer::fpar(officer::ftext("受托人：杭州铂赛生物科技有限公司", fp_tsmall), fp_p = fp_pcenter))
  )
  content <- unlist(lapply(content, assis_docx_officer))
  c("", content, "", "\\newpage", "")
}

assis_docx_officer <- function(x) {
  c("```{=openxml}", officer::to_wml(x), "```")
}

get_tex_cover.bosai_analysis <- function() {
  strwrap(paste0("\\begin{titlepage}
      \\newgeometry{top=6.5cm}
      \\ThisCenterWallPaper{1.12}{[[[coverpage]]]}
      \\begin{center}
      \\textbf{\\huge [[[title]]]}
      \\vspace{4em}
      \\begin{textblock}{10}(3,4.85)
      \\Large \\textbf{\\textcolor{black}{[[[info$id]]]}}
      \\end{textblock}
      \\begin{textblock}{10}(3,5.8)
      \\Large \\textbf{\\textcolor{black}{[[[author]]]}}
      \\end{textblock}
      \\begin{textblock}{10}(3,6.75)
      \\Large \\textbf{\\textcolor{black}{[[[info$type]]]}}
      \\end{textblock}
      \\begin{textblock}{10}(3,7.7)
      \\Large \\textbf{\\textcolor{black}{[[[info$client]]]}}
      \\end{textblock}
      \\end{center}
      \\end{titlepage}
      \\restoregeometry
      "
      ), 50)
}

get_tex_cover.bosai_idea <- function() {
  strwrap(paste0("\\begin{titlepage}
      \\newgeometry{top=7.5cm}
      \\ThisCenterWallPaper{1.12}{[[[coverpage]]]}
      \\begin{center}
      \\textbf{\\huge [[[title]]]}
      \\vspace{4em}
      \\begin{textblock}{10}(3.2,9.25)
      \\huge \\textbf{\\textcolor{black}{[[[date]]]}}
      \\end{textblock}
      \\end{center}
      \\end{titlepage}
      \\restoregeometry
      "
      ), 50) 
}
