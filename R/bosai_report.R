
formatName.bosai <- function(file, info) {
  if (missing(info)) {
    info <- get("info", envir = parent.frame(1))
  }
  nfile <- paste0(info$id, "-", info$client, "-",
    info$type, "-", info$title, "-", format.Date(info$finish, "%Y.%m.%d"), ".docx")
  file.copy(file, nfile)
  browseURL(nfile)
  tools::file_path_sans_ext(nfile)
}

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

write_bosaiDocx <- function(report, savename, title,
  change_include_fun = "inclu.fig",
  yml = readLines(file.path(.expath, "bosai.yml")),
  origin_include_fun = "knitr::include_graphics", ...)
{
  write_biocStyle(report, savename, title, change_include_fun, yml, origin_include_fun, ...)
}

order_publish.bosai <- function(file = "index.Rmd", output = "output.Rmd", title = "",
  funs = list(write_bosaiDocx))
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
      content <- get_tex_cover.bosai_idea()
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.bosai_idea(env)
    }
  } else {
    if (knitr::is_latex_output()) {
      content <- get_tex_cover.bosai_analysis()
    } else if (knitr::pandoc_to("docx")) {
      content <- get_docx_cover.bosai_analysis(env)
    }
  }
  content <- vapply(content, glue::glue, character(1),
    .open = "[[[", .close = "]]]", .envir = env)
  writeLines(content)
}

get_docx_cover.bosai_idea <- function(env) {
  items <- c(
    "  研究方向：", "[[[title]]]",
    "  委托人：", "[[[info$client]]]",
    "  受托人：", "杭州铂赛生物科技有限公司"
  )
  .get_docx_coverStyle.bosai("生物医药合作项目开发", items, 6, 3, env = env)
}

get_docx_cover.bosai_analysis <- function(env) {
  items <- c(
    "  项目标题：", "[[[title]]]",
    "  单    号：", "[[[info$id]]]",
    "  分析人员：", "[[[author]]]",
    "  分析类型：", "[[[info$type]]]",
    "  委 托 人：", "[[[info$client]]]",
    "  受 托 人：", "杭州铂赛生物科技有限公司"
  )
  .get_docx_coverStyle.bosai("生信分析报告", items, 4, 3, env = env)
}

.get_docx_coverStyle.bosai <- function(header, items, header.before = 4, header.after = 3,
  env, lineNChar = 45)
{
  fp_header <- officer::fp_par("center", line_spacing = 2)
  fp_items <- officer::fp_par(line_spacing = 2.5)
  fp_tlarge <- officer::fp_text(font.size = 26, font.family = "Microsoft YaHei", bold = T)
  fp_tsmall <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = T)
  fp_tsmall_und <- officer::fp_text(font.size = 14, font.family = "SimSun", bold = T, underlined = T)
  blank <- officer::fpar()
  par_title <- list(officer::fpar(officer::ftext(header, fp_tlarge), fp_p = fp_header))
  N <- 0L
  NEnd <- length(items) / 2
  par_items <- mapply(seq(1, length(items), 2), seq(2, length(items), 2), SIMPLIFY = F,
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

assis_docx_par <- function(x) {
  c("```{=openxml}", officer::to_wml(x), "```")
}

assis_docx_img <- function(x) {
  paste0("`", officer::to_wml(x), "`{=openxml}")
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
