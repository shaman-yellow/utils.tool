
summary_month.bosai <- function(data, time = Sys.Date(),
  back_dir = .prefix("summary", "wd"),
  .time = if (is(time, "Date")) time else as.Date(time, tryFormats = "%Y-%m-%d"), 
  upd = FALSE, ...)
{
  data <- dplyr::filter(data, finish <= !!time)
  data <- register_data(data, file.path(back_dir, "month_records.rds"))
  summary_week.bosai(
    time, orders = data, templ_dir = .prefix("summary/"), templ_file = "summary_month.xlsx",
    dir_prefix = "summary_month_", dir_suffix = paste0(
      lubridate::year(time), "_", lubridate::month(time)
    ), mode = "month"
  )
  register_data <- upd
  return(data)
}


summary_assess.bosai <- function(data, time = Sys.Date(),
  back_dir = .prefix("summary", "wd"),
  .time = if (is(time, "Date")) time else as.Date(time, tryFormats = "%Y-%m-%d"))
{
  data <- dplyr::filter(data, finish <= !!time)
  data <- register_data(data, file.path(back_dir, "records.rds"))
  export <- file.path(
    back_dir, paste0("summary_", Sys.Date(), ".csv")
  )
  data <- dplyr::mutate(data, cost = ifelse(type == "思路设计", 3, NA))
  data <- dplyr::select(data, type, id, title, cost, finish)
  data <- dplyr::mutate(data, name = paste0(id, ", ", title))
  data <- dplyr::select(data, type, name, cost, finish)
  data.table::fwrite(data, export)
  register_data <- TRUE
  browseURL(normalizePath(export))
}

backup_jobs.bosai <- function(data, time = Sys.Date(),
  month = lubridate::month(.time - 12),
  year = lubridate::year(.time - 12), back_dir = .prefix("backup", "db"),
  .time = if (is(time, "Date")) time else as.Date(time, tryFormats = "%Y-%m-%d"))
{
  dir.create(back_dir, FALSE)
  file_records <- file.path(back_dir, "records.rds")
  data <- dplyr::filter(data, finish < !!time)
  if (file.exists(file_records)) {
    records <- readRDS(file_records)
    data <- dplyr::anti_join(data, records)
  } else {
    records <- tibble::tibble()
    data <- data
  }
  data <- dplyr::filter(data, finish < time)
  if (!nrow(data)) {
    message("No jobs need to backup.")
    print(records)
    return()
  }
  thedir <- .thedir_job(month, year, back_dir)
  data$backup_dir <- thedir
  num <- 0L
  dir.create(thedir, FALSE)
  res <- mapply(data$.dir, data$id, data$client, data$type, data$title, data$save,
     FUN = function(dir, id, client, type, title, rds_item) {
       pattern <- paste0(".*", id, ".*", client, ".*", type, ".*", title, ".*")
       files <- list.files(dir, pattern, full.names = TRUE)
       if (!length(files)) {
         message("No zip files found in dir: ", pattern)
         isThat <- usethis::ui_yeah("Select other files?")
         if (isThat) {
           files <- list.files(dir, "*.zip|*.docx|*.pdf", full.names = TRUE)
           files <- files[ menu(files, title = "All available files") ]
           if (!length(files)) {
             return()
           }
         } else {
           return()
         }
       }
       script <- paste0(gs(rds_item, "^\\.items_|\\.rds$", ""), ".Rmd")
       script <- file.path(dir, script)
       if (!file.exists(script)) {
         stop(glue::glue("No such file: {script}"))
       }
       files <- c(files, script)
       if (length(files)) {
         num <<- num + 1L
         to <- paste0(thedir, "/", num, "_", basename(tools::file_path_sans_ext(files[1])))
         .cp_job(c(files, script), to)
       }
     })
  records <- dplyr::bind_rows(records, data)
  saveRDS(records, file_records)
  unlist(res)
  browseURL(thedir)
}

register_data <- function(data, rds_records, no_to_stop = TRUE,
  fun_mutate = function(x) dplyr::mutate(x, summary_month = Sys.Date(), .before = 1))
{
  if (file.exists(rds_records)) {
    records <- readRDS(rds_records)
    data <- dplyr::anti_join(data, records)
  } else {
    records <- tibble::tibble()
  }
  if (!nrow(data)) {
    print(records)
    if (no_to_stop) {
      stop("No data need to backup.")
    }
  }
  records <- dplyr::bind_rows(records, fun_mutate(data))
  env <- parent.frame(1)
  withr::defer({
    isThat <- try(get("register_data", envir = env), TRUE)
    if (!inherits(isThat, "try-error") && is.logical(isThat) && isThat) {
      saveRDS(records, rds_records)
    }
  }, env)
  return(data)
}

formatName.bosai <- function(file, info) {
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

summary_week.bosai <- function(
  time = Sys.Date(),
  orders = get_orders(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  week = lubridate::week(time),
  path = .prefix(),
  templ_dir = .prefix("summary/"),
  templ_file = "summary.xlsx",
  rm = FALSE, must_include = "__must_include__", dir_prefix = "summary_", 
  dir_suffix = week, mode = c("week", "month"))
{
  mode <- match.arg(mode)
  dir <- file.path(path, paste0(dir_prefix, dir_suffix))
  targets <- c(ass = templ_file)
  if (rm) {
    unlink(list.files(dir, full.names = TRUE, all.files = TRUE, recursive = TRUE), TRUE, TRUE)
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
  if (mode == "week") {
    orders <- dplyr::filter(orders,
      is.na(lubridate::week(finish)) | (lubridate::week(finish) == !!week) |
        id %in% must_include
    )
  }
  ## modify the assess table
  data_ass <- dplyr::mutate(orders, expect = finish)
  data_ass <- dplyr::select(data_ass,
    id, client, type, title, start, end, expect, finish, note
  )
  fun <- function(wb, data) {
    openxlsx2::wb_add_data(wb, 1, data, col_names = FALSE,
      dims = do.call(openxlsx2::wb_dims, pos.data_ass), na.strings = "")
  }
  pos.data_ass <- list(3, 1)
  wb <- fun(wb, dplyr::filter(data_ass, !is.na(finish) | id %in% must_include))
  if (mode == "week") {
    pos.data_ass <- list(15, 1)
    wb <- fun(wb, dplyr::filter(data_ass, is.na(finish)))
  }
  openxlsx2::wb_save(wb, targets[[ "ass" ]])
  ## check
  gett(dirname(normalizePath(targets[[ "ass" ]])))
  browseURL(normalizePath(targets[[ "ass" ]]))
}

td <- function(character_date) {
  if (grpl(character_date, "-")) {
    as.Date(character_date, tryFormats = "%Y-%m-%d")
  } else {
    as.Date(character_date, tryFormats = "%Y%m%d")
  }
}

write_bosaiDocx <- function(report, savename, title,
  yml = file.path(.expath, "bosai.yml"), ...)
{
  write_biocStyle(report, savename, title, yml, ...)
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

extract_anno <- function(file, type = tools::file_ext(file), ...) {
  # fun <- function(file) {
  #   name <- basename(file)
  #   if (name != make.names(name)) {
  #     file.rename(file, file <- file.path(dirname(file), name))
  #   }
  #   return(file)
  # }
  if (type == "pdf") {
    cdRun(glue::glue("pdfannots {shQuote(file)} -o {save}"))
  } else if (type == "docx") {
    tmp <- tempfile("anno", fileext = ".md")
    cli::cli_alert_info(glue::glue("tempfile: {tmp}"))
    # file <- fun(file)
    cdRun(glue::glue("pandoc {shQuote(file)} -o {tmp} --track-changes=all"))
    postModify_annoDocx(tmp, save = file.path(dirname(file), "comment_reply.md"), ...)
  }
}

postModify_annoDocx <- function(file, prefix = "comment",
  save = file.path(dirname(file), paste0(prefix, "_reply.md")), add_reply = TRUE)
{
  if (file.exists(save)) {
    if (!usethis::ui_yeah('file.exists(save), overwrite that?')) {
      return(message("Do nothing."))
    }
  }
  lines <- readLines(file)
  chapter <- sep_list(lines, "^[-]+$|^[=]+$", -1)
  names(chapter) <- vapply(chapter, function(x) x[1], character(1))
  chapter <- lapply(chapter, function(text) {
    text <- paste(text, collapse = " ")
    text <- stringr::str_extract_all(text,
      "\\[[^\\]]*?\\]\\{\\.comment-start.*?\\}.*?\\[[^\\]]*?\\]\\{\\.comment-end.*?\\}")
    text <- unlist(text)
    text <- text[ !is.na(text) ]
    text <- gs(text, "¶", ", ")
  })
  chapter <- lst_clear0(chapter)
  chapter <- lapply(chapter, function(x) {
    content <- gs(x, ".*.comment-start[^}]*\\}(.*)\\[\\]\\{.comment-end.*", "\\1")
    comment <- strx(x, "(?<=\\[)[^\\]]*(?=\\])")
    meta <- strx(x, "(?<=\\{)[^}]*(?=\\})")
    meta <- gs(meta, '.*author="([^\"]*?)".*', "\\1")
    glue::glue(paste0("* Content: {content}\n* {meta} Comment: {comment}"))
  })
  if (add_reply) {
    chapter <- lapply(chapter,
      function(x) {
        unlist(lapply(x, function(i) c(i, "", "> Reply: ", "", "---------------", "")))
      })
  }
  lines <- lapply(seq_along(chapter),
    function(n) {
      head <- names(chapter)[n]
      level <- length(stringr::str_extract_all(head, "\\.")[[1]]) + 1
      prefix <- paste0(rep("#", level), collapse = "")
      c(paste0(prefix, " ", head), "", chapter[[ n ]], "")
    })
  writeLines(unlist(lines), save)
}
