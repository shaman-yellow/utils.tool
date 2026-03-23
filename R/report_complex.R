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

order_packaging.complex <- function(target = "output.pdf",
  register = autoRegisters, idname = gidn(), external_file = NULL,
  ..., extras = NULL,
  codes = c(.getRICodes(), s(target, "_out.docx", ".Rmd")))
{
  extras <- c(extras, codes)
  order_packaging(
    target, register, idname, extras = extras, external_file = external_file, ...
  )
}

.getRICodes <- function(file_install = normalizePath("~/.vim/others/install.R"),
  dir_package = normalizePath("~/utils.tool"), package_name = basename(dir_package))
{
  archive_package <- paste0(package_name, ".tar.gz")
  cdRun(
    "git ls-files -z | xargs -0 tar -czf ", archive_package,
    path = dir_package
  )
  file.copy(file_install, ".", TRUE)
  file.copy(
    file.path(dir_package, archive_package), ".", TRUE
  )
  return(c(basename(file_install), archive_package))
}

write_complexDocx <- function(report, savename, title,
  yml = file.path(.expath, "complex.yml"), ...)
{
  play_overture(report, savename, title, yml, ...)
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
  coverpage = NULL, institution = "...")
{
  if (missing(info)) {
    file_info <- ".items_analysis.rds"
    if (file.exists(file_info)) {
      info <- readRDS(file_info)
    } else {
      stop('file.exists(file_info), and missing `info`...')
    }
  }
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
      content <- get_docx_cover.complex_analysis(
        info$title, info, env
      )
    }
  }
  content <- vapply(content, glue::glue, character(1),
    .open = "[[[", .close = "]]]", .envir = env)
  writeLines(content)
}


# get_docx_cover.complex_idea <- function(env) {
#   items <- c(
#     "  研究方向：", "[[[title]]]",
#     "  委托人：", "[[[info$client]]]",
#     "  受托人：", "..."
#   )
#   .get_docx_coverStyle.complex("生物医药合作项目开发", items, 6, 3, env = env)
# }

get_docx_cover.complex_analysis <- function(title, items, env) {
  .get_docx_coverStyle.complex(title, items, env = env)
}

.get_docx_coverStyle.complex <- function(title, items, header.before = 8, header.after = 12,
  env, lineNChar = 45)
{
  par_title <- officer::fpar(
    officer::ftext(title,
      officer::fp_text(
        eastasia.family = "SimSun", font.family = "Times New Roman", font.size = 18, bold = TRUE
        )), fp_p = officer::fp_par(line_spacing = 1.5, text.align = "center")
  )
  fp_text_id <- officer::fp_text(
    bold = TRUE, font.size = 16, font.family = "Times New Roman", eastasia.family = "SimSun"
  )
  par_id <- officer::fpar(
    officer::ftext(glue::glue("项目编号："), fp_text_id),
    officer::run_linebreak(),
    officer::ftext(glue::glue("Project-{s(items$id, '[a-z]+', '')}"), fp_text_id),
    fp_p = officer::fp_par(text.align = "center", line_spacing = 1.5)
  )
  blank <- officer::fpar()
  content <- c(
    rep(list(blank), header.before), list(par_title),
    rep(list(blank), header.after), list(par_id)
  )
  content <- unlist(lapply(content, assis_docx_par))
  c("", content, "", "\\newpage", "")
}

.md_p_significant <- "\\\\*, P < 0.05; \\\\*\\\\*, p < 0.01; \\\\*\\\\*\\\\*, p < 0.001"

