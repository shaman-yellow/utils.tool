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
      content <- get_docx_cover.complex_analysis(info$id, env)
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

get_docx_cover.complex_analysis <- function(title, env) {
  .get_docx_coverStyle.complex(title, 10, 3, env = env)
}

.get_docx_coverStyle.complex <- function(header, header.before = 10, header.after = 3,
  env, lineNChar = 45)
{
  header <- paste0("项目编号: ", header)
  fp_header <- officer::fp_par("center", line_spacing = 2)
  fp_tlarge <- officer::fp_text(font.size = 25, font.family = "SimSun", bold = TRUE)
  par_title <- list(officer::fpar(officer::ftext(header, fp_tlarge), fp_p = fp_header))
  blank <- officer::fpar()
  content <- c(
    rep(list(blank), header.before), par_title
  )
  content <- unlist(lapply(content, assis_docx_par))
  c("", content, "", "\\newpage", "")
}

