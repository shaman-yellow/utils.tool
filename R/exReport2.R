# ==========================================================================
# for conversion of 'report'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

as_report.rough <- function(lines) {
  yaml.pos <- grep("^---", lines)
  yaml.pos <- 1:(yaml.pos[2])
  yaml <- lines[ yaml.pos ]
  lines <- lines[-yaml.pos]
  new_report(lines, yaml = yaml)
}

write_biocStyle <- function(
  report, savename, title, change_include_fun = "inclu.fig",
  bioyml = readLines(paste0(.expath, "/", "biocstyle.yml")),
  origin_include_fun = "knitr::include_graphics")
{
  if (length(x <- grep("^title:", bioyml)) > 0) {
    bioyml <- bioyml[-x]
  }
  if (!grepl("^title: ", title)) {
    title <- paste0("title: ", title)
  }
  bioyml <- c(title, bioyml)
  yaml(report) <- bioyml
  lines <- call_command(report)
  if (!is.null(change_include_fun)) {
    lines <- gsub(origin_include_fun, change_include_fun, lines)
  }
  writeLines(lines, savename)
}

write_thesisDocx <- function(report, savename, title,
  change_include_fun = "inclu.fig",
  yml = readLines(paste0(.expath, "/", "ch_thesis.yml")),
  origin_include_fun = "knitr::include_graphics")
{
  write_biocStyle(report, savename, title, change_include_fun, yml, origin_include_fun)
}

write_articlePdf <- function(report, savename, title,
  change_include_fun = NULL,
  yml = readLines(paste0(.expath, "/", "articleWithCode.yml")),
  origin_include_fun = "knitr::include_graphics")
{
  write_biocStyle(report, savename, title, change_include_fun, yml, origin_include_fun)
}

kable_less <- function(x, ...) {
  if (nrow(x) > 50) {
    x <- head(x, n = 25)
  }
  knitr::kable(x, ...)
}
