# ==========================================================================
# for conversion of 'report'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

exclude_yaml <- function(lines) {
  yaml.pos <- grep("^---$", lines)
  yaml.pos <- 1:(yaml.pos[2])
  lines[-yaml.pos]
}

write_biocStyle <- function(
  report, savename, title, bioyml = file.path(.expath, "biocstyle.yml"),
  bib = NULL, ...)
{
  bioyml <- readLines(bioyml)
  if (length(x <- grep("^title:", bioyml))) {
    bioyml <- bioyml[-x]
  }
  if (!is.null(title)) {
    if (!grepl("^title: ", title)) {
      title <- paste0("title: ", title)
    }
    bioyml <- c(title, bioyml)
  }
  if (!is.null(bib)) {
    bioyml <- gsub("(bibliography:).*", paste0("\\1 ", bib), bioyml)
  }
  lines <- exclude_yaml(readLines(report))
  if (grpl(bioyml[1], "---")) {
    lines <- c(bioyml, "", lines)
  } else {
    lines <- c("---", bioyml, "---", "", lines)
  }
  writeLines(lines, savename)
  rmarkdown::render(savename)
}

write_thesisDocx <- function(report, savename, title,
  yml = file.path(.expath, "ch_thesis.yml"), ...)
{
  write_biocStyle(report, savename, title, yml, ...)
}

write_thesisDocxEn <- function(report, savename, title,
  yml = file.path(.expath, "en_thesis.yml"), ...)
{
  write_biocStyle(report, savename, title, yml, ...)
}

write_articlePdf <- function(report, savename = "output.Rmd", title = "",
  yml = file.path(.expath, "articleWithCode.yml"), ...)
{
  write_biocStyle(report, savename, title, yml, ...)
}

kable_less <- function(x, ...) {
  if (nrow(x) > 50) {
    x <- head(x, n = 25)
  }
  knitr::kable(x, ...)
}
