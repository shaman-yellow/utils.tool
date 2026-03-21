# ==========================================================================
# for conversion of 'report'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

exclude_yaml <- function(lines) {
  yaml.pos <- grep("^---$", lines)
  yaml.pos <- 1:(yaml.pos[2])
  lines[-yaml.pos]
}

play_overture <- function(
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
  chunk_location <- parse_chunk_location(lines)
  options(chunk_location = chunk_location)
  autoRegisters_relocate <- new.env()
  options(autoRegisters_relocate = autoRegisters_relocate)
  # see `autor`
  options(autor_legend_env = new.env())
  if (grpl(bioyml[1], "---")) {
    lines <- c(bioyml, "", lines)
  } else {
    lines <- c("---", bioyml, "---", "", lines)
  }
  writeLines(lines, savename)
  rmarkdown::render(savename, clean = FALSE)
}

parse_overture <- function(lines, env = .GlobalEnv) {
  isQuos <- grpl(lines, "^```")
  inchunk <- FALSE
  isInChunk <- vapply(seq_along(lines), function(n) {
    if (isQuos[n]) {
      inchunk <<- !inchunk
    }
    inchunk
  }, logical(1))
  whichInChunk <- which(isInChunk & !isQuos)
  isOverture <- grpl(lines, "<!-- OVERTURE_[A-Z]+ -->")
  inOverture <- FALSE
  groupOv <- rep(0, length(lines))
  num <- 0L
  isInOverture <- vapply(seq_along(lines), FUN.VALUE = logical(1),
    function(n) {
      if (isOverture[n]) {
        if (!inOverture) {
          num <<- num + 1L
        }
        inOverture <<- !inOverture
      }
      if (inOverture) {
        groupOv[n] <<- num
      }
      inOverture
    })
  ovlnums <- split(seq_along(lines), groupOv)
  ovlnums <- ovlnums[ names(ovlnums) != "0" ]
  res <- lapply(seq_along(ovlnums),
    function(n) {
      lnums <- ovlnums[[ as.character(n) ]]
      fields <- c(lnums, max(lnums) + 1L)
      lnumCodes <- lnums[ lnums %in% whichInChunk ]
      res <- eval(parse(text = lines[lnumCodes]), envir = env)
      if (is.null(res)) {
        res <- ""
      } else if (!is.character(res)) {
        stop('!is.character(res), the eval results of codes in Overture should be character.')
      }
      namel(fields, res)
    })
  lapply(rev(seq_along(res)),
    function(n) {
      args <- res[[n]]
      lines <<- lines[ -args$fields ]
      lines <<- append(lines, args$res, min(args$fields) - 1L)
    })
  lines
}

parse_chunk_location <- function(lines) {
  isQuos <- grpl(lines, "^```")
  inchunk <- FALSE
  isInChunk <- vapply(seq_along(lines), function(n) {
    if (isQuos[n]) {
      inchunk <<- !inchunk
    }
    inchunk
  }, logical(1))
  allHeaders <- stringr::str_extract(lines, "^[#]+(?=\\s)")
  headerLevels <- nchar(allHeaders)
  headerLevels <- ifelse(isInChunk, NA, headerLevels)
  lowestLevel <- max(headerLevels, na.rm = TRUE)
  nowFields <- rep(0L, lowestLevel)
  lineFields <- rep(list(nowFields), length(lines))
  isHeaders <- !is.na(headerLevels)
  lastLevel <- 0L
  for (i in seq_along(headerLevels)) {
    if (isHeaders[i]) {
      nowFields[ headerLevels[i] ] <- nowFields[ headerLevels[i] ] + 1L
      if (headerLevels[i] < lastLevel) {
        nowFields[ lastLevel ] <- 0L
      }
      lastLevel <- headerLevels[i]
    }
    lineFields[[ i ]] <- nowFields
  }
  allChunks <- ifelse(isInChunk, lines, NA)
  chunkFields <- lineFields[ grpl(allChunks, "^```") ]
  allChunks <- sep_list(allChunks[ !is.na(allChunks) ], "^```", TRUE)
  if (length(chunkFields) != length(allChunks)) {
    stop("length(chunkFields) != length(allChunks)")
  }
  allChunks[-1] <- lapply(allChunks[-1], function(x) x[-1])
  chunkBelong <- vapply(chunkFields, function(level) {
    for (i in seq_along(lineFields)) {
      if (identical(level, lineFields[[i]])) {
        return(i)
      }
    }
  }, integer(1))
  chunkBelong <- lines[chunkBelong]
  list(lineFields = lineFields, allChunks = allChunks,
    chunkFields = chunkFields, chunkBelong = chunkBelong)
}

write_thesisDocx <- function(report, savename, title,
  yml = file.path(.expath, "ch_thesis.yml"), ...)
{
  play_overture(report, savename, title, yml, ...)
}

write_thesisDocxEn <- function(report, savename, title,
  yml = file.path(.expath, "en_thesis.yml"), ...)
{
  play_overture(report, savename, title, yml, ...)
}

write_articlePdf <- function(report, savename = "output.Rmd", title = "",
  yml = file.path(.expath, "articleWithCode.yml"), ...)
{
  play_overture(report, savename, title, yml, ...)
}

kable_less <- function(x, ...) {
  if (nrow(x) > 50) {
    x <- head(x, n = 25)
  }
  knitr::kable(x, ...)
}
