# ==========================================================================
# for conversion of 'report'
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

exclude_yaml <- function(lines) {
  yaml.pos <- grep("^---$", lines)
  yaml.pos <- 1:(yaml.pos[2])
  lines[-yaml.pos]
}

.set_overture_on <- function() {
  options("__OVERTURE__" = TRUE)
}

.set_overture_off <- function() {
  options("__OVERTURE__" = FALSE)
}

.is_in_overture <- function() {
  getOption("__OVERTURE__", FALSE)
}

play_overture <- function(
  report, savename, title, bioyml = file.path(.expath, "biocstyle.yml"),
  bib = NULL, dir_jobsComments = "R_jobsComments", ...)
{
  .set_overture_on()
  on.exit(.set_overture_off())
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
  if (dir.exists(dir_jobsComments)) {
    upd_all_job_comments(dir_jobsComments)
  }
  lines <- parse_overture(lines, ...)
  # location
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

parse_overture <- function(lines, env = .GlobalEnv, ...)
{
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
  ovRes <- lapply(seq_along(ovlnums),
    function(n) {
      lnums <- ovlnums[[ as.character(n) ]]
      fields <- c(lnums, max(lnums) + 1L)
      lnumCodes <- lnums[ lnums %in% whichInChunk ]
      codes <- lines[lnumCodes]
      asis <- eval(parse(text = codes), envir = env)
      if (is.null(asis)) {
        asis <- ""
      } else if (!is.character(asis)) {
        stop('!is.character(asis), the eval results of codes in Overture should be character.')
      }
      namel(fields, asis, codes)
    })
  append_chunk_after_overture <- list()
  if (!is.null(setLevel <- getOption("autor_locate_dir_after_final_overture", 2L))) {
    if (!is.integer(setLevel)) {
      stop('!is.integer(setLevel).')
    }
    rds_db <- ".job_locate_in_script.rds"
    if (!file.exists(rds_db)) {
      # isOverture <- grpl()
      warning(
        crayon::red("autor_locate_dir_after_final_overture is not NULL, ",
          "but db file not found: .job_locate_in_script.rds")
      )
    } else {
      chLoc <- parse_chunk_location(lines)
      ovField <- lapply(ovlnums,
        function(ns) {
          n <- max(ns)
          chLoc$lineFields[[n]]
        })
      groupBy <- vapply(ovField, FUN.VALUE = character(1),
        function(field) {
          paste0(head(field, n = setLevel), collapse = "_")
        })
      if (setLevel >= length(ovField[[1]])) {
        slevel <- length(ovField[[1]])
      } else {
        slevel <- setLevel
      }
      smallestField <- lapply(ovField, function(x) x[slevel])
      # all overtures under the field
      ovGroupByField <- split(seq_along(ovField), groupBy)
      # the same order:
      whichLast <- vapply(ovGroupByField, max, integer(1))
      job_belong <- readRDS(rds_db)
      chunk_start <- '```{r eval = TRUE, echo = FALSE, results = "asis"}'
      chunk_end <- '```'
      lst_cmdAppend <- lapply(whichLast,
        function(n) {
          codes <- ovRes[[n]]$codes
          mayBeJobName <- unlist(stringr::str_extract_all(codes, "[a-zA-Z0-9_.]+"))
          mayBeJobName <- mayBeJobName[ !is.na(mayBeJobName) ]
          jobs <- mayBeJobName[ mayBeJobName %in% names(job_belong) ]
          if (!length(jobs)) {
            warning("No any 'job' matched in overture {n}, codes, skip:\n{bind(codes, co = '\n')}.")
            cmdAppend <- NULL
            dir <- NULL
          } else {
            dirs <- vapply(jobs, FUN.VALUE = character(1),
              function(job) {
                job_belong[[job]]$dir
              })
            if (length(unique(dirs)) > 1) {
              warning(glue::glue(crayon::red("A overture (n = {n}) belong to Multiple directory: {bind(dirs)}")))
              freq <- table(dirs)
              dir <- names(freq)[ which.max(freq) ]
            } else {
              dir <- dirs[[1]]
            }
            cmd <- glue::glue("autor_locate_dir('{dir}')")
            cmdAppend = c("", chunk_start, cmd, chunk_end, "")
          }
          list(dir = dir, cmdAppend = cmdAppend)
        })
      dirs_autor_ovLoc <- lapply(lst_cmdAppend, function(x) x$dir)
      # the same order: dirs_autor_ovLoc, ovCodes_and_locDir
      ovCodes_and_locDir <- lapply(
        seq_along(dirs_autor_ovLoc),
          function(n) {
            ovInGroups <- ovGroupByField[[n]]
            lapply(ovInGroups,
              function(ovNum) {
                codes <- ovRes[[ovNum]]$codes
                namel(codes, dir = dirs_autor_ovLoc[[n]], seq_of_overture = ovNum)
              })
          }
      )
      ovCodes_and_locDir <- unlist(ovCodes_and_locDir, recursive = FALSE, use.names = FALSE)
      order <- order(vapply(ovCodes_and_locDir, function(x) x$seq_of_overture, integer(1)))
      ovCodes_and_locDir <- ovCodes_and_locDir[ order ]
      # set global for overture usage.
      options(overture_codes_and_location = ovCodes_and_locDir)
      # prepare the appending
      cmdAppend <- lapply(lst_cmdAppend, function(x) x$cmdAppend)
      append_chunk_after_overture <- setNames(cmdAppend, unname(whichLast))
    }
  }
  lapply(rev(seq_along(ovRes)),
    function(n) {
      args <- ovRes[[n]]
      lines <<- lines[ -args$fields ]
      others <- args$asis
      if (!is.null(ex <- append_chunk_after_overture[[ as.character(n) ]])) {
        others <- c(others, ex)
      }
      lines <<- append(lines, others, min(args$fields) - 1L)
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
