# ==========================================================================
# workflow of geo
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_geo <- setClass("job_geo", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    method = "GEO <https://www.ncbi.nlm.nih.gov/geo/> used for expression dataset aquisition",
    tag = "raw:geo",
    analysis = "GEO 数据获取",
    params = list(rna = FALSE)
    ))

.job_publish <- setClass("job_publish", 
  contains = c("job"),
  representation = representation(),
  prototype = prototype(
    method = "Supplementary file from article refer to"))

setClassUnion("job_PUBLISH", c("job_geo", "job_publish"))

job_geo <- function(id)
{
  .job_geo(object = id)
}

setMethod("step0", signature = c(x = "job_geo"),
  function(x){
    step_message("Prepare your data with function `job_geo`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_geo"),
  function(x, getGPL = TRUE){
    step_message("Get GEO metadata and information.")
    about <- e(GEOquery::getGEO(object(x), getGPL = getGPL))
    metas <- get_metadata.geo(about)
    prods <- get_prod.geo(metas)
    ## add GSE number
    prods <- .lich(c(list("Data Source ID" = object(x)), prods))
    prods <- .set_lab(prods, sig(x), object(x))
    x@params$about <- about
    x@params$metas <- metas
    x@params$prods <- prods
    guess <- metas$res[[1]]
    guess <- dplyr::rename_all(guess, make.names)
    guess <- dplyr::select(guess, 1:2, dplyr::ends_with(".ch1"))
    guess <- .set_lab(guess, sig(x), object(x), "metadata")
    x@params$guess <- guess
    x@params$test <- list(
      genes = try(as_tibble(about[[ 1 ]]@featureData@data), TRUE),
      counts = try(as_tibble(about[[ 1 ]]@assayData$exprs), TRUE),
      anno = about[[ 1 ]]@annotation,
      anno.db = try(.get_biocPackage.gpl(about[[ 1 ]]@annotation), TRUE)
    )
    x <- methodAdd(x, "以 R 包 `GEOquery` ({packageVersion('GEOquery')}) 获取 {object(x)} 数据集。")
    x <- snapAdd(x, "以 `GEOquery` 获取 {object(x)} 的数据信息。")
    return(x)
  })

.get_biocPackage.gpl <- function(gpl) {
  if (!is.character(gpl) | identical(gpl, character(0))) {
    stop("Incorrect character of gpl ID.")
  } else {
    data <- RCurl::getURL("https://gist.githubusercontent.com/seandavi/bc6b1b82dc65c47510c7/raw/b2027d7938896dce6145a1ebbcea75826813f6e1/platformMap.txt")
    data <- ftibble(data)
    dplyr::filter(data, gpl == !!gpl)$bioc_package
  }
}

setMethod("step2", signature = c(x = "job_geo"),
  function(x, filter_regex = NULL, baseDir = .prefix("GEO", "db"), rna = TRUE)
  {
    step_message("Download geo datasets or yellow{{RNA seq data}}.")
    if (rna) {
      if (dim(x$about[[1]])[1] > 0) {
        message("Is this a Microarray dataset?")
        return()
      }
      x$about[[1]] <- e(GEOquery::getRNASeqData(object(x)))
      message("Replace data in `x$about[[1]]`.")
      x <- methodAdd(x, "以 `GEOquery::getRNASeqData` 获取 RNA count 数据以及基因注释。")
      x$rna <- TRUE
    } else {
      if (!dir.exists(baseDir)) {
        dir.create(baseDir)
      }
      dir <- file.path(baseDir, object(x))
      continue <- 1L
      if (dir.exists(dir)) {
        continue <- usethis::ui_yeah(glue::glue("File exists ({dir}), continue?"))
      }
      if (continue) {
        e(GEOquery::getGEOSuppFiles(object(x), filter_regex = filter_regex,
            baseDir = baseDir))
        x$dir <- dir
      }
    }
    return(x)
  })

setMethod("meta", signature = c(x = "job_geo"),
  function(x, use = 1L){
    counts <- as_tibble(x@params$about[[ use ]]@assayData$exprs)
    genes <- as_tibble(x@params$about[[ use ]]@featureData@data)
    namel(counts, genes)
  })

setMethod("asjob_limma", signature = c(x = "job_geo"),
  function(x, metadata, use = 1L, normed = FALSE, use.col = NULL)
  {
    rna <- x$rna
    project <- object(x)
    if (rna) {
      counts <- as_tibble(data.frame(x@params$about[[ use ]]@assays@data$counts, check.names = FALSE))
      genes <- as_tibble(data.frame(x@params$about[[ use ]]@elementMetadata, check.names = FALSE))
      genes <- dplyr::mutate(genes, rownames = as.character(GeneID), .before = 1)
      if (missing(use.col)) {
        use.col <- "Symbol"
      } else {
        if (!is.character(genes[[ use.col ]])) {
          stop('`use.col` should character (symbol).')
        }
      }
    } else {
      counts <- as_tibble(x@params$about[[ use ]]@assayData$exprs)
      genes <- as_tibble(x@params$about[[ use ]]@featureData@data)
    }
    if (is.null(use.col)) {
      if (any(colnames(genes) == "rownames")) {
        if (any(grpl(colnames(genes), "Gene Symbol"))) {
          message("Rename 'Gene Symbol' as: hgnc_symbol")
          genes <- dplyr::relocate(genes, rownames, hgnc_symbol = `Gene Symbol`)
        }
      } else {
        guess <- grpf(colnames(genes), "^gene", ignore.case = TRUE)
        message("All available:\n\t", paste0(guess, collapse = ", "))
        message("Use the firist.")
        guess <- guess[[ 1 ]]
        keep <- !is.na(genes[[ guess ]]) & genes[[ guess ]] != "" & !duplicated(genes[[ guess ]])
        ## format
        genes <- dplyr::mutate(genes, rownames = !!rlang::sym(guess))
        genes <- dplyr::relocate(genes, rownames)
        message("Col.names of the data.frame:\n\t",
          stringr::str_trunc(paste0(colnames(counts), collapse = ", "), 40))
        counts <- dplyr::mutate(counts, rownames = !!genes[[ guess ]])
        counts <- dplyr::relocate(counts, rownames)
        genes <- genes[ keep, ]
        counts <- counts[ keep,  ]
      }
    } else {
      genes <- dplyr::relocate(genes, rownames, hgnc_symbol = !!rlang::sym(use.col))
    }
    if (normed) {
      counts <- dplyr::select(counts, -1)
      counts <- data.frame(counts)
      rownames(counts) <- genes$rownames
      x <- job_limma_normed(counts, metadata)
      x$genes <- genes
      x$rna <- rna
      x$project <- project
      x
    } else {
      cli::cli_alert_info("new_dge")
      res <- try(job_limma(new_dge(metadata, counts, genes)))
      if (inherits(res, "try-error")) {
        Terror <<- namel(metadata, counts, genes)
        stop("Error. Check `Terror`.")
      } else {
        res$rna <- rna
        res$project <- project
        return(res)
      }
    }
  })

setMethod("expect", signature = c(x = "job_geo", ref = "ANY"),
  function(x, ref, force = FALSE, id = x@object, ret = c("meta", "job")){
    ret <- match.arg(ret)
    if (missing(ref)) {
      ref <- geo_cols()
    }
    x$guess <- expect(x$guess, ref, force = force, id = id)
    if (ret == "job") {
      return(x)
    } else if (ret == "meta") {
      return(x$guess)
    }
  })

showStrings <- function(x, stat = TRUE) {
  if (any(duplicated(x))) {
    if (stat) {
      freq <- table(x)
      x <- paste0(names(freq), " (n=", unname(freq), ")")
    }
  }
  if (length(x) > 10) {
    x <- c(head(x, n = 10), "...")
  }
  stringr::str_wrap(bind("'", x, "'"), indent = 4, exdent = 4)
}

preset_group_string <- function(x) {
  knit_strings <- function(x) {
    x <- gs(x, "[()]", "")
    gs(x, "^[^a-zA-Z0-9]|[^a-zA-Z0-9]$", "")
  }
  if (!is.character(x)) {
    content <- showStrings(x)
    writeLines(paste0(crayon::yellow("Target not character:\n"), content))
    if (usethis::ui_yeah("Use that?")) {
      x <- as.character(x)
    } else {
      stop("...")
    }
  }
  if (isUnique <- all(table(x) == 1)) {
    message('all(table(x) == 1), try eliminating numbering.')
    elim <- function(ch) gs(ch, "[0-9]+[^0-9]*$", "")
    strings <- x
    while (isUnique && any(grpl(strings, "[0-9]"))) {
      strings <- knit_strings(elim(strings))
      isUnique <- all(table(strings) == 1)
    }
    if (min(nchar(strings)) && length(unique(strings)) != 1) {
      x <- strings
    } else {
      stop(
        'min(nchar(strings)) && length(unique(strings)) != 1, illegal final results:\n',
        showStrings(strings)
      )
    }
  }
  if (length(x) < 2) {
    stop('length(x) < 2, two few string to process.')
  }
  x <- strsplit(x, "\\s")
  len <- min(lengths(x))
  if (len > 2) {
    ## cut leading
    for (i in seq_len(len)) {
      chs <- vapply(x, function(ch) ch[i], character(1))
      if (length(unique(chs)) > 1) {
        break
      }
    }
    if (i > 1) {
      cut <- bind(x[[1]][seq_len(i - 1)], co = " ")
      message(glue::glue("Leading string cutoff: '{cut}'"))
      x <- lapply(x, function(x) x[ -seq_len(i - 1) ])
      len <- min(lengths(x))
    }
    ## cut trailing
    if (len > 2) {
      for (i in seq_len(len)) {
        chs <- vapply(x, function(ch) rev(ch)[i], character(1))
        if (length(unique(chs)) > 1) {
          break
        }
      }
      if (i > 1) {
        cut <- bind(rev(rev(x[[1]])[seq_len(i - 1)]), co = " ")
        message(glue::glue("Trailing string cutoff: '{cut}'"))
        x <- lapply(x, function(x) rev(rev(x)[ -seq_len(i - 1) ]))
      }
    }
  }
  x <- vapply(
    x, function(ch) {
      knit_strings(paste0(ch, collapse = " "))
    }, character(1)
  )
  x <- gs(x, " ", "_")
  message(crayon::yellow("[Function: preset_group_string] Final results:\n"), showStrings(x))
  return(x)
}

.pattern_geo_group <- function() {
  c("drug", "disease", "treatment", "protocol", "tissue", "title")
}

.expect_col_geo_group <- setClass("expect_col_geo_group",
  contains = c("expect_col"),
  prototype = prototype(
    name = "group", pattern_find = .pattern_geo_group(),
    fun_mutate = preset_group_string,
    fun_check = function(x) {
      isThat <- any(duplicated(x)) && length(unique(x)) > 1
      if (!isThat) {
        message("[Function: fun_check] At least one quantity must be greater than 1.")
      }
      isThat
    }
    ))

.expect_col_geo_sample <- setClass(".expect_col_geo_sample",
  contains = c("expect_col"),
  prototype = prototype(
    name = "sample", pattern_find = c("rownames"),
    fun_check = function(x) {
      isThat <- all(!duplicated(x))
      if (!isThat) {
        message("[Function: fun_check] All should be unique.")
      }
      isThat
    }
    ))

## do recode first, then mutate.
geo_cols <- function(
  group_mutate = preset_group_string, sample_mutate = NULL,
  group = .pattern_geo_group(), sample = "rownames",
  group_recode = character(0), sample_recode = character(0),
  db_file = .prefix("expect_cols_geo.rds", "db"),
  global = function(x) dplyr::relocate(x, sample, group),
  uniqueness = "rownames")
{
  cols <- list(
    .expect_col_geo_sample(
      pattern_find = sample, fun_mutate = sample_mutate, 
      pattern_recode = sample_recode
      ),
    .expect_col_geo_group(
      pattern_find = group, fun_mutate = group_mutate, 
      pattern_recode = group_recode
    )
  )
  .expect_cols_geo(
    cols, db_file = db_file, global = global, uniqueness = uniqueness
  )
}

.expect_cols_geo <- setClass("expect_cols_geo",
  contains = c("expect_cols"),
  prototype = prototype(
    list(.expect_col_geo_sample(), .expect_col_geo_group()),
    uniqueness = "rownames",
    global = function(x) {
      dplyr::relocate(x, sample, group)
    }))

setMethod("initialize", "expect_cols_geo", 
  function(.Object, ...) {
    .Object <- callNextMethod()
    if (!length(.Object@db_file)) {
      .Object@db_file <- .prefix("expect_cols_geo.rds", "db")
    }
    .Object
  })


get_metadata.geo <- function(lst,
  select = rlang::quos(rownames, title),
  pattern = c("diagnosis", "Sex", "^age", "^time point", "data_processing"),
  abbrev = c("data_processing"))
{
  res <- lapply(lst,
    function(eset){
      as_tibble(eset@phenoData@data)
    })
  if (!is.null(select)) {
    main <- lapply(res,
      function(data){
        cols <- colnames(data)
        extra <- unlist(.find_and_sort_strings(cols, pattern), use.names = FALSE)
        namesSelect <- vapply(select, rlang::as_label, character(1))
        select <- select[ which(namesSelect %in% cols) ]
        dplyr::select(data, !!!select, dplyr::all_of(extra))
      })
    if (!is.null(abbrev)) {
      abbrev <- lapply(res,
        function(data){
          cols <- colnames(data)
          extra <- unlist(.find_and_sort_strings(cols, abbrev), use.names = FALSE)
          dplyr::distinct(data, dplyr::pick(extra))
        })
    } else abbrev <- NULL
    res <- namel(main, abbrev, res)
  }
  lst_clear0(res)
}

get_prod.geo <- function(lst) {
  res <- as.list(dplyr::distinct(do.call(rbind, lst$abbrev)))
  .lich(res)
}
