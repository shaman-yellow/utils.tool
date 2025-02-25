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
  .job_geo(object = strx(id, "GSE[0-9]+"))
}

setMethod("step0", signature = c(x = "job_geo"),
  function(x){
    step_message("Prepare your data with function `job_geo`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_geo"),
  function(x, getGPL = TRUE, dir_cache = .prefix("GEO", "db"))
  {
    step_message("Get GEO metadata and information.")
    if (!is.null(dir_cache)) {
      dir.create(dir_cache, FALSE)
      file_cache <- file.path(
        dir_cache, paste0(object(x), "_", getGPL, ".rds")
      )
      if (file.exists(file_cache)) {
        message(glue::glue('file.exists(file_cache): {file_cache}'))
        about <- readRDS(file_cache)
      } else {
        message(glue::glue("Download {object(x)}..."))
        about <- e(GEOquery::getGEO(object(x), getGPL = getGPL))
        saveRDS(about, file_cache)
      }
    } else {
      about <- e(GEOquery::getGEO(object(x), getGPL = getGPL))
    }
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
      if (TRUE) {
        if (TRUE) {
          quantifications <- getRNASeqQuantResults_custom(object(x))
        } else {
          quantifications <- e(GEOquery:::getRNASeqQuantResults(object(x)))
        }
        se <- e(SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = quantifications$quants),
            rowData = quantifications$annotation
            ))
        supp <- ""
        if (any(isThat <- !x$guess$rownames %in% colnames(se))) {
          x$ncbiNotGot <- notGot <- x$guess$rownames[isThat]
          message(
            glue::glue('any(!x$guess$rownames %in% colnames(se)): Not got: {bind(notGot)}')
          )
          supp <- glue::glue("缺失样本: {bind(notGot)} ('NCBI-generated data' 缺失样本计数数据的原因包括运行未通过 50% 的对齐率或由于技术原因处理失败)")
        }
        x <- methodAdd(
          x, "以 `GEOquery:::getRNASeqQuantResults` 获取 RNA count 数据
          (NCBI-generated data, 参考 <https://www.ncbi.nlm.nih.gov/geo/info/rnaseqcounts.html>)
          {supp} 以及基因注释。"
        )
        x$about[[1]] <- se
      } else {
        ## as this some times error, due to sample missing.
        x$about[[1]] <- e(GEOquery::getRNASeqData(object(x)))
        message("Replace data in `x$about[[1]]`.")
        x <- methodAdd(x, "以 `GEOquery::getRNASeqData` 获取 RNA count 数据 (NCBI-generated data) 以及基因注释。")
      }
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
        files <- list.files(dir, "\\.tar", full.names = TRUE)
        if (length(files)) {
          lapply(files, 
            function(file) {
              utils::untar(file, exdir = normalizePath(dir))
            })
        }
        files <- list.files(dir, ".", full.names = TRUE)
        x$dir_files <- files
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

setMethod("map", signature = c(x = "job_limma", ref = "job_geo"),
  function(x, ref, which = 1L, type = "pipe", ...)
  {
    message("Get extra annotation (Symbol) if possible.")
    if (!any(colnames(object(x)$genes) == "ID")) {
      stop('!any(colnames(object(x)$genes) == "ID").')
    }
    gpl <- ref@params$about[[which]]@annotation
    file_anno <- paste0(gpl, "_", type, ".rda")
    if (!file.exists(file_anno)) {
      anno <- e(AnnoProbe::idmap(gpl, type = type))
    } else {
      name <- load(file_anno)
      anno <- get(name)
    }
    object(x)$genes <- map(
      object(x)$genes, "ID", anno, "probe_id", "symbol", col = "symbol"
    )
    return(x)
  })

setMethod("asjob_limma", signature = c(x = "job_geo"),
  function(x, metadata, use = 1L, normed = "guess", use.col = NULL)
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
        if (any(isThat <- grpl(colnames(genes), "Gene Symbol"))) {
          oname <- colnames(genes)[ isThat ][1]
          message(glue::glue("Rename '{oname}' as: hgnc_symbol:\n{showStrings(genes[[oname]])}"))
          genes <- dplyr::relocate(genes, rownames, hgnc_symbol = !!rlang::sym(oname))
        }
      } else {
        guess <- grpf(colnames(genes), "^gene", ignore.case = TRUE)
        if (length(guess)) {
          message("All available:\n\t", paste0(guess, collapse = ", "))
          message("Use the firist.")
          guess <- guess[[ 1 ]]
        } else {
          which <- menuThat(colnames(genes), "Use which as gene ID?")
          guess <- colnames(genes)[which]
        }
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
    if (any(colnames(genes) == "gene_assignment")) {
      genes <- dplyr::mutate(
        genes, GENE_SYMBOL = strx(
          gene_assignment, "(?<=// )[^/ ]+(?= //)"
        ), .before = 2
      )
    } else if (any(colnames(genes) == "SPOT_ID.1")) {
      genes <- dplyr::mutate(
        genes, GENE_SYMBOL = gs(
          SPOT_ID.1, "[^/]+// RefSeq // [^(]+\\(([^)]+)\\).*", "\\1"
        ),
        .before = 2
      )
    }
    message(
      glue::glue("Gene annotation:\n{showStrings(colnames(genes), trunc = FALSE)}")
    )
    if (!rna && identical(normed, "guess")) {
      message("Guess whether data has been normed.")
      if (all(range(counts[, -1], na.rm = TRUE) >= 0)) {
        message('all(range(counts[, -1]) > 0), the data has not been normed.')
        normed <- FALSE
      } else {
        message('all(range(counts[, -1]) > 0) == FALSE, the data has been normed.')
        normed <- TRUE
      }
    } else {
      normed <- FALSE
    }
    if (!is.null(x$ncbiNotGot)) {
      message(glue::glue("Filter out NCBI generated data (Missing samples): {bind(x$ncbiNotGot)}"))
      metadata <- dplyr::filter(metadata, !rownames %in% x$ncbiNotGot)
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

showStrings <- function(x, stat = TRUE, trunc = TRUE) {
  if (any(duplicated(x))) {
    if (stat) {
      freq <- table(x)
      x <- paste0(names(freq), " (n=", unname(freq), ")")
    }
  }
  if (length(x) > 10 && trunc) {
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
  x <- gs(x, " |-", "_")
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
  res <- as.list(dplyr::distinct(do.call(dplyr::bind_rows, lst$abbrev)))
  if (any(lengths(res) > 1)) {
    res <- lapply(res, bind, co = " ### ")
  }
  .lich(res)
}

getRNASeqQuantResults_custom <- function(gse) {
  fun_download <- function(gse) {
    if (length(gse) != 1 || !grpl(gse, "^GSE")) {
      stop('length(gse) != 1 || !grpl(gse, "^GSE").')
    }
    strs <- e(RCurl::getURL(glue::glue("https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse}")))
    xmls <- XML::htmlParse(strs)
    fun_get <- function(x, xpath) {
      res <- XML::xpathApply(x, xpath)
      paste0("https://www.ncbi.nlm.nih.gov", XML::xmlGetAttr(res[[1]], "href"))
    }
    url_data <- fun_get(xmls, "//a[@id='download_raw_counts']")
    url_anno <- fun_get(xmls, "//h2[text()='Human gene annotation table']/..//a")
    quants <- e(GEOquery:::readRNAQuantRawCounts(url_data))
    annotation <- e(GEOquery:::readRNAQuantAnnotation(url_anno))
    list(quants = quants, annotation = annotation)
  }
  dir.create(dir <- .prefix("GEO", "db"), FALSE)
  expect_local_data(
    dir, paste0("rnaQuant_", gse), fun_download, list(gse = gse)
  )
}
