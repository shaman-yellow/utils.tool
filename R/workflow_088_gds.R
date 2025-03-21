# ==========================================================================
# workflow of gds
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gds <- setClass("job_gds", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "gds",
    info = c("https://www.ncbi.nlm.nih.gov/books/NBK179288/"),
    cite = "",
    method = "",
    tag = "gds",
    analysis = "GSE 数据搜索"
    ))

# setMethod("snap", signature = c(x = "job_gds"),
  # function(x, patterns)
  # {
  #   if (is.null(names(patterns))) {
  #     stop("names should not be NULL")
  #   }
  #   texts <- lapply(names(patterns), function(name) {
  #     data <- dplyr::filter(object(x),
  #       grpl(summary, patterns[[name]], TRUE) | grpl(title, patterns[[name]], TRUE)
  #     )
  #     if (nrow(data)) {
  #       alls <- paste0(data$GSE, collapse = ", ")
  #       glue::glue("与 {name} 相关的数据共 {nrow(data)} 个，分别为: {alls}。")
  #     } else {
  #       glue::glue("未找到与 {name} 相关的数据。")
  #     }
  #   })
  #   paste0(paste0("- ", texts), collapse = "\n")
  # })

setMethod("focus", signature = c(x = "job_gds"),
  function(x, patterns)
  {
    if (is.null(names(patterns))) {
      stop("names should not be NULL")
    }
    lapply(names(patterns), function(name) {
      dplyr::filter(object(x),
        grpl(summary, patterns[[name]], TRUE) | grpl(title, patterns[[name]], TRUE)
      )
    })
  })

job_gds <- function(keys,
  n = "6:300",
  # Mus musculus, Rattus norvegicus
  org = "Homo sapiens",
  elements = c("GSE", "title", "summary", "taxon", "gdsType",
    "n_samples", "PubMedIds", "BioProject"),
  ...)
{
  query <- add_field(keys, "[Description]")
  if (is.integer(n)) {
    n <- paste0(min(n), ":", max(n))
  }
  extra <- glue::glue(
    "({n}[Number of Samples]) AND (GSE[Entry Type])"
  )
  if (!is.null(org)) {
    extra <- paste(extra, glue::glue("AND ({org}[Organism])"))
  }
  object <- try(
    edirect_db("gds", c(query, extra), elements, ...), TRUE
  )
  if (!inherits(object, "try-error")) {
    object <- dplyr::mutate(object, GSE = paste0("GSE", GSE))
    object <- .set_lab(object, keys[1], "EDirect query")
  } else {
    object <- NULL
    snapAdd_onExit("x", "未发现合适数据集。")
  }
  x <- .job_gds(object = object, params = list(query = query, elements = elements, org = org, n = n))
  x <- methodAdd(x, "使用 Entrez Direct (EDirect) <https://www.ncbi.nlm.nih.gov/books/NBK3837/> 搜索 GEO 数据库 (`esearch -db gds`)，查询信息为: ({paste0(paste0('(', query, ''), collapse = ' AND ')}) AND ({extra})，转化为数据表格。")
  x <- snapAdd(x, "以 Entrez Direct (EDirect) 搜索 GEO 数据库 (检索条件见方法章节) 。")
  x
}

setMethod("step1", signature = c(x = "job_gds"),
  function(x, ..., single_cell = FALSE, clinical = TRUE, 
    rna_or_array = TRUE, mRNA = TRUE,
    excludes = paste0(
      "\\bKO\\b|\\bWT\\b|_WT_|_KO_|wildtype|mutant|knock|deficien|absen",
      "|SuperSeries|transgenic|CD[0-9]+"
      ),
    extras = NULL, deep_inspect = TRUE,
    control = "normal|control|healthy|ctrl|adjacent|\\bN[0-9]+\\b|\\bC[0-9]+\\b",
    expects = NULL, cl = 5
  )
  {
    step_message("Statistic.")
    args <- rlang::enquos(...)
    fun_expects <- function() {
      if (!is.null(expects)) {
        isThat <- expects %in% object(x)$GSE
        glue::glue("Got expects: {try_snap(isThat)} {if (any(!isThat)) bind(expects[ !isThat ]) else ''}")
      } else {
        ""
      }
    }
    fun_elements <- function(x) {
      gs(strsplit(x, "\\|")[[1]], "\\\\b", "")
    }
    if (!is.null(single_cell)) {
      if (single_cell) {
        object(x) <- dplyr::filter(
          object(x), grpl(paste(title, summary), "single[ -]cell|scRNA", TRUE)
        )
      } else {
        object(x) <- dplyr::filter(
          object(x), !grpl(paste(title, summary), "single[ -]cell|scRNA", TRUE)
        )
      }
      snap <- if (single_cell) "仅保留" else "滤除"
      x <- methodAdd(x, "以正则匹配，{snap} 'summary' 或 'title' 中包含 'single cell' 或 'scRNA' 的数据例。")
      message(glue::glue("dim: {bind(dim(object(x)))}, single_cell == {single_cell}. {fun_expects()}"))
    }
    if (!is.null(clinical)) {
      if (clinical) {
        keyClinical <- "in vitro|cell line|CD[0-9]+|vehicle|vector|DMSO|/ml|\\bnm\\b"
        object(x) <- dplyr::filter(
          object(x), !grpl(summary, keyClinical, TRUE)
        )
        x <- methodAdd(x,
          "仅查询临床数据，因此滤除匹配到关键词 {bind(fun_elements(keyClinical))} 的数据例。
          (注：以上仅为查找合适的 GEO 数据所做的数据筛选，与实际分析无关) 。"
        )
      }
      message(glue::glue("dim: {bind(dim(object(x)))}, clinical == {clinical}. {fun_expects()}"))
    }
    if (!is.null(rna_or_array)) {
      if (rna_or_array) {
        object(x) <- dplyr::filter(object(x), 
          grpl(gdsType, "Expression profiling by high throughput sequencing|Expression profiling by array"))
        x <- methodAdd(x, "仅获取类型包含 'Expression profiling by high throughput sequencing' 或 'Expression profiling by array' 的数据例。")
      }
      message(glue::glue("dim: {bind(dim(object(x)))}, rna_or_array == {rna_or_array}. {fun_expects()}"))
    }
    if (!is.null(excludes)) {
      if (!is.null(extras)) {
        excludes <- paste0(excludes, "|", extras)
      }
      object(x) <- dplyr::filter(
        object(x),
        !grpl(summary, excludes, TRUE) & !grpl(title, excludes, TRUE)
      )
      x$excludes <- excludes
      EX <- fun_elements(excludes)
      x <- methodAdd(x, "此外，排除 summary 或 title 中匹配到字符集 EX ({bind(EX)}) 的数据例。")
      message(glue::glue("dim: {bind(dim(object(x)))}, excludes .... {fun_expects()}"))
    }
    if (deep_inspect) {
      if (nrow(object(x)) > 200) {
        stop(glue::glue('nrow(object(x)) > 200, too many data to fetch ({nrow(object(x))}).'))
      } else if (!nrow(object(x))) {
        stop('nrow(object(x)) == FALSE')
      }
      others <- expect_geo_extra(object(x)$GSE, cl = cl)
      object(x) <- tbmerge(object(x), others, by = "GSE", all.x = TRUE)
      if (!is.null(excludes)) {
        x <- methodAdd(x, "上述得到 {nrow(object(x))} 个 GSE 数据集。以 R 抓取网页 (例如，`RCurl::getURL` 抓取 <{object(x)$GSE[1]}>)，解析 'Overall design' 和 'Samples' 模块，匹配字符集 EX，排除匹配到的数据例。")
        object(x) <- dplyr::filter(
          object(x), !grpl(design, excludes, TRUE) & !grpl(samples, excludes, TRUE)
        )
        message(glue::glue("dim: {bind(dim(object(x)))}, deep_inspect, excludes .... {fun_expects()}"))
      }
      if (!is.null(clinical) && clinical) {
        object(x) <- dplyr::filter(
          object(x), 
          !grpl(paste0(design, samples), keyClinical, TRUE)
        )
        x <- methodAdd(x, "排除 'Overall design' 中包含 {bind(fun_elements(keyClinical))} 的数据例。")
        message(glue::glue("dim: {bind(dim(object(x)))}, deep_inspect, clinical .... {fun_expects()}"))
      }
      if (!is.null(control)) {
        message(crayon::red("Match control sample by pattern."))
        object(x) <- dplyr::filter(
          object(x), 
          grpl(paste(samples, summary, design), !!control, TRUE) |
            nchar(paste(summary, design)) < 300
        )
        x <- methodAdd(
          x, "为了获得包含对照样本的数据集，(描述字符过少的例外) 以正则匹配 ('Samples', 'Summary', 'Overall design') 包含 {bind(fun_elements(control))} 的数据例。共 {nrow(object(x))} 个。"
        )
        message(glue::glue("dim: {bind(dim(object(x)))}, deep_inspect, match control .... {fun_expects()}"))
      }
      if (mRNA) {
        otherRNA <- "siRNA|miRNA|\\bmiR\\b|lncRNA"
        object(x) <- dplyr::filter(
          object(x), !grpl(paste(design, samples), otherRNA, TRUE)
        )
        x <- methodAdd(x, "仅获取包含 'protein coding' 测序的数据集，排除 'Samples' 和 'Overall design' 中包含 {bind(fun_elements(otherRNA))} 字符的数据例。余下共 {nrow(object(x))} 个。")
        message(glue::glue("dim: {bind(dim(object(x)))}, deep_inspect, mRNA .... {fun_expects()}"))
      }
    }
    if (length(args)) {
      object(x) <- trace_filter(object(x), quosures = args)
      x <- snapAdd(x, "{snap(object(x))}")
      gses <- bind(head(object(x)$GSE, n = 30))
      x <- snapAdd(
        x, "这些数据为：{less(gses, 10)}等。"
      )
      message(glue::glue("dim: {bind(dim(object(x)))}, args .... {fun_expects()}"))
    }
    p.AllGdsType <- wrap(new_pie(object(x)$gdsType), 7, 7)
    x <- plotsAdd(x, p.AllGdsType)
    return(x)
  })

setMethod("step2", signature = c(x = "job_gds"),
  function(x, which = "all", force = FALSE, ...)
  {
    step_message("Try get metadata of GSE dataset.")
    if (identical(which, "all")) {
      gses <- object(x)$GSE
    } else {
      gses <- object(x)$GSE[ which ]
    }
    if (length(gses) > 50 && !force) {
      stop('length(gse) > 50 && !force, too many items for query.')
    }
    x$res <- batch_geo(gses, ...)
    names(x$res$metas) <- gses
    if (TRUE) {
      x$res$res <- NULL
    } else {
      names(x$res$res) <- gses
    }
    x$querys <- gses
    snap <- length(which(vapply(x$res$metas, is, logical(1), "df")))
    x <- methodAdd(x, "以 `GEOquery` 获取 GSE 数据集 (n={snap})。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_gds"),
  function(x, ref, mode = c("survival", "control"), from_backup = TRUE, not = FALSE)
  {
    step_message("Search in metadata.")
    if (missing(ref)) {
      mode <- match.arg(mode)
      if (mode == "survival") {
        ref <- "Survival|Event|Dead|Alive|Status|Day|Time"
      } else if (mode == "control") {
        ref <- "normal|control|healthy|ctrl|adjacent"
      }
    }
    if (from_backup && !is.null(x$backup)) {
      object(x) <- x$backup$object
      x$res$metas <- x$backup$metas
    }
    which <- which(hunt(x, ref, "which", not = not))
    x <- methodAdd(x, "从元数据中匹配{if (not) '并排除' else ''}包含关键词的数据：'{ref}'，共得到 {length(which)} 个数据集。")
    x$backup <- namel(object = object(x), metas = x$res$metas)
    object(x) <- object(x)[which, ]
    x$res$metas <- x$res$metas[ which ]
    return(x)
  })

setMethod("expect", signature = c(x = "job_gds", ref = "ANY"),
  function(x, ref, force = FALSE, ids = names(x$res$metas),
    sleep = NULL, return_type = c("self", "data.frame"), ...)
  {
    if (x@step < 2L) {
      stop('x@step < 2L, should have metadata of all dataset.')
    }
    if (missing(ref)) {
      ref <- geo_cols()
    }
    return_type <- match.arg(return_type)
    metas <- x$res$metas
    if (length(metas) != length(ids)) {
      stop('length(metas) != length(ids), not match?')
    }
    if (length(force) == 1 && is.null(names(force))) {
      force <- rep(force, length(ids))
      names(force) <- ids
    } else if (length(force) < length(ids)) {
      if (is.null(names(force))) {
        stop('is.null(names(force)), do not know how to match force.')
      }
      if (any(!names(force) %in% ids)) {
        stop('any(!names(force) %in% ids), some not match?')
      }
      extraForce <- ids[ !ids %in% names(force) ]
      extraForce <- nl(extraForce, rep(FALSE, length(extraForce)), FALSE)
      force <- c(force, extraForce)
    }
    metas <- mapply(metas, ids, SIMPLIFY = FALSE,
      FUN = function(meta, id) {
        expect(meta, ref, force = force[[id]], id = id, sleep = sleep, ...)
      })
    if (return_type == "self") {
      x$res$metas <- metas
      return(x)
    } else {
      return(metas)
    }
  })

setMethod("anno", signature = c(x = "job_gds"),
  function(x, col = "group", group_limit = 8) {
    if (x@step < 2L) {
      stop('x@step < 2L, should have metadata of all dataset.')
    }
    metas <- x$res$metas
    if (!identical(names(metas), x@object$GSE)) {
      stop('!identical(names(metas), x@object$GSE), can not match GSE id in object(x)')
    }
    types <- x@object$gdsType
    summaries <- x@object$summary
    snaps <- mapply(
      metas, names(metas), types, summaries, SIMPLIFY = FALSE,
      FUN = function(data, id, type, summary) {
        scPattern <- "single|scRNA"
        if (grpl(summary, scPattern) || any(grpl(unlist(data), scPattern))) {
          extype <- " (scRNA-seq) "
        } else {
          extype <- ""
        }
        string <- data[[ col ]]
        if (is.null(string)) {
          stop(glue::glue('is.null(string), can not found "{col}" in metadata of "{id}".'))
        }
        x <- table(string)
        if (length(x) < group_limit) {
          snap_freqs <- glue::glue("{names(x)} (n = {unname(x)})")
          snap_freqs <- stringr::str_wrap(
            paste0("- ", snap_freqs), 100, 4, 4
          )
          snap_freqs <- bind(snap_freqs, co = "\n")
        } else {
          snap_freqs <- "..."
        }
        type <- .mutate_gdsType(type)
        glue::glue("- **{id}**, **Type**: {extype}{type}\n{snap_freqs}")
      }
    )
    snaps <- bind(unlist(snaps), co = "\n")
    x <- snapAdd(x, paste0("可用数据，及其组别为：\n\n", snaps), step = "a", add = FALSE)
    return(x)
  })

.mutate_gdsType <- function(x) {
  x <- gs(x, "profiling by array", "Microarray")
  x <- gs(x, "profiling by high throughput sequencing", "RNA-seq")
  x <- gs(x, "Expression", "")
  x <- gs(x, "coding RNA", "coding")
  x
}

add_field <- function(keys, field) {
  if (any(which <- !grpl(keys, "\\["))) {
    keys[which] <- paste0(keys[which], field)
  }
  keys
}

edirect_db <- function(db, query, elements,
  cache = .prefix(paste0("query_", db), "db"),
  force = FALSE)
{
  dir.create(cache, FALSE)
  if (length(query)) {
    query <- paste0("(", query, ")")
    query <- paste0(query, collapse = " AND ")
  }
  target <- file.path(cache, make.names(query))
  if (!file.exists(target) || force) {
    script <- glue::glue("esearch -db {db} -query '{query}' |\\
      efetch -format docsum |\\
      xtract -pattern DocumentSummary \\
      -sep '|' -element {paste0(elements, collapse = ' ')} \\
      > {target}"
    )
    cdRun(script)  
  }
  object <- as_tibble(read.csv(target, header = FALSE, sep = "\t"))
  colnames(object) <- elements
  object
}

setMethod("active", signature = c(x = "job_gds"),
  function(x, ids = NULL, which = NULL, n = 10)
  {
    if (is.null(ids)) {
      if (!is.null(which)) {
        ids <- object(x)$GSE[ which ]
      } else {
        ids <- head(object(x)$GSE, n = n)
      }
    }
    urls <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ids)
    lapply(urls, browseURL)
  })

setMethod("hunt", signature = c(x = "job_gds", ref = "character"),
  function(x, ref, get = c("meta", "gse", "summary", "which"), not = FALSE){
    if (is.null(x$res$metas)) {
      stop('is.null(x$res$metas), has not run "step2" ?')
    }
    get <- match.arg(get)
    isThat <- vapply(x$res$metas, 
      function(data) {
        if (is(data, "df")) {
          any(vapply(data, 
            function(col) {
              any(grpl(col, ref, TRUE))
            }, logical(1)) | grpl(colnames(data), ref, TRUE))
        } else FALSE
      }, logical(1))
    if (not) {
      isThat <- !isThat
    }
    if (get == "which") {
      return(isThat)
    } else if (get == "meta") {
      return(x$res$metas[ isThat ])
    } else {
      if (length(x$res$metas) != nrow(object(x))) {
        stop('length(x$res$metas) != object(x), you have filter the "object(x)" manually?')
      }
      if (get == "summary") {
        return(object(x)[isThat, ])
      } else if (get == "gse") {
        return(object(x)$GSE[isThat])
      }
    }
  })

setMethod("vis", signature = c(x = "job_gds"),
  function(x, ref = NULL, n = 10){
    data <- object(x)
    if (!is.null(ref)) {
      data <- dplyr::filter(data, grpl(summary, !!ref))
    }
    print(data, n = n)
  })

setMethod("step0", signature = c(x = "job_gds"),
  function(x){
    step_message("Prepare your data with function `job_gds`.")
  })

batch_geo <- function(gses, getGPL = FALSE, cl = 1L, 
  db_dir = .prefix("gds_geoMeta", "db"))
{
  if (!dir.exists(db_dir)) {
    dir.create(db_dir)
  }
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop('!requireNamespace("GEOquery", quietly = TRUE), no package of "GEOquery".')
  }
  usedPackages <- c(
    "GEOquery", "dplyr", "RCurl", "usethis", "stringr",
    "cli", "crayon", "glue", "rlang"
  )
  lapply(usedPackages, requireNamespace, quietly = TRUE)
  res <- pbapply::pblapply(gses,
    function(x) {
      file_cache <- file.path(db_dir, paste0(x, ".rds"))
      if (file.exists(file_cache)) {
        return(readRDS(file_cache))
      }
      times <- 0L
      res <- NULL
      while((is.null(res) || inherits(res, "try-error")) && times <= 5L) {
        times <- times + 1L
        res <- try(
          suppressMessages(step1(job_geo(x), getGPL = getGPL, dir_cache = NULL)), FALSE
        )
        if (times > 1) {
          cli::cli_alert_info("Got error, try again in 3 second.")
          Sys.sleep(3)
        }
      }
      if (is.null(res)) {
        cli::cli_alert_danger(x)
      } else {
        cli::cli_alert_success(x)
      }
      data <- .job_geo(object = x,
        step = 1L, params = list(guess = try(res$guess, TRUE))
      )
      saveRDS(data, file_cache)
      return(res)
    }, cl = cl
  )
  metas <- lapply(res, function(x) try(x$guess))
  namel(res, metas)
}

expect_geo_extra <- function(gses, cl = NULL, nGroup = 20,
  dir_db = .prefix("GEOquery", "db"))
{
  if (any(!grpl(gses, "^GSE[0-9]+"))) {
    stop('any(!grpl(gses, "^GSE[0-9]+")).')
  }
  # IDs: your query, col: the ID column, res: results table
  dir.create(dir_db, F)
  db <- new_db(file.path(dir_db, "GEO_extras.rdata"), "GSE")
  db <- not(db, gses)
  query <- db@query
  if (length(query)) {
    if (length(query) > nGroup) {
      query <- grouping_vec2list(query, nGroup, TRUE)
      n <- 0L
      lapply(query, 
        function(query) {
          n <<- n + 1L
          message(glue::glue("Group: {n}"))
          res <- .get_geo_extras(query, cl = cl)
          db <<- upd(db, res, ids = query)
        })
    } else {
      res <- .get_geo_extras(query, cl = cl)
      db <- upd(db, res)
    }
  }
  res <- dplyr::filter(db@db, GSE %in% !!gses)
}

.get_geo_extras <- function(gses, cl = NULL, maxTry = 10) {
  urls <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gses)
  htmls <- pbapply::pblapply(urls, cl = cl, 
    function(url) {
      res <- NULL
      n <- 0L
      while (is.null(res) ||
        (n <= maxTry && inherits(res, "try-error") &&
          grpl(res, "0A000126:SSL routines")))
      {
        if (!is.null(res)) {
          message(glue::glue("\nTry again ({n}): {url}\n"))
        }
        n <- n + 1L
        res <- try(RCurl::getURL(url), TRUE)
      }
      if (inherits(res, "try-error")) {
        print(res)
        stop('inherits(res, "try-error"), see above. (', n, ')')
      }
      res
    })
  lst <- lapply(htmls, 
    function(html) {
      xml <- XML::htmlParse(html)
      node_design <- XML::xpathApply(xml, "//table//td[text()='Overall design']/../td[@style]")
      if (is.null(node_design)) {
        design <- ""
      } else if (length(node_design) != 1) {
        stop('length(node_design) != 1.')
      } else {
        design <- XML::xmlValue(node_design)
      }
      node_samples <- XML::xpathApply(
        xml, "//tr[@valign='top']/td[matches(text(), '^Samples')]/..//tr"
      )
      samples <- XML::xmlValue(node_samples)
      samples <- samples[ grpl(samples, "GSM") ]
      samples <- lapply(
        strsplit(samples, "\n"), setNames, nm = c("sample", "title")
      )
      samples <- dplyr::bind_rows(samples)
      namel(design, samples = paste0(samples$title, collapse = " // "))
    })
  Terror <<- data <- dplyr::bind_rows(lst)
  dplyr::mutate(data, GSE = !!gses, .before = 1)
}

