# ==========================================================================
# get data (pdb) from protein database
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## hgnc_symbol
## pdb ids
## pdb files
## visualization (use pymol)

## Active site (queryup)
## molecular weight
## other infomation

#' @import methods

.query_pdb <- setClass("query_pdb", 
  contains = c("job"),
  representation = representation(),
  prototype = NULL)

#' @export new_pdb
#' @export new_query_pdb
new_pdb <- new_query_pdb <- function(wd = "protein_pdb") {
  x <- .query_pdb()
  # x$mart <- new_biomart()
  if (!dir.exists(wd)) {
    dir.create(wd)
  }
  x$wd <- wd
  return(x)
}

#' @exportMethod via_symbol
setMethod("via_symbol", signature = c(x = "query_pdb"),
  function(x, symbol, max_pdb = 3, organism_id = "9606", cl = 1, fun_download = download_pdb.RCurl)
  {
    object(x) <- e(UniProt.ws::mapUniProt(
      from = "Gene_Name",
      to = "UniProtKB-Swiss-Prot",
      columns = c("accession", "id", "mass"),
      query = list(taxId = 9606, ids = symbol)
    ))
    lst <- touch_uniprotkb(object(x)$Entry)
    names(lst) <- make.names(object(x)$From)
    used <- lapply(lst, .get_needed)
    x$used <- used
    pdb_ids <- lapply(used, function(x) head(x$pdbs, max_pdb))
    pdb_files <- get_pdb(unlist(pdb_ids), cl = cl, mkdir.pdb = x$wd,
      fun_download = download_pdb.RCurl)
    names(pdb_files) <- toupper(names(pdb_files))
    x$pdb_files <- pdb_files
    x$pdb_ids <- pdb_ids
    return(x)
  })

#' @exportMethod show
setMethod("show", signature = c(object = "query_pdb"),
  function(object){
    if (!is.null(object$pdb_ids)) {
      lst <- lapply(object$pdb_ids, paste0, collapse = " ")
      show(.lich(lst))
    } else {
      callNextMethod(object)
    }
  })

#' @exportMethod vis
setMethod("vis", signature = c(x = "query_pdb"),
  function(x, symbol, id = NULL){
    id <- match.arg(id, x$pdb_ids[[ symbol ]])
    data <- readLines(x$pdb_files[[ id ]])
    require(r3dmol)
    obj <- r3dmol::r3dmol(viewer_spec = m_viewer_spec(cartoonQuality = 10,
        lowerZoomLimit = 50, upperZoomLimit = 350))
    obj <- m_add_model(obj, data = data, format = "pdb")
    obj <- m_zoom_to(obj)
    obj <- m_set_style(obj, style = m_style_cartoon(color = '#00cc96'))
    obj <- m_set_style(obj, sel = m_sel(ss = 's'),
      style = m_style_cartoon(color = '#636efa', arrows = FALSE))
    obj <- m_set_style(obj, sel = m_sel(ss = 'h'),
      style = m_style_cartoon(color = '#ff7f0e'))
    obj <- m_rotate(obj, angle = 90, axis = 'y')
    obj <- m_spin(obj)
    attr(obj, ".PDB_id") <- id
    obj
  })

#' @exportMethod anno
setMethod("anno", signature = c(x = "query_pdb"),
  function(x, symbols, ids = NULL, title = "Proteins Infomation",
    .name = "annotation", .style = "BiocStyle::html_document")
  {
    objs <- lapply(seq_along(symbols),
      function(n) {
        vis(x, symbols[[n]], ids[[n]])
      })
    ids <- vapply(objs, function(x) attr(x, ".PDB_id"), character(1))
    names(objs) <- ids
    meta <- apply(object(x), 1, c, simplify = FALSE)
    names(meta) <- symbols
    names(ids) <- symbols
    des_all <- lapply(symbols,
      function(symbol) {
        info_protein <- baseInfo_protein(meta[[ symbol ]])
        des_protein <- descrip_protein(x$used[[ symbol ]])
        structure_protein <- structure_protein(ids[[ symbol ]])
        c(paste0("# ", symbol), "", info_protein, des_protein, structure_protein, "")
      })
    lines <- unlist(des_all)
    header <- c("---",
      paste0("title: ", title),
      paste0("output: ", .style), "---", "")
    writeLines(c(header, lines), md <- paste0(x$wd, "/", .name, ".Rmd"))
    file <- rmarkdown::render(md)
    utils::browseURL(file)
  })

structure_protein <- function(ids) {
  lines <- unlist(lapply(ids, function(id) {
      c(paste0("PDB ID: ", ids), "",
        "```{r, echo = FALSE}",
        paste0("objs[[\"", id, "\"]]"),
        "```", "")
      }))
  c("## Structure View", "", lines)
}

baseInfo_protein <- function(lst) {
  c("## Information", "",
    paste0("- Symbol: ", lst[[ "From" ]]),
    paste0("- UniProtKB Swiss Prot: ", lst[[ "Entry" ]]),
    paste0("- Mass: ", lst[[ "Mass" ]]),
    ""
  )
}

descrip_protein <- function(lst.xml) {
  fun_site <- function(lst.xml, name, fun) {
    xmls <- lst.xml[[ name ]]
    if (length(xmls) >= 1) {
      unlist(lapply(xmls, fun), use.names = FALSE)
    } else "No records"
  }
  active <- fun_site(lst.xml, "active_sites", text_activeSite)
  binding <- fun_site(lst.xml, "binding_sites", text_bindingSite)
  c("## Active sites", "", active,
    "", "",
    "## Binding sites", "", binding,
    ""
  )
}

text_bindingSite <- function(xml) {
  rExtractChildren(xmlToList2(xml), fix = mdlist,
    function(x) x$attrs[["type"]],
    function(x) x$name,
    function(x) {
      if (x$name == "dbReference") {
        paste0(x$name, ": ", x$attrs[[ "id" ]])
      } else {
        if (length(x$attrs) > 0) {
          paste0(x$name, " ", names(x$attrs), ": ", unname(x$attrs))
        } else {
          paste0(x$name, ":")
        }
      }
    },
    function(x) {
      paste0(x$value)
    }
  )
}

text_activeSite <- function(xml) {
  rExtractChildren(xmlToList2(xml), fix = mdlist,
    function(x) x$attrs[["description"]],
    function(x) x$name,
    function(x) {
      paste0(names(x$attrs), ": ", unname(x$attrs))
    }
  )
}

xmlToList2 <- function(xml) {
  name <- XML::xmlName(xml)
  attrs <- XML::xmlAttrs(xml)
  value <- XML::xmlValue(xml)
  if (identical(value, character(0))) {
    value <- NULL
  }
  if (length(xml) > 0L) {
    children <- XML::xmlSApply(xml, xmlToList2, simplify = FALSE)
  } else {
    children <- NULL
  }
  list(name = name, attrs = attrs, value = value, children = children)
}

mdlist <- function(str, n) {
  indent <- paste0(rep(" ", (n - 1) * 4), collapse = "")
  paste0(indent, "- ", str)
}

rExtractChildren <- function(lst, ..., fix = NULL) {
  args <- list(...)
  dep <- length(args)
  n <- 0L
  hasSet <- 0L
  lines <- c()
  rfun <- function(lst) {
    if (!hasSet) {
      n <<- n + 1L
      hasSet <<- 1L
    }
    fun_get <- args[[ n ]]
    str <- fun_get(lst)
    if (is.function(fix)) {
      str <- fix(str, n)
    }
    lines <<- c(lines, str)
    if (n <= dep) {
      if (length(lst$children) >= 1) {
        hasSet <<- 0L
        n.before <- n
        lapply(lst$children, rfun)
        n <<- n.before
      }
    }
  }
  rfun(lst)
  lines
}

#' @export touch_uniprotkb
touch_uniprotkb <- function(ids) {
  pbapply::pblapply(ids,
    function(id) {
      url <- paste0("https://rest.uniprot.org/uniprotkb/", id, ".xml")
      res <- RCurl::getURL(url)
      res
    })
}

.get_needed <- function(lines) {
  x <- XML::xmlParse(lines, useInternalNodes = FALSE)
  x <- x$doc$children$uniprot$children$entry$children
  x.feature <- x[ names(x) == "feature" ]
  x.dbReference <- x[ names(x) == "dbReference" ]
  getThat <- function(x, types) {
    isThat <- vapply(x, FUN.VALUE = logical(1),
      function(x) {
        XML::xmlGetAttr(x, "type") %in% types
      })
    x[ isThat ]
  }
  getId <- function(x) {
    vapply(x, FUN.VALUE = character(1), USE.NAMES = FALSE,
      function(x) {
        XML::xmlGetAttr(x, "id")
      })
  }
  list(active_sites = getThat(x.feature, "active site"),
    binding_sites = getThat(x.feature, "binding site"),
    pdbs = getId(getThat(x.dbReference, "PDB"))
  )
}

get_pdb <- function(ids, cl = 3, mkdir.pdb = "protein_pdb", fun_download = download_pdb.RCurl) {
  cli::cli_alert_info("Obtain PDB files")
  dir.create(mkdir.pdb, FALSE)
  ids <- tolower(ids)
  if (length(ids) == 0) {
    stop("length(ids) == 0")
  }
  exists <- gs(list.files(mkdir.pdb, ".*\\.pdb$"), "\\.pdb$", "")
  ids <- ids[ !ids %in% exists ]
  if (length(ids) > 0) {
    fun_download(ids, cl = cl, mkdir.pdb = mkdir.pdb)
  }
  files <- list.files(mkdir.pdb, ".*\\.pdb$", full.names = TRUE)
  names(files) <- basename(tools::file_path_sans_ext(files))
  files
}

download_pdb.RCurl <- function(ids, cl, mkdir.pdb) {
  pbapply::pblapply(ids, cl = cl,
    function(id) {
      url <- paste0("https://files.rcsb.org/download/", id, ".pdb")
      utils::download.file(url, file.path(mkdir.pdb, paste0(id, ".pdb")))
    })
}

download_pdb.shell <- function(ids, cl, mkdir.pdb) {
  ids <- grouping_vec2list(ids, round(length(ids) / cl), TRUE)
  pbapply::pblapply(ids, cl = cl,
    function(ids) {
      tmp <- tempfile(fileext = ".txt")
      cat(ids, sep = ", ", file = tmp)
      system(paste0("get_pdb.sh", " -f ", tmp, " -p -o ", mkdir.pdb))
    })
  lapply(list.files(mkdir.pdb, ".*\\.gz$", full.names = TRUE), R.utils::gunzip)
}


