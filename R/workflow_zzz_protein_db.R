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

## package queryup
setGeneric("via_symbol", 
  function(x, ...) standardGeneric("via_symbol"))

#' @exportMethod via_symbol
setMethod("via_symbol", signature = c(x = "query_pdb"),
  function(x, symbol, max_pdb = 3, organism_id = "9606", cl = 1, run = "get_pdb.sh")
  {
    object(x) <- e(UniProt.ws::mapUniProt(
      from = "Gene_Name",
      to = "UniProtKB-Swiss-Prot",
      columns = c("accession", "id", "mass"),
      query = list(taxId = 9606, ids = 'pik3r1')
    ))
    lst <- touch_uniprotkb(object(x)$Entry)
    names(lst) <- make.names(object(x)$From)
    used <- lapply(lst, .get_needed)
    pdb_ids <- lapply(used, function(x) head(x$pdbs, max_pdb))
    pdb_files <- get_pdb(unlist(pdb_ids), cl = cl, mkdir.pdb = x$wd)
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
  x <- XML::xmlParse(lines, useInternalNodes = F)
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
    vapply(x, FUN.VALUE = character(1), USE.NAMES = F,
      function(x) {
        XML::xmlGetAttr(x, "id")
      })
  }
  list(active_sites = getThat(x.feature, "active site"),
    binding_sites = getThat(x.feature, "binding site"),
    pdbs = getId(getThat(x.dbReference, "PDB"))
  )
}

