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

new_query_pdb <- function(wd = ".") {
  x <- .query_pdb()
  x$mart <- new_biomart()
  x$wd <- wd
  return(x)
}

## package queryup
setGeneric("via_symbol", 
  function(x, ...) standardGeneric("via_symbol"))

setMethod("via_symbol", signature = c(x = "query_pdb"),
  function(x, symbol){
  })

# query <- list("gene_exact" = "Pik3r1", "reviewed" = "true", "organism_id" = "9606")
# df <- query_uniprot(query, columns = c("id", "mass"), show_progress = FALSE)



