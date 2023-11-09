# ==========================================================================
# workflow of biomart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_biomart <- setClass("job_biomart", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.bioconductor.org/packages/release/bioc/html/biomaRt.html"),
    cite = "[@MappingIdentifDurinc2009]",
    method = "Package biomaRt used for gene annotation"
    ))

job_biomart <- function(mart_dataset, global = T)
{
  x <- .job_biomart()
  x$mart <- new_biomart(mart_dataset)
  if (global) {
    options(mart_dataset = "mmu")
    message("Set 'mmu' for global `new_biomart` use.")
  }
  x$mart_dataset <- mart_dataset
  return(x)
}

setMethod("step0", signature = c(x = "job_biomart"),
  function(x){
    step_message("Prepare your data with function `job_biomart`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_biomart"),
  function(x, values,
    filters = c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
    attrs = NULL, mode = x$mart_dataset)
  {
    step_message("Get annotation.")
    filters <- match.arg(filters)
    if (is.null(attrs)) {
      if (mode == "hsa") {
        attrs <- general_attrs()
      } else if (mode == "mmu") {
        attrs <- c("mgi_symbol", "ensembl_transcript_id", "entrezgene_id", "hgnc_symbol", "description")
      }
    }
    x$anno <- filter_biomart(x$mart, attrs, filters, values)
    return(x)
  })
