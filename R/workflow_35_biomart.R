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
    method = "R package `biomaRt` used for gene annotation",
    tag = "gene:anno",
    analysis = "Biomart 基因注释"
    ))

job_biomart <- function(mart_dataset, global = T, clear = F)
{
  x <- .job_biomart()
  if (clear) {
    options(biomart = NULL)
  }
  x$mart <- new_biomart(mart_dataset)
  if (global) {
    options(mart_dataset = mart_dataset)
    message("Set '", mart_dataset, "' for global `new_biomart` use.")
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
    filters = c("ensembl_transcript_id", "ensembl_gene_id",
      "hgnc_symbol", "mgi_symbol", "entrezgene_id", "rgd_symbol"),
    attrs = NULL, mode = x$mart_dataset)
  {
    step_message("Get annotation.")
    filters <- match.arg(filters)
    if (is.null(attrs)) {
      if (mode == "hsa") {
        attrs <- general_attrs()
      } else if (mode == "mmu") {
        attrs <- c("mgi_symbol", "ensembl_transcript_id", "ensembl_gene_id", "entrezgene_id", "description")
      } else if (mode == "rno") {
        attrs <- c("rgd_symbol", "ensembl_transcript_id", "ensembl_gene_id", "entrezgene_id", "description")
      }
    }
    attrs <- unique(c(attrs, filters))
    anno <- filter_biomart(x$mart, attrs, filters, values)
    anno <- .set_lab(anno, sig(x), "genes", "annotation")
    x$anno <- anno
    return(x)
  })

setMethod("map", signature = c(x = "job_biomart", ref = "job_edqtl"),
  function(x, ref, use = "signif_variant_site_pairs", label.factor = 1){
    ## prepare for editing site
    db <- object(ref)[[ use ]]
    ed.site <- unique(db$gene_id)
    group_chr <- function(x) {
      chr <- gs(x, "^chr([0-9]+).*", "\\1")
      split(x, chr)
    }
    ed.site <- group_chr(ed.site)
    ed.site <- lapply(ed.site,
      function(x) {
        site <- gs(x, ".*_([0-9]+)$", "\\1")
        site <- as.integer(site)
        names(site) <- x
        site
      })
    ## prepare for main data (annotation)
    genes <- dplyr::filter(x$anno, grpl(chromosome_name, "^[0-9]+$"))
    genes <- split(genes, ~ chromosome_name)
    ## match
    res <- sapply(names(ed.site), simplify = F,
      function(chr) {
        site <- ed.site[[ chr ]]
        genes <- genes[[ chr ]]
        res <- lapply(site,
          function(si) {
            isThat <- genes$start_position <= si & si <= genes$end_position
            genes[ isThat, ]$hgnc_symbol
          })
        lst_clear0(res)
      })
    res <- lst_clear0(res)
    names(res) <- NULL
    res <- unlist(res)
    ## get the db data
    db <- dplyr::filter(db, gene_id %in% dplyr::all_of(names(res)))
    db <- dplyr::mutate(db, ref_gene = res[match(gene_id, names(res))])
    db <- dplyr::relocate(db, ref_gene, gene_id)
    data <- dplyr::select(db, Editing_Site = gene_id, Variant = variant_id, Symbol = ref_gene)
    data <- dplyr::mutate(data, Editing_Site = gs(Editing_Site, "_", "\n"))
    p.edqtl <- new_allu(data, col.fill = 3, axes = 1:3, label.auto = T, label.factor = label.factor)
    p.edqtl <- .set_lab(p.edqtl, sig(x), "The matched RNA editing site")
    x$p.edqtl <- p.edqtl
    db <- .set_lab(db, sig(x), "The matched RNA editing site DATA")
    x$edqtl <- db
    return(x)
  })
