# ==========================================================================
# workflow of gtex
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_edqtl <- setClass("job_edqtl", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://gtexportal.org/home/downloads/adult-gtex#qtl"),
    cite = "[@TheGtexConsorNone2020]",
    method = "The edQTL data were abtained from GTEx database"
    ))

job_edqtl <- function(mode = c("edqtl", "eqtl"))
{
  mode <- match.arg(mode)
  lst <- list(
    edqtl = list(
      path = "../gtex/edqtl",
      file = "bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.tar",
      patterns = c("edsite", "signif_variant_site_pairs"),
      anno_file = "bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.README.txt"),
    eqtl = list(
      path = "../gtex/eqtl",
      file = "bulk-qtl_v8_single-tissue-cis-qtl_GTEx_Analysis_v8_eQTL.tar",
      patterns = c("egenes", "signif_variant_gene_pairs"),
      anno_file = "bulk-qtl_v8_single-tissue-cis-qtl_README_eQTL_v8.txt"
    )
  )
  lst <- lst[[ mode ]]
  .check_untar(paste0(lst$path, "/", lst$file))
  x <- .job_edqtl()
  files <- list.files(lst$path, "txt.gz$", full.names = T, recursive = T)
  x$metadata <- tibble::tibble(
    tissue = stringr::str_extract(get_filename(files), "^[^.]*"),
    files = files
  )
  x$db_path <- lst$path
  x$anno_file <- paste0(lst$path, "/", lst$anno_file)
  x$patterns <- lst$patterns
  return(x)
}

setMethod("step0", signature = c(x = "job_edqtl"),
  function(x){
    step_message("Prepare your data with function `job_edqtl`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_edqtl"),
  function(x, tissue){
    step_message("Read the tissue edQRL files.")
    db <- sapply(x$patterns,
      function(pattern) {
        file <- list.files(x$db_path, paste0(tissue, ".*", pattern, ".*"), full.names = T, recursive = T)
        if (length(file) == 0) {
          stop("Check the input tissue as no file retrieved.")
        } else if (length(file) > 1) {
          stop("Too many files found.")
        } else {
          ftibble(file)
        }
      })
    object(x) <- db
    return(x)
  })

setMethod("step2", signature = c(x = "job_edqtl"),
  function(x, job_biomart = NULL, use = "signif_variant_gene_pairs"){
    step_message("Get hgnc_symbol for `gene_id`.")
    if (is.null(job_biomart)) {
      bt <- job_biomart("hsa")
      bt <- step1(bt, unique(gname(object(x)[[ use ]]$gene_id)), "ensembl_gene_id")
      x$job_biomart <- bt
    }
    use.eq <- object(x)[[ use ]]
    use.eq <- dplyr::mutate(use.eq, symbol = gname(gene_id))
    use.eq <- map(use.eq, "symbol", bt$anno, "ensembl_gene_id", "hgnc_symbol")
    use.eq <- dplyr::filter(use.eq, !is.na(hgnc_symbol))
    x$use.eq <- use.eq
    return(x)
  })

setMethod("map", signature = c(x = "job_edqtl", ref = "job_publish"),
  function(x, ref, filter, use.filter = "hgnc_symbol"){
    if (ref@cite == "[@MendelianRandoLiuX2022]") {
    }
    x <- .job_publish()
    return(x)
  })

.check_untar <- function(file, path = get_path(file)) {
  recordfile <- paste0(path, "/.check_untar")
  if (!file.exists(recordfile)) {
    untar(file, exdir = path)
    writeLines("", recordfile)
  }
}

setMethod("anno", signature = c(x = "job_edqtl"),
  function(x, file = x$anno_file){
    system(paste0("vim ", file))
  })


