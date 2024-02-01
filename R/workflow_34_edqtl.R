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
    info = c("https://gtexportal.org/home/downloads/adult-gtex#qtl",
      "https://www.nature.com/scitable/topicpage/quantitative-trait-locus-qtl-analysis-53904/"),
    cite = "[@TheGtexConsorNone2020]",
    method = "The QTL data were abtained from GTEx database"
    ))

job_edqtl <- function(mode = c("edqtl", "eqtl"))
{
  mode <- match.arg(mode)
  lst <- list(
    edqtl = list(
      path = .prefix("gtex/edqtl", "db"),
      file = "bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.tar",
      patterns = c("edsite", "signif_variant_site_pairs"),
      anno_file = "bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.README.txt"),
    eqtl = list(
      path = .prefix("gtex/eqtl", "db"),
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
  x$mode <- gs(mode, "qtl", "QTL")
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
    step_message("Read the tissue QRL files.")
    db <- sapply(x$patterns,
      function(pattern) {
        file <- list.files(x$db_path, paste0(tissue, ".*", pattern, ".*"), full.names = T, recursive = T)
        message("Search get files:\n\t", paste0(file, collapse = ", "))
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
    use.eq <- .set_lab(use.eq, sig(x), paste("all used", x$mode, "data"))
    x$use.eq <- use.eq
    return(x)
  })

setMethod("step3", signature = c(x = "job_edqtl"),
  function(x, filter, label = "DEGs", use.filter = "hgnc_symbol", data = x@params$use.eq){
    step_message("Filter the edqtl data.")
    lst <- c(nl(label, list(filter)), list(DB_edQTL = data[[ use.filter ]]))
    p.venn <- new_venn(lst = lst)
    p.venn <- .set_lab(p.venn, sig(x), paste("database of", x$mode, "intersect with"), label)
    data <- dplyr::filter(data, !!rlang::sym(use.filter) %in% !!filter)
    data <- .set_lab(data, sig(x), paste0("database of ", x$mode, " intersect with ", label, " DATA"))
    x@plots[[ 3 ]] <- namel(p.venn)
    x@tables[[ 3 ]] <- list(filtered.eq = data)
    return(x)
  })

setMethod("map", signature = c(x = "job_edqtl", ref = "job_publish"),
  function(x, ref, label = "Filtered_edQTL", data = x@tables$step3$filtered.eq){
    if (ref@cite == "[@MendelianRandoLiuX2022]") {
      data <- dplyr::select(data, variant_id, hgnc_symbol)
      filter_fun <- function(refd, ref_label) {
        lst <- c(nl(label, list(data$variant_id)), nl(ref_label, list(refd$variant_id)))
        p.venn <- new_venn(lst = lst)
        p.venn <- .set_lab(p.venn, sig(x), paste0("filtered ", x$mode, " data intersect with"), ref_label)
        refd <- dplyr::select(refd, 1:2)
        data <- tbmerge(refd, data, by = "variant_id")
        data <- .set_lab(data, sig(x), paste0("filtered ", x$mode, " data intersect with ", ref_label), "DATA")
        namel(data, p.venn)
      }
      lst.micro <- filter_fun(object(ref)$snp_microbiota, "microbiota related")
      attr(lst.micro$data, "pattern") <- ref$get_pattern(lst.micro$data$Microbiome.features)
      lst.metab <- filter_fun(object(ref)$snp_metabolite, "metabolite related")
      object <- namel(lst.micro, lst.metab)
    }
    x <- .job_publish(object = object, cite = ref@cite)
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


