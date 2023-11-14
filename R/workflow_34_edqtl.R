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

job_edqtl <- function(path = "../gtex/edqtl", file = "bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.tar")
{
  .check_untar(paste0(path, "/", file))
  x <- .job_edqtl()
  files <- list.files(path, "txt.gz$", full.names = T)
  x$metadata <- tibble::tibble(
    tissue = stringr::str_extract(get_filename(files), "^[^.]*"),
    files = files
  )
  x$db_path <- path
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
    db <- sapply(c("edsite", "signif_variant_site_pairs"),
      function(name) {
        file <- list.files(x$db_path, paste0(tissue, ".*", name, ".*"), full.names = T)
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

.check_untar <- function(file, path = get_path(file)) {
  recordfile <- paste0(path, "/.check_untar")
  if (!file.exists(recordfile)) {
    untar(file, exdir = path)
    writeLines("", recordfile)
  }
}

setMethod("anno", signature = c(x = "job_edqtl"),
  function(x, file = "../gtex/edqtl/bulk-qtl_v8_editing-qtl_GTEx_Analysis_v8_edQTL.README.txt"){
    system(paste0("vim ", file))
  })


