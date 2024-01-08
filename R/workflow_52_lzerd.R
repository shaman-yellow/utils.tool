# ==========================================================================
# workflow of lzerd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_lzerd <- setClass("job_lzerd", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://lzerd.kiharalab.org/"),
    cite = "[@LzerdWebserverChrist2021]",
    method = "`LZerD` webserver used for protein–protein docking"
    ))

job_lzerd <- function(symbols)
{
  .job_lzerd(object = symbols)
}

setMethod("step0", signature = c(x = "job_lzerd"),
  function(x){
    step_message("Prepare your data with function `job_lzerd`.")
  })

setMethod("step1", signature = c(x = "job_lzerd"),
  function(x, unique = T){
    step_message("Download pdb files for protein docking.")
    mart <- new_biomart()
    x$anno <- filter_biomart(mart, c("hgnc_symbol", "pdb"), "hgnc_symbol",
      object(x), distinct = F)
    x$anno <- dplyr::filter(x$anno, pdb != "")
    data <- dplyr::distinct(x$anno, hgnc_symbol, .keep_all = T)
    x$layouts <- combn(data$pdb, 2, simplify = F)
    x$pdb_files <- get_pdb(unique(unlist(x$layouts)))
    return(x)
  })
