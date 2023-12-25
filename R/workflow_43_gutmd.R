# ==========================================================================
# workflow of gutmd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gutmd <- setClass("job_gutmd", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: http://bio-annotation.cn/gutMDisorder"),
    cite = "[@GutmdisorderACheng2019]",
    method = "Database `gutMDisorder` used for finding associations between gut microbiota and metabolites"
    ))

job_gutmd <- function(db = "../Gut Microbe and Metabolite-human.txt")
{
  if (!file.exists(db)) {
    data <- e(RCurl::getURL("http://bio-annotation.cn/gutmgene/public/res/Gut%20Microbe%20and%20Metabolite-human.txt"))
    write_tsv(data, db)
    db_data <- data.table::fread(text = data)
  } else {
    db_data <- ftibble(db)
  }
  db_data <- dplyr::rename_all(db_data, make.names)
  .job_gutmd(object = db_data)
}

setMethod("step0", signature = c(x = "job_gutmd"),
  function(x){
    step_message("Prepare your data with function `job_gutmd`.")
  })

setMethod("step1", signature = c(x = "job_gutmd"),
  function(x, micro_patterns = NULL, metab_patterns = NULL, label = NULL){
    step_message("Match microbiota or metabolite in database")
    tables <- list()
    if (!is.null(micro_patterns)) {
      res <- list()
      res$matches <- matchThats(object(x)[[1]], micro_patterns)
      res$db_matched <- dplyr::filter(object(x),
        Gut.Microbiota %in% unlist(!!res$matches, use.names = F)
      )
      x$res_micro <- res
    }
    if (!is.null(metab_patterns)) {
      res <- list()
      alls <- rm.no(c(object(x)$Substrate, object(x)$Metabolite))
      res$matches <- matchThats(alls, metab_patterns)
      names <- unlist(res$matches, use.names = F)
      res$db_matched <- dplyr::filter(object(x),
        Substrate %in% !!names | Metabolite %in% !!names
      )
      x$res_metab <- res
      metab <- dplyr::select(res$db_matched, Metabolite, Substrate, Gut.Microbiota, Classification)
      metab <- .set_lab(metab, label, "gutMDisorder", "Matched metabolites and their related microbiota")
      tables$metab <- namel(metab)
    }
    x@tables <- tables
    return(x)
  })
