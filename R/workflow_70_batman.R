# ==========================================================================
# workflow of batman
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_batman <- setClass("job_batman", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "batman",
    info = c("http://bionet.ncpsb.org.cn/batman-tcm/#/home"),
    cite = "[@BatmanTcm20Kong2024]",
    method = "Database `BATMAN-TCM` was used as source data of TCM ingredients and target proteins"
    ))

job_batman <- function()
{
  .job_batman()
}

setMethod("step0", signature = c(x = "job_batman"),
  function(x){
    step_message("Prepare your data with function `job_batman`.")
  })

setMethod("step1", signature = c(x = "job_batman"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })

get_batman_data <- function(savedir = .prefix("batman", "db")) {
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  rdata <- paste0(savedir, "/all.rds")
  if (!file.exists(rdata)) {
    url_base <- "http://batman2.cloudna.cn/downloadApiFile/data/browser/"
    files <- c(
      "herb_browse.txt",
      "known_browse_by_ingredients.txt.gz",
      "predicted_browse_by_ingredients.txt.gz"
    )
    lst <- pbapply::pblapply(files,
      function(file) {
        url <- paste0(url_base, file)
        target <- paste0(savedir, "/", file)
        if (!file.exists(target)) {
          x <- RCurl::getURLContent(url)
          writeBin(as.vector(x), target)
        }
        data <- ftibble(target, sep = "\t")
        if (ncol(data) == 1 & file == "predicted_browse_by_ingredients.txt.gz") {
          colnames(data) <- "Name"
          data <- dplyr::mutate(data,
            PubChem_CID = strx(Name, "^[0-9]+"),
            IUPAC_name = gs(Name, "^[^ ]* (.*) [^ ]*$", "\\1"),
            predicted_target_proteins = strx(Name, "[^ ]+$")
          )
          data <- dplyr::select(data, - Name)
        }
        data
      })
    names(lst) <- get_realname(files)
    return(lst)
  }
}
