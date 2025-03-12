# ==========================================================================
# workflow of gtopdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gtopdb <- setClass("job_gtopdb", 
  contains = c("job"),
  prototype = prototype(
    pg = "gtopdb",
    info = c("https://www.guidetopharmacology.org/"),
    cite = "[@The_IUPHAR_BPS_Hardin_2024]",
    method = "",
    tag = "gtopdb",
    analysis = "GtoPdb 药理学靶点及实验配体数据"
    ))

job_gtopdb <- function(dir_save = .prefix("gtopdb", "db"))
{
  dir.create(dir_save, FALSE)
  url <- "https://www.guidetopharmacology.org/DATA/targets_and_families.csv"
  file <- file.path(dir_save, basename(url))
  if (!file.exists(file)) {
    download.file(url, file)
  }
  db <- ftibble(file)
  .job_gtopdb(object = db)
}

setMethod("step0", signature = c(x = "job_gtopdb"),
  function(x){
    step_message("Prepare your data with function `job_gtopdb`.")
  })

setMethod("step1", signature = c(x = "job_gtopdb"),
  function(x){
    step_message("Filter ...")
    return(x)
  })
