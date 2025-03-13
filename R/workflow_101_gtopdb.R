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
  colnames(db) <- gs(colnames(db), " ", "_")
  db <- set_lab_legend(db, "GtoPdb targets and families", "为 GtoPdb 数据库所有 targets and families 条目。")
  x <- .job_gtopdb(object = db)
  x <- methodAdd(x, "获取 GtoPdb 数据库 'target and family' {cite_show('The_IUPHAR_BPS_Hardin_2024')} (<https://www.guidetopharmacology.org/download.jsp>)。")
  return(x)
}

setMethod("step0", signature = c(x = "job_gtopdb"),
  function(x){
    step_message("Prepare your data with function `job_gtopdb`.")
  })

setMethod("step1", signature = c(x = "job_gtopdb"),
  function(x, mode = c("ic")){
    step_message("Filter ...")
    mode <- match.arg(mode)
    if (mode == "ic") {
      types <- c("vgic", "lgic", "other_ic")
      x$db <- dplyr::filter(object(x), Type %in% !!types)
      x$db <- set_lab_legend(
        x$db, "GtoPdb family of ion channels", "为 GtoPdb 数据库所有 ion channels family。"
      )
      p.pie <- wrap(new_pie(x$db$Family_name), 8, 7)
      p.pie <- setLegend(p.pie, "为 GtoPdb 所有 Ion channels 分布饼图。")
      x <- plotsAdd(x, p.distribution_of_ion_channels = p.pie)
      feature(x) <- rm.no(x$db$HGNC_symbol)
      x <- snapAdd(x, "获取 `GtoPdb` 数据库所有 Ion channels 类型靶点。")
    }
    return(x)
  })
