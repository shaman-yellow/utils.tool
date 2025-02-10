# ==========================================================================
# workflow of lnctard
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_lnctard <- setClass("job_lnctard", 
  contains = c("job"),
  prototype = prototype(
    pg = "lnctard",
    info = c("https://lnctard.bio-database.com/Download"),
    cite = "[@LncTarD_2_0_an_Zhao_2023]",
    method = "",
    tag = "lnctard",
    analysis = "LncTarD LncRNA调控靶点"
    ))

job_lnctard <- function(dir_db = .prefix("lnctard", "db"))
{
  dir.create(dir_db, FALSE)
  file_url <- "https://lnctard.bio-database.com/downloadfile/lnctard2.0.zip"
  # --no-check-certificate
  file <- file.path(dir_db, basename(file_url))
  file_data <- paste0(tools::file_path_sans_ext(file), ".txt")
  if (!file.exists(file_data) && !file.exists(file)) {
    download.file(file_url, file)
    utils::unzip(file, exdir = dir_db)
  }
  db <- ftibble(file_data)
  x <- .job_lnctard()
  x <- methodAdd(x, "从 `LncTarD` 数据库{cite_show('LncTarD_2_0_an_Zhao_2023')}获取 LncRNA 调控数据 (<https://lnctard.bio-database.com/downloadfile/lnctard2.0.zip>)。")
  x$db <- db
  return(x)
}

setMethod("step0", signature = c(x = "job_lnctard"),
  function(x){
    step_message("Prepare your data with function `job_lnctard`.")
  })

setMethod("step1", signature = c(x = "job_lnctard"),
  function(x, lnc){
    step_message("Find target genes.")
    t.regulate <- dplyr::filter(x@params$db, Regulator %in% !!lnc)
    t.regulate <- setLegend(t.regulate, "为 LncRNA ({less(lnc)}) 调控的靶点基因附表。")
    x <- tablesAdd(x, t.lncRNA_regulation_data = t.regulate)
    x <- snapAdd(x, "从 `LncTarD` 数据库检索 {less(lnc)} 的调控靶基因。")
    return(x)
  })
