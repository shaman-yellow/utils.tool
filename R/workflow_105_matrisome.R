# ==========================================================================
# workflow of matrisome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_matrisome <- setClass("job_matrisome", 
  contains = c("job"),
  prototype = prototype(
    pg = "matrisome",
    info = c("https://doi.org/10.1074/mcp.M111.014647"),
    cite = "[@The_matrisome_Naba_2012]",
    method = "",
    tag = "matrisome",
    analysis = "Matrisome 基质体相关基因获取"
    ))

job_matrisome <- function()
{
  .job_matrisome()
}

setMethod("step0", signature = c(x = "job_matrisome"),
  function(x){
    step_message("Prepare your data with function `job_matrisome`.")
  })

setMethod("step1", signature = c(x = "job_matrisome"),
  function(x, type = c("human", "mouse")){
    step_message("Get data of matrisome.")
    type <- match.arg(type)
    cli::cli_alert_info("get_matrisome_data")
    db <- get_matrisome_data(type)
    x$db <- set_lab_legend(
      db, glue::glue("{x@sig} all {type} matrisome data"),
      glue::glue("为获取自文献[@The_matrisome_Naba_2012]的基质体数据。")
    )
    x$.feature <- split(
      s(x$db$Gene.symbol, "gene_", ""), x$db$Category
    )
    x <- snapAdd(x, "从 {cite_show('The_matrisome_Naba_2012')} (补充材料) 获取基质体 (matrisome) 数据。")
    return(x)
  })

get_matrisome_data <- function(type = c("human", "mouse"),
  file_db = .prefix("../data/matrisome_small.rds"))
{
  type <- match.arg(type)
  if (!file.exists(file_db)) {
    stop('!file.exists(file_db).')
  }
  db <- readRDS(file_db)
  db <- db[[type]]
  return(db)
}

