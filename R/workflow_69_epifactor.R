# ==========================================================================
# workflow of epifactor
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_epifactor <- setClass("job_epifactor", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "epifactor",
    info = c("https://epifactors.autosome.org/description"),
    cite = "[@Epifactors2022Maraku2023]",
    method = "Database `EpiFactors` used for screening epigenetic regulators",
    tag = "epi",
    analysis = "EpiFactors 表观遗传调控因子数据获取"
    ))

job_epifactor <- function(use = c("protein"), version = "v2.0")
{
  x <- .job_epifactor()
  x$use <- match.arg(use)
  x$version <- version
  return(x)
}

setMethod("step0", signature = c(x = "job_epifactor"),
  function(x){
    step_message("Prepare your data with function `job_epifactor`.")
  })

setMethod("step1", signature = c(x = "job_epifactor"),
  function(x, dir = .prefix("epifactor", "db")){
    step_message("Obtain epigenetic regulators.")
    file <- switch(x$use, protein = "EpiGenes_main.csv")
    dir.create(dir, FALSE)
    db_file <- file.path(dir, paste0(file, "_", x$version, ".rds"))
    if (file.exists(db_file)) {
      data <- readRDS(db_file)
    } else {
      url <- paste0("https://epifactors.autosome.org/public_data/", x$version, "/", file)
      data <- ftibble(RCurl::getURL(url))
      saveRDS(data, db_file)
    }
    data <- .set_lab(data, sig(x), paste0("all ", x$use), "of epigenetic regulators")
    x@tables[[ 1 ]] <- nl(x$use, list(data))
    x <- methodAdd(x, "从数据库 `EpiFactors` {cite_show('Epifactors2022Maraku2023')} 获取表观遗传调控蛋白的数据。")
    x <- snapAdd(x, "获取 `EpiFactors` 中的表观遗传调控蛋白。")
    x$.feature <- as_feature(data$HGNC_symbol, x)
    return(x)
  })

setMethod("filter", signature = c(x = "job_epifactor"),
  function(x, types = c("RNA methylation"), ...)
  {
    types <- paste0(types, collapse = "|")
    x@tables$step1$protein <- trace_filter(
      x@tables$step1$protein, grpl(Modification, !!types), ...
    )
    x$.feature <- as_feature(x@tables$step1$protein$HGNC_symbol, x)
    x <- snapAdd(x, snap(x@tables$step1$protein), add = TRUE)
    return(x)
  })


