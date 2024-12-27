# ==========================================================================
# workflow of m6a
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_m6a <- setClass("job_m6a", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://rnamd.org/m6a/api.php"),
    cite = "[@M6aAtlasV20Liang2024]",
    method = "The API of `m6A-Atlas` used for obtaining m6A related data from the website",
    tag = "m6a:site",
    analysis = "m6A-Atlas m6A 数据获取"
    ))

job_m6a <- function(species = c("HomoSapiens", "MusMusculus"))
{
  .job_m6a(object = match.arg(species))
}

setMethod("step0", signature = c(x = "job_m6a"),
  function(x){
    step_message("Prepare your data with function `job_m6a`.")
  })

setMethod("step1", signature = c(x = "job_m6a"),
  function(x, genes, save = .prefix("m6a_Atlas", "db"), sleep = .1)
  {
    step_message("Get data from the website.")
    dir.create(save, FALSE)
    types <- c("LowResolution", "HighResolution")
    t.data <- sapply(types, simplify = FALSE,
      function(type) {
        db <- new_db(paste0(save, "/", object(x), "_", type, ".rds"), ".id")
        db <- not(db, genes)
        if (length(db@query)) {
          cli::cli_alert_info(paste0("Obtain data of ", type))
          res <- pbapply::pbsapply(db@query, simplify = FALSE,
            function(gene) {
              url <- paste0("http://www.rnamd.org/m6a/api", type, ".php?species=", object(x), "&gene=", gene)
              content <- RCurl::getURL(url)
              Sys.sleep(sleep)
              if (identical(content, "")) {
                NULL
              } else {
                data <- try(data.table::fread(content))
                if (inherits(data, "try-error")) {
                  Terror <<- data
                  stop("Can not resolve the unexpected response.")
                } else {
                  as_tibble(data)
                }
              }
            })
          res <- frbind(res, idcol = TRUE, fill = TRUE)
          db <- upd(db, res, db@query)
        }
        res(db, what = genes)
      })
    t.data <- .set_lab(t.data, sig(x), names(t.data), "m6A-Atlas search results")
    lab(t.data) <- glue::glue("{sig(x)} m6A-Atlas search results")
    x@tables[[ 1 ]] <- namel(t.data)
    x <- methodAdd(x, "以 `m6A-Atlas` 数据库 {cite_show('M6aAtlasV20Liang2024')} 提供的 API 获取所需基因的 m6A 修饰靶点 (<http://rnamd.org/m6a/api.php>)。")
    return(x)
  })


