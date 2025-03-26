# ==========================================================================
# workflow of htf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_htf <- setClass("job_htf", 
  contains = c("job"),
  prototype = prototype(
    pg = "htf",
    info = c("http://bioinfo.life.hust.edu.cn/hTFtarget#!/"),
    cite = "[@hTFtarget_A_Co_Zhang_2020]",
    method = "",
    tag = "tf",
    analysis = "hTFtarget 转录因子数据"
    ))

job_htf <- function(genes)
{
  x <- .job_htf(object = rm.no(genes))
  return(x)
}

setMethod("step0", signature = c(x = "job_htf"),
  function(x){
    step_message("Prepare your data with function `job_htf`.")
  })

setMethod("step1", signature = c(x = "job_htf"),
  function(x, dir_db = .prefix("hTFtarget", "db"), cl = NULL)
  {
    step_message("Get data...")
    anno <- tibble::tibble(symbol = object(x))
    anno <- map_gene(anno, "symbol", "SYMBOL", "ENSEMBL")
    ids <- anno[[ "ENSEMBL" ]]
    ids <- rm.no(ids)
    if (any(notThat <- !grpl(ids, "^ENSG"))) {
      message(glue::glue("Not ENSEMBL: {try_snap(notThat)}"))
      ids <- ids[ !notThat ]
    }
    # IDs: your query, col: the ID column, res: results table
    dir.create(dir_db, FALSE)
    db <- new_db(file.path(dir_db, "tf.rdata"), "target_ensembl")
    db <- not(db, ids)
    query <- db@query
    if (length(query)) {
      link <- start_drive()
      link$open()
      res <- pbapply::pbsapply(query, simplify = FALSE,
        function(id) {
          url <- glue::glue("http://bioinfo.life.hust.edu.cn/hTFtarget#!/targets/gene?gene={id}")
          link$navigate(url)
          message(glue::glue("Navigate: {url}"))
          Sys.sleep(1.5)
          content <- link$getPageSource()[[1]]
          html <- e(rvest::read_html(content))
          data <- e(rvest::html_table(html)[[1]])
        }, cl = cl)
      res <- rbind_list(res, "target_ensembl")
      link$close()
      end_drive()
      db <- upd(db, res)
    }
    res <- dplyr::filter(db@db, target_ensembl %in% !!ids)
    res <- map(
      res, "target_ensembl", anno, "ENSEMBL", "symbol", col = "target_symbol"
    )
    res <- dplyr::relocate(res, target_symbol, target_ensembl)
    res <- setLegend(res, "从 `hTFtarget` 数据库获取的 {less(object(x))} 结合的转录因子附表。")
    x <- tablesAdd(x, t.transcription_factor_data = res)
    ## feature
    res <- dplyr::distinct(res, target_symbol, TF)
    feature(x) <- split(res$TF, res$target_symbol)
    x <- methodAdd(x, "以 `RSelenium` 抓取 `hTFtarget` 数据库  (<http://bioinfo.life.hust.edu.cn/hTFtarget#!/>) {cite_show('hTFtarget_A_Co_Zhang_2020')} 基因与转录因子结合数据。")
    x <- snapAdd(x, "从 `hTFtarget` 数据库获取与 {less(object(x))} 结合的转录因子数据。")
    return(x)
  })

