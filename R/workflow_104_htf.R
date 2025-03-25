# ==========================================================================
# workflow of htf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_htf <- setClass("job_htf", 
  contains = c("job"),
  prototype = prototype(
    pg = "htf",
    info = c(""),
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
  function(x){
    step_message("Get data...")

    # link <- start_drive()
    # link$open()
    # link$navigate("http://bioinfo.life.hust.edu.cn/hTFtarget#!/targets/gene?gene=ENSG00000263761")
    # content <- link$getPageSource()
    # end_drive()
    #
    # html <- rvest::read_html(content[[1]])
    # tbl_htf <- rvest::html_table(html)[[1]]
    # tbl_htf

    return(x)
  })

