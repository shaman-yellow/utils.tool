# ==========================================================================
# workflow of genecard
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_genecard <- setClass("job_genecard", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.genecards.org/"),
    cite = "[@TheGenecardsSStelze2016]",
    method = "The Human Gene Database `GeneCards` used for disease related genes prediction"
    ))

job_genecard <- function(disease)
{
  .job_genecard(object = disease)
}

setMethod("step0", signature = c(x = "job_genecard"),
  function(x){
    step_message("Prepare your data with function `job_genecard`.")
  })

setMethod("step1", signature = c(x = "job_genecard"),
  function(x, score = 5){
    step_message("Get from GeneCards website.")
    t.genecards <- get_from_genecards(object(x), score = score)
    t.genecards <- .set_lab(t.genecards, sig(x), "disease related targets from GeneCards")
    x@tables[[ 1 ]] <- namel(t.genecards)
    return(x)
  })

get_from_genecards <- function(query, score = 5, keep_drive = F) {
  link <- start_drive(browser = "firefox")
  Sys.sleep(3)
  link$open()
  link$navigate(url)
  query <- gs(query, " ", "%20")
  url <- paste0('https://www.genecards.org/Search/Keyword?queryString=', query,
    '&pageSize=25000&startPage=0')
  html <- link$getPageSource()[[1]]
  link$close()
  end_drive()
  html <- XML::htmlParse(html)
  table <- XML::readHTMLTable(html)
  table <- as_tibble(data.frame(table[[1]]))
  colnames(table) %<>% gs("\\.+", "_")
  colnames(table) %<>% gs("X_|_$", "")
  table <- select(table, -1, -2)
  table <- filter(table, Score > !!score)
  return(table)
}
