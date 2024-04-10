# ==========================================================================
# workflow of uniprotkb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_uniprotkb <- setClass("job_uniprotkb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.uniprot.org/help/api_queries"),
    cite = "",
    method = paste0("The API of `UniProtKB` (<https://www.uniprot.org/help/api_queries>) ",
      "used for mapping of names or IDs of proteins")
    ))

job_uniprotkb <- function(query, db_file = "uniprotkb/query.rds", organism_id = 9606)
{
  x <- .job_uniprotkb(object = rm.no(query))
  x$organism_id <- organism_id
  x$db_file <- db_file
  return(x)
}

setMethod("step0", signature = c(x = "job_uniprotkb"),
  function(x){
    step_message("Prepare your data with function `job_uniprotkb`.")
  })

setMethod("step1", signature = c(x = "job_uniprotkb"),
  function(x, cl = 1L){
    step_message("Query via UniProtKB API.")
    db <- new_db(x$db_file, ".id")
    db <- not(db, object(x))
    res <- pbapply::pbsapply(db@query, simplify = F, cl = cl,
      function(query) {
        url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=",
          gs(query, " ", "%20"),
          "+reviewed:true+AND+organism_id:", x$organism_id,
          "&format=tsv")
        res <- try(data.table::fread(text = RCurl::getURL(url)), T)
        if (!inherits(res, "try-error")) {
          res
        } else NULL
      })
    res <- frbind(res, idcol = T, fill = T)
    db <- upd(db, res)
    raw_results <- as_tibble(dplyr::filter(db@db, .id %in% object(x)))
    fun_format <- function() {
      data <- dplyr::select(raw_results, query = .id, protein.names = `Protein names`,
        gene.names = `Gene Names`)
      message("Only keep the first row for formated results.")
      data <- dplyr::distinct(data, query, .keep_all = T)
      data <- split_lapply_rbind(data, 1:nrow(data),
        function(x) {
          data.frame(x, symbols = strsplit(x$gene.names, " ")[[1]])
        })
      data <- dplyr::select(data, -gene.names)
      data
    }
    format_results <- fun_format()
    p.query <- new_pie(object(x) %in% format_results$query, title = "Covered in UniprotKB")
    p.query <- .set_lab(p.query, sig(x), "Covered in UniprotKB")
    x@tables[[ 1 ]] <- namel(raw_results, format_results)
    x@plots[[ 1 ]] <- namel(p.query)
    return(x)
  })
