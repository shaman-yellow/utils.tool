# ==========================================================================
# workflow of esearch
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_esearch <- setClass("job_esearch", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://www.ncbi.nlm.nih.gov/books/NBK179288/"),
    cite = ""
    ))

job_esearch <- function(key, jour = .jour_bioinf())
{
  data <- esearch.mj(key, jour)
  .job_esearch(object = data)
}

setMethod("vis", signature = c(x = "job_esearch"),
  function(x, ref = NULL, n = 100, col = c(".id", "Title")){
    data <- object(x)
    if (!is.null(ref)) {
      data <- dplyr::filter(data, grpl(.id, !!ref))
    }
    data <- dplyr::select(data, dplyr::all_of(col))
    print(data, n = n)
  })

setMethod("step0", signature = c(x = "job_esearch"),
  function(x){
    step_message("Prepare your data with function `job_esearch`.
      "
    )
  })
