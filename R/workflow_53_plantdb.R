# ==========================================================================
# workflow of plantdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_plantdb <- setClass("job_plantdb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://plantaedb.com/"),
    cite = "",
    method = "The Database `PlantaeDB` <https://plantaedb.com/> used for collating data of herbal ingredients"
    ))

job_plantdb <- function(url)
{
  .job_plantdb(object = url)
}

setMethod("step0", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Prepare your data with function `job_plantdb`.")
  })

setMethod("step1", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Collating data from URLs.")
    t.data <- lapply(object(x), function(x) get_plantaedb_data(x)$compounds)
    names(t.data) <- get_realname(object(x))
    t.data <- frbind(t.data, idcol = T)
    x@tables[[ 1 ]] <- namel(t.data)
    return(x)
  })

setMethod("step2", signature = c(x = "job_plantdb"),
  function(x){
    step_message("Obtaining targets.")
    # https://go.drugbank.com/releases/latest
    return(x)
  })

get_plantaedb_data <- function(url) {
  html <- RCurl::getURL(url)
  tables <- get_table.html(html)
  tables <- lapply(tables,
    function(x) {
      if (any(colnames(x) %in% c("Internal ID", "PubChem ID"))) {
        as_tibble(x)
      } else NULL
    })
  tables <- lst_clear0(tables)
  if (length(tables) == 2)
    names(tables) <- c("info", "compounds")
  else
    names(tables) <- "info"
  fun_format <- function(x) {
    seps <- grp(x$Name, "^>")
    group <- 0L
    index <- 0L
    seps <- vapply(1:nrow(x), FUN.VALUE = double(1),
      function(n) {
        if (is.na(seps[ index + 1L ]) || n > seps[ index + 1L ]) {
          group <<- group + 1L
          index <<- index + 1L
        } else if (n == seps[ index + 1L ]) {
          group <- 0L
        }
        group
      })
    x <- split(x, factor(seps, levels = unique(seps)))
    names <- x[[ "0" ]]$Name
    x <- x[-1]
    names(x) <- names
    x <- dplyr::rename(frbind(x, idcol = T), classes = .id)
    x <- dplyr::mutate(x, `Canonical SMILES` = gs(`Canonical SMILES`, "^Click to see\n\\s*", ""))
    x
  }
  if (nrow(tables$compounds))
    tables$compounds <- fun_format(tables$compounds)
  tables
}

get_bindingdb_data <- function(url = "https://www.bindingdb.org/bind/downloads/BindingDB_All_202401_tsv.zip",
  save = paste0("../", get_filename(url)),
  file_unzip = gs(gs(save, ".zip$", ""), "_tsv$", ".tsv"))
{
  if (!file.exists(file_unzip)) {
    if (!file.exists(save)) {
      cdRun("wget ", url, " -O ", save)
    }
    utils::unzip(save, exdir = get_path(save))
  }
}
