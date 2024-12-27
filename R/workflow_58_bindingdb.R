# ==========================================================================
# workflow of bindingdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_bindingdb <- setClass("job_bindingdb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp"),
    cite = "[@BindingdbIn20Gilson2016]",
    method = "The `BindingDB` database was used for discovering association between Ligands and Receptors",
    tag = "target",
    analysis = "BindingDB 药物靶点数据获取"
    ))

job_bindingdb <- function(file_db = .prefix("BindingDB_All_202401.tsv", "db"))
{
  if (!file.exists(file_db)) {
    stop("file.exists(file_db) == FALSE")
  }
  index <- ld_cutRead(file_db, c("PubChem CID", "UniProt (SwissProt) Primary ID of Target Chain"))
  colnames(index) %<>% strx("[^ ]+")
  .job_bindingdb(object = index)
}

setMethod("step0", signature = c(x = "job_bindingdb"),
  function(x){
    step_message("Prepare your data with function `job_bindingdb`.")
  })

setMethod("step1", signature = c(x = "job_bindingdb"),
  function(x, cids)
  {
    step_message("Filter to find targets.")
    res <- dplyr::filter(object(x), PubChem %in% !!cids, UniProt != "")
    info <- e(UniProt.ws::mapUniProt(
        from = "UniProtKB_AC-ID", to = "Gene_Name",
        columns = c("accession", "id"),
        query = list(ids = unique(res$UniProt))
        ))
    .add_internal_job(.job(method = "R package `UniProt.ws` used for querying Gene or Protein information"))
    res <- map(res, "UniProt", info, "From", "To", col = "symbol")
    res <- dplyr::filter(res, !is.na(symbol))
    x@tables[[1]] <- namel(res)
    object(x) <- NULL
    return(x)
  })


