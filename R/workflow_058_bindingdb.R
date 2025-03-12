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

job_bindingdb <- function(dir = .prefix("bindingdb", "db"),
  file_db = file.path(dir, "BindingDB_All.tsv.gz"))
{
  if (!file.exists(file_db)) {
    stop("file.exists(file_db) == FALSE")
  }
  file_head <- add_filename_suffix(tools::file_path_sans_ext(file_db), "head")
  if (!file.exists(file_head)) {
    db_head <- ftibble(file_db, nrows = 5)
    data.table::fwrite(db_head, file_head, sep = "\t")
  } else {
    db_head <- ftibble(file_head)
  }
  index <- ld_cutRead(
    file_db, index_cols <- c(
      "BindingDB Reactant_set_id",
      "PubChem CID", "UniProt (SwissProt) Primary ID of Target Chain"
    ), 
    tmp = file.path(dir, "pubchem_uniprot.tsv")
  )
  index <- dplyr::mutate(index, .index = seq_len(nrow(index)), .before = 1)
  colnames(index) <- strx(colnames(index), "[^ ]+")
  x <- .job_bindingdb(object = index)
  x$file_db <- file_db
  x$db_head <- db_head
  x$index_cols <- index_cols
  return(x)
}

setMethod("step0", signature = c(x = "job_bindingdb"),
  function(x){
    step_message("Prepare your data with function `job_bindingdb`.")
  })

setMethod("step1", signature = c(x = "job_bindingdb"),
  function(x, cids = NULL, symbol = NULL)
  {
    step_message("Filter to find targets.")
    if (!is.null(cids) && is.null(symbol)) {
      res <- dplyr::filter(object(x), PubChem %in% !!cids, UniProt != "")
      info <- e(UniProt.ws::mapUniProt(
          from = "UniProtKB_AC-ID", to = "Gene_Name",
          columns = c("accession", "id"),
          query = list(ids = unique(res$UniProt))
          ))
      .add_internal_job(.job(method = "R package `UniProt.ws` used for querying Gene or Protein information"))
      res <- map(res, "UniProt", info, "From", "To", col = "symbol")
      res <- dplyr::filter(res, !is.na(symbol))
      x$res_filter <- res
      x@tables[[1]] <- namel(res)
      object(x) <- NULL
    } else if (!is.null(symbol)) {
      x$map_uniprot <- e(UniProt.ws::mapUniProt(
          from = "Gene_Name", to = "UniProtKB-Swiss-Prot",
          columns = c("accession", "id"),
          query = list(taxId = 9606, ids = unique(symbol))
          ))
      res <- dplyr::filter(object(x), UniProt %in% x$map_uniprot$Entry)
      x$res_filter <- res
      object(x) <- NULL
      methodAdd_onExit("x", "以 R 包 `UniProt.ws` ({packageVersion('UniProt.ws')}) 将基因名 (symbol) 获取 UniprotKB 数据库的 Entry ID。根据 Entry ID 检索 BindingDB 数据中，与该蛋白结合的化合物信息。")
      snapAdd_onExit("x", "在 BindingDB 数据中，检索与 {less(symbol)} 相关的化合物结合数据。")
    }
    x <- methodAdd(x, "获取 BindingDB 数据库化合物与靶点蛋白结合数据 (<https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp>) {cite_show('BindingdbIn20Gilson2016')}。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_bindingdb"),
  function(x, dir = "bindingdb"){
    step_message("Get all bindingdb infomation.")
    dir.create(dir, FALSE)
    file_rows <- file.path(dir, "rows.tsv")
    if (file.exists(file_rows)) {
      if (sureThat("Remove exists file?")) {
        file.remove(file_rows)
      }
    }
    x$annotation <- ld_slice_read(x$file_db, x$res_filter$.index + 1, file_rows)
    colnames(x$annotation) <- colnames(x$db_head)[ seq_len(ncol(x$annotation)) ]
    cols <- c(x$index_cols, "PDB ID(s) of Target Chain")
    x$annotation <- dplyr::relocate(x$annotation, dplyr::all_of(cols))
    colnames(x$annotation)[seq_along(cols)] <- strx(cols, "[^ ]+")
    x$annotation <- dplyr::select(x$annotation, -tidyselect::where(is.logical))
    x$annotation <- map(
      x$annotation, "UniProt", x$map_uniprot, "Entry", "From", col = "symbol"
    )
    x$annotation <- .set_lab(x$annotation, sig(x), "bindingdb annotation data")
    x$annotation <- setLegend(x$annotation, glue::glue("为 BindingDB 记录的化合物与蛋白结合信息附表。"))
    return(x)
  })

ld_slice_read <- function(file, rows, tmp = "/tmp/ldtmp.txt") {
  if (tools::file_ext(file) != "gz") {
    stop('tools::file_ext(file) != "gz".')
  }
  fileRows <- tempfile("rows", fileext = ".txt")
  writeLines(as.character(rows), fileRows)
  if (!file.exists(tmp)) {
    cdRun(
      glue::glue(
        "zcat {{{file}}} | awk 'NR==FNR {a[$1]; next} FNR in a' {{{fileRows}}} - > {{{tmp}}}",
        .open = "{{{", .close = "}}}"
      )
    )
  }
  ftibble(tmp)
}

setMethod("asjob_vina", signature = c(x = "job_bindingdb"),
  function(x){
    if (x@step < 2L) {
      stop('x@step < 2L.')
    }
    data <- x$annotation
    synonyms <- try_get_syn(data$PubChem)
    synonyms <- dplyr::mutate(synonyms, as.integer(synonyms$CID))
    data <- map(
      data, "PubChem", synonyms, "CID", "Synonym", col = "Synonym"
    )
    targets_annotation <- dplyr::select(data, symbol, PDB)
    targets_annotation <- reframe_col(
      targets_annotation, "PDB", function(x) unlist(strsplit(x, ","))
    )
    targets_annotation <- dplyr::distinct(targets_annotation)
    targets_annotation <- dplyr::rename(
      targets_annotation, hgnc_symbol = symbol, pdb = PDB
    )
    data <- dplyr::select(data, Synonym, symbol, PubChem)
    data <- dplyr::mutate(
      data, Synonym = ifelse(
        is.na(Synonym), paste0("CID:", PubChem), Synonym
      )
    )
    data <- dplyr::distinct(data)
    x <- job_vina(.layout = data)
    x$targets_annotation <- targets_annotation
    x <- snapAdd(x, "将 BindingDB 预筛选的化合物与靶点蛋白，用于分子对接。")
    return(x)
  })

