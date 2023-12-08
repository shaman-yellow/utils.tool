# ==========================================================================
# work and function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# asNamespace("base")

setClass("aplot")
setClassUnion("can_not_be_draw", c("recordedplot", "aplot"))

files <- setClass("files", 
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

fig <- setClass("fig", 
  contains = c("files"),
  representation = representation(),
  prototype = NULL)

f <- function(file) {
  .file_fig(file)
}

file_fig <- function(name, file) {
  nl(name, list(.file_fig(file)))
}

.file_fig <- setClass("file_fig", 
  contains = c("fig"),
  representation = representation(),
  prototype = NULL)

setMethod("show", signature = c(object = "files"),
  function(object){
    if (length(object) > 1) {
      lapply(object, function(x) browseURL(as.character(x)))
    } else {
      browseURL(as.character(object))
    }
  })

setClassUnion("figs", c("gg.obj", "fig"))
# setClassUnion("data.frame_or_matrix", c("data.frame", "matrix"))
setClassUnion("easywrite", c("data.frame", "matrix", "character", "factor", "numeric"))

# setClassUnion("maybe_numeric_or_character", c("numeric", "character", "missing"))
# setClassUnion("maybe_character", c("NULL", "character", "missing"))
# setClassUnion("maybe_logical", c("NULL", "logical", "missing"))

#' @importClassesFrom data.table data.table
#' @importClassesFrom tibble tbl_df
.df_like <- c("tbl_df", "data.table")
# lapply(.df_like, setClass, where = topenv())
.df <- c("data.frame", "matrix", "matrix", "array", .df_like)
setClassUnion("df", .df)

# <https://cran.r-project.org/web/packages/RSelenium/index.html>
format_bindingdb.tsv <- function(file,
  select = c("PubChem CID", "PDB ID(s) of Target Chain"), cl = 4, lines = NULL)
{
  .message_info("Read", file)
  data <- readLines(file)
  if (!is.null(lines)) {
    data <- data[ lines ]
  }
  .message_info("Split", "data")
  data <- strsplit(data, "\t")
  pos <- which(data[[1]] %in% select)
  pos <- pos[ 1:length(select) ]
  data <- pbapply::pblapply(data, function(ch) ch[ pos ], cl = cl)
  colnames <- data[[1]]
  .message_info("as.data.frame", "")
  data <- data.frame(t(dplyr::bind_cols(data[ -1 ])))
  data <- tibble::as_tibble(data)
  colnames(data) <- colnames
  namel(data, pos)
}

get_realname <- function(filename) {
  name <- get_filename(filename)
  gsub("\\..*$", "", name)
}

cdRun <- function(..., path = ".", sinkFile = NULL)
{
  owd <- getwd()
  setwd(normalizePath(path))
  expr <- paste0(unlist(list(...)), collapse = "")
  script <- tempfile("Script_", fileext = ".sh")
  writeLines(expr, script)
  writeLines(crayon::yellow(paste0("\nThe script file is: ", script)))
  if (is.null(sinkFile)) {
    tryCatch(system(expr), finally = setwd(owd))
  } else {
    tryCatch(
      capture.output(system(expr),
        file = sinkFile, split = T),
      finally = setwd(owd))
  }
}

# ==========================================================================
# autodock vina
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

vina_limit <- function(lig, recep, timeLimit = 120, dir = "vina_space", ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit, dir = dir), T)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoring = "ad4", stout = "/tmp/res.log", timeLimit = 60)
{
  if (!file.exists(dir)) {
    dir.create(dir, F)
  } 
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
  if (!file.exists(paste0(wd, "/", subdir, "_out.pdbqt"))) {
    dir.create(wd, F)
    file.copy(c(recep, lig), wd)
    .message_info("Generating affinity maps", subdir)
    .cdRun <- function(...) cdRun(..., path = wd)
    files <- get_filename(c(lig, recep))
    .cdRun("prepare_gpf.py -l ", files[1], " -r ", files[2], " -y")
    .cdRun("autogrid4 -p ", reals[2], ".gpf ", " -l ", reals[2], ".glg")
    cat("\n$$$$\n", date(), "\n", subdir, "\n\n", file = stout, append = T)
    try(.cdRun("timeout ", timeLimit, 
        " vina  --ligand ", files[1],
        " --maps ", reals[2],
        " --scoring ", scoring,
        " --exhaustiveness ", exhaustiveness,
        " --out ", subdir, "_out.pdbqt",
        " >> ", stout), T)
  }
}

vinaShow <- function(Combn, recep, subdir = Combn, dir = "vina_space",
  timeLimit = 3, backup = NULL)
{
  if (!file.exists(path <- paste0(dir, "/", subdir))) {
    stop("file.exists(path <- paste0(dir, \"/\", subdir))")
  } 
  wd <- paste0(dir, "/", subdir)
  out <- paste0(Combn, "_out.pdbqt")
  recep <- paste0(recep, ".pdbqt")
  res <- paste0(Combn, ".png")
  .cdRun <- function(...) cdRun(..., path = wd)
  try(.cdRun("timeout ", timeLimit, 
      " pymol ", out, " ", recep,
      " -g ", res), T)
  if (is.character(backup)) {
    dir.create(backup, F)
    file.copy(paste0(wd, "/", res), backup, T)
  }
}

summary_vina <- function(space = "vina_space", pattern = "_out\\.pdbqt$")
{
  files <- list.files(space, pattern, recursive = T, full.names = T)
  res_dock <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (length(lines) >= 1) {
        name <- gsub("_out\\.pdbqt", "", get_filename(file))
        top <- stringr::str_extract(lines[2], "[\\-0-9.]{1,}")
        top <- as.double(top)
        names(top) <- name
        top
      }
    })
  res_dock <- unlist(res_dock)
  res_dock <- tibble::tibble(
    Combn = names(res_dock), Affinity = unname(res_dock)
  )
  res_dock <- dplyr::mutate(
    res_dock, PubChem_id = stringr::str_extract(Combn, "^[^_]{1,}"),
    PDB_ID = stringr::str_extract(Combn, "[^_]{1,}$"),
    dir = paste0(space, "/", Combn),
    file = paste0(dir, "/", Combn, "_out.pdbqt")
  )
  dplyr::select(res_dock, PubChem_id, PDB_ID, Affinity, dir, file, Combn)
}

smiles_as_sdfs.obabel <- function(smiles) {
  lst.sdf <- pbapply::pbsapply(smiles, simplify = F,
    function(smi) {
      ChemmineOB::convertFormat("SMI", "SDF", source = test)
    })
  lst.sdf
}

tbmerge <- function(x, y, ...) {
  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)
  tibble::as_tibble(data.table::merge.data.table(x, y, ...))
}

ld_cols <- function(file, sep = "\t") {
  line <- readLines(file, n = 1)
  strsplit(line, sep)[[ 1 ]]
}

lst_clear0 <- function(lst, len = 0) {
  lst[ vapply(lst, function(v) if (length(v) > len) T else F, logical(1)) ]
}

ld_cutRead <- function(file, cols, abnum = T, sep = "\t", tmp = "/tmp/ldtmp.txt") {
  if (is.character(cols)) {
    names <- ld_cols(file, sep)
    pos <- which( names %in% cols )
    if (abnum) {
      pos <- head(pos, n = length(cols))
    }
  } else {
    pos <- cols
  }
  cdRun("cut -f ", paste0(pos, collapse = ","), " ", file, " > ", tmp)
  ftibble(tmp)
}

sdf_as_pdbqts <- function(sdf_file, mkdir.pdbqt = "pdbqt", check = F) {
  dir.create(mkdir.pdbqt, F)
  check_sdf_validity <- function(file) {
    lst <- sep_list(readLines(file), "^\\${4,}$")
    lst <- lst[ - length(lst) ]
    osum <- length(lst)
    lst <- lapply(lst,
      function(line) {
        pos <- grep("^\\s*-OEChem|M\\s*END", line)
        if (pos[2] - pos[1] > 5)
          line
      })
    lst <- lst[ !vapply(lst, is.null, logical(1)) ]
    nsum <- length(lst)
    list(data = lst, osum = osum, nsum = nsum, dec = osum - nsum)
  }
  if (check) {
    lst <- check_sdf_validity(sdf_file)
    lines <- unlist(lst$data)
    lst <- lst[ -1 ]
    writeLines(lines, sdf_file <- gsub("\\.sdf", "_modified.sdf", sdf_file))
  } else {
    lst <- list(file = sdf_file)
  }
  system(paste0("mk_prepare_ligand.py -i ", sdf_file, " --multimol_outdir ", mkdir.pdbqt))
  lst$file <- sdf_file
  lst$pdbqt <- list.files(mkdir.pdbqt, "\\.pdbqt$", full.names = T)
  lst$pdbqt.num <- length(lst$pdbqt)
  lst$pdbqt.cid <- stringr::str_extract(lst$pdbqt, "(?<=/|^)[0-9]{1,}")
  return(lst)
}

select_files_by_grep <- function(files, pattern){
  res <- lapply(files,
    function(file) {
      line <- readLines(file)
      if (any(grepl(pattern, line, T)))
        return(file)
      else NULL
    })
  res <- res[ !vapply(res, is.null, logical(1)) ]
  res
}

filter_pdbs <- function(files, pattern = "ORGANISM_SCIENTIFIC: HOMO SAPIENS") {
  select_files_by_grep(files, pattern)
}

prepare_receptor <- function(files, mkdir.pdbqt = "protein_pdbqt") {
  dir.create(mkdir.pdbqt, F)
  file <- lapply(files,
    function(file) {
      if (!is.null(file)) {
        newfile <- gsub("\\.pdb$", ".pdbqt", get_filename(file))
        newfile <- paste0(mkdir.pdbqt, "/", newfile)
        system(paste0("prepare_receptor -r ", file, " -o ", newfile, " -A hydrogens"))
        return(newfile)
      }
    })
  file <- file[ !vapply(file, is.null, logical(1)) ]
  file[ vapply(file, file.exists, logical(1)) ]
}

# sdfs_as_pdbqts <- function(lst.sdf, mkdir.sdf = "sdf", mkdir.pdbqt = "pdbqt")
# {
#   lapply(c(mkdir.sdf, mkdir.pdbqt), dir.create, showWarnings = F)
#   all_sdfs <- paste0(lst.sdf, collapse = "\n\n\n")
#   writeLines(all_sdfs, paste0(mkdir.sdf, "/all_compounds.sdf"))
#   n <- 1
#   pbapply::pbsapply(lst.sdf, 
#     function(sdf){
#       writeLines(sdf, file <- paste0(mkdir.sdf, "/", n, ".sdf"))
#       # system(paste0("mk_prepare_ligand.py -i ", file, " -o ", ))
#       n <<- n + 1
#     })
#   meta <- data.frame(
#     smiles = names(lst.sdf),
#     filename = paste0(1:length(smiles), ".sdf")
#   )
#   data.table::fwrite(meta, paste0(mkdir.sdf, "/metadata_of_filename.csv"))
#   data.table::fwrite(
#     dplyr::mutate(meta, filename = gsub("\\.sdf$", "\\.pdbqt$", filename)),
#     paste0(mkdir.pdbqt, "/metadata_of_filename.csv")
#   )
#   return(meta)
# }

# ==========================================================================
# crawl data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

start_drive <- function(command = "java -jar ~/operation/selenium.jar",
  port = 4444, extra = NULL, browser = c("chrome", "firefox"), ...)
{
  system(paste(command, "-port", port, extra), wait = F)
  new_link(port = port, browser = match.arg(browser))
}

end_drive <- function(pattern = "[0-9]\\s*java.*-jar") {
  tmp <- tempfile(fileext = ".txt")
  system(paste0("ps aux > ", tmp))
  text <- readLines(tmp)
  text <- text[ grepl(pattern, text) ]
  pid <- stringr::str_extract(text, "(?<=\\s)[0-9]{1,}(?=\\s)")
  system(paste("kill", paste0(pid, collapse = " ")))
}

download_herbTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Gene Targets")
{
  download_herbCompounds(link, ids, heading = heading)
}

download_compoundTargets <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Ingredient"),
  heading = "Related Gene Targets")
{
  download_herbCompounds(link, ids, urls, heading = heading)
}

download_herbCompounds <- function(link, ids,
  urls = paste0("http://herb.ac.cn/Detail/?v=", ids, "&label=Herb"),
  heading = "Related Ingredients")
{
  get_all <- function() list.files("~/Downloads", "\\.xlsx$", full.names = T)
  checks <- get_all()
  if (length(checks) > 0)
    file.remove(checks)
  link$open()
  jobn <- 0
  res <- lapply(urls,
    function(url) {
      link$navigate(url)
      jobn <<- jobn + 1
      if (jobn > 1) {
        link$refresh()
      }
      Sys.sleep(.5)
      res <- try(silent = T, {
        ele <- link$findElement(
          using = "xpath",
          value = paste0("//main//div//h4[text()='", heading, "']/../div//button")
        )
        Sys.sleep(.5)
        ele$sendKeysToElement(list("Download", key = "enter"))
      })
      Sys.sleep(.5)
      if (inherits(res, "try-error")) {
        writeLines("", paste0("~/Downloads/empty_", jobn, ".xlsx"))
      }
      checks <- get_all()
      cli::cli_alert_info(paste0("Job number ", jobn))
      while (length(checks) < jobn) {
        cli::cli_alert_info("Retry remotes response ...")
        ele$sendKeysToElement(list("Download", key = "enter"))
        Sys.sleep(1)
        checks <- get_all()
      }
    })
  link$close()
}

get_from_genecards <- function(query, score = 5, keep_drive = F) {
  query <- gs(query, " ", "%20")
  url <- paste0('https://www.genecards.org/Search/Keyword?queryString=', query,
    '&pageSize=25000&startPage=0')
  link <- start_drive(browser = "firefox")
  Sys.sleep(3)
  link$open()
  link$navigate(url)
  html <- link$getPageSource()[[1]]
  html <- XML::htmlParse(html)
  table <- XML::readHTMLTable(html)
  table <- as_tibble(data.frame(table[[1]]))
  colnames(table) %<>% gs("\\.+", "_")
  colnames(table) %<>% gs("X_|_$", "")
  table <- select(table, -1, -2)
  table <- filter(table, Score > !!score)
  if (!keep_drive)
    end_drive()
  return(table)
}

get_table.html <- function(file) {
  ht <- e(XML::htmlParse(readLines(file)))
  XML::readHTMLTable(ht)
}

moveToDir_herbs <- function(ids,
  file.pattern = "\\.xlsx$", index.pfun = file_seq.by_time,
  from = "~/Downloads", to = "herbs_ingredient", suffix = ".xlsx", .id = "herb_id", readFrom = NULL)
{
  if (is.null(readFrom)) {
    if (!file.exists(to)) {
      args <- as.list(environment())
      args$readFrom <- NULL
      files <- do.call(moveToDir, args)
      isOrdered <- T
    } else {
      files <- list.files(to, file.pattern, full.names = T)
      isOrdered <- F
    }
  } else {
    files <- readFrom
    isOrdered <- F
  }
  data <- lapply(files,
    function(file) {
      res <- try(openxlsx::read.xlsx(file), T)
      if (inherits(res, "try-error")) {
        tibble::tibble()
      } else res
    })
  if (isOrdered) {
    names(data) <- ids
  } else {
    names(data) <- get_realname(files)
  }
  data <- tibble::as_tibble(data.table::rbindlist(data, idcol = T, fill = T))
  cols <- colnames(data)
  data <- dplyr::relocate(data, .id)
  colnames(data)[ which(cols == ".id") ] <- .id
  data
}

moveToDir <- function(ids,
  file.pattern, index.pfun = file_seq.by_time,
  from, to, suffix, ...)
{
  files <- list.files(from, file.pattern, full.names = T)
  index <- index.pfun(files)
  files <- files[order(index)][1:length(ids)]
  files <- rev(files)
  dir.create(to, F, T)
  files <- lapply(1:length(ids),
    function(n) {
      file.copy(files[ n ], file <- paste0(to, "/", ids[n], suffix), T)
      file
    })
  files
}

file_seq.by_time <- function(files) {
  info <- fs::file_info(files)
  info <- dplyr::mutate(info, .INDEX = 1:nrow(info))
  info <- dplyr::arrange(info, dplyr::desc(modification_time))
  info <- dplyr::mutate(info, .SEQ = 1:nrow(info))
  info <- dplyr::arrange(info, .INDEX)
  info$.SEQ
}

format_index.by_num <- function(index) {
  index <- gsub("\\(|\\)", "", index)
  index <- ifelse(is.na(index), "0", index)
  as.integer(index)
}

new_link <- function(port = 4444L, browser = c("chrome", "firefox"), addr = "localhost")
{
  browser <- match.arg(browser)
  if (browser == "firefox") {
    prof <- RSelenium::makeFirefoxProfile(list('permissions.default.image' = 2L))
    RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = prof
    )
  } else if (browser == "chrome") {
    RSelenium::remoteDriver(
      remoteServerAddr = addr, port = port,
      browserName = browser,
      extraCapabilities = list()
    )
  }
}

merge.componentsGenes <- function(data, genes.lst) {
  names(genes.lst) <- data[[1]]
  genes.lst <- lapply(genes.lst,
    function(lst) {
      if (nrow(lst$ids) == 0) {
        return(NULL)
      } else {
        lst$ids$genes <- lst$data
        return(lst$ids)
      }
    })
}

get_tcm.components <- function(id, key = "Components",
  link_prefix = "http://www.tcmip.cn/ETCM/index.php/Home/Index/yc_details.html?id=",
  src = NULL)
{
  if (is.null(src)) {
    src <- get_tcm.base(id, link_prefix)
  }
  pos <- grep(paste0("^\\s*\\{\"ID*.*", key), src)
  data <- src[pos]
  ids.data <- extract_id(data)
  data <- gsub("<[^<^>]*>", "##", data)
  data <- gsub("^.*?####|####[^#]*$", "", data)
  data <- strsplit(data, "##, ##")[[1]]
  list(data = data, ids = ids.data, src = src)
}

extract_id <- function(ch) {
  links <- stringr::str_extract_all(ch, "(?<=href=\').*?(?=\')")[[1]]
  ids <- stringr::str_extract(links, "[^=]*$")
  data.frame(links = links, ids = ids)
}

get_tcm.componentToGenes <- function(ids, key = "Candidate Target Genes", dep = T) 
{
  link_prefix <- "http://www.tcmip.cn/ETCM/index.php/Home/Index/cf_details.html?id="
  lst <- pbapply::pblapply(ids,
    function(id) {
      data <- get_tcm.components(id, key, link_prefix)
      if (dep) {
        data$src <- NULL
      }
      data
    })
  lst
}

get_tcm.base <- function(id, link_prefix) {
  link <- paste0(link_prefix, id)
  src <- xml2::read_html(link)
  data <- rvest::html_text(src)
  data <- strsplit(data, "\r\n")[[1]]
  data
}

get_c2_data <- function(pattern = NULL,
  mode = c("hsa", "mmu"))
{
  mode <- match.arg(mode)
  path <- if (mode == "hsa") {
    "../human_c2_v5p2.rdata"
  } else if (mode == "mmu"){
    "../mouse_c2_v5p2.rdata"
  }
  if (!file.exists(path)) {
    ## download_url <- "https://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"
    stop("file.exists(path) == F")
  }
  name <- load(path)
  db <- get(name)
  names(db) <- tolower(names(db))
  if (!is.null(pattern))
    db <- db[ grep(pattern, names(db), ignore.case = T) ]
  return(db)
}

writeWraps <- function(lst, dir, width = 7, height = 7, ..., postfix = ".pdf") 
{
  fun <- function(p, file) write_graphics(p, file, mkdir = dir)
  writeDatas(lst, dir, fun = fun, postfix = postfix)
}

writePlots <- function(lst, dir, width = 7, height = 7, ..., postfix = ".pdf") 
{
  fun <- function(p, file) ggsave(file, p, width = width, height = height)
  writeDatas(lst, dir, fun = fun, postfix = postfix)
}

writeDatas <- function(lst, dir, ..., fun = data.table::fwrite, postfix = ".csv") 
{
  if (is.null(names(lst))) {
    stop("is.null(names(lst)) == T")
  }
  if (postfix %in% c(".pdf", ".png", ".jpg")) {
    super.dir <- get_savedir("figs")
  } else {
    super.dir <- get_savedir("tabs")
  }
  dir <- paste0(super.dir, "/", dir)
  dir.create(dir, F)
  n <- 1
  lapply(names(lst),
    function(name) {
      data <- lst[[name]]
      file <- paste0(dir, "/", n, "_", name, postfix)
      n <<- n + 1
      if (is.null(data) | is.character(data)) {
        writeLines(if (is.character(data)) data else "", gs(file, paste0(postfix, "$"), ".txt"))
      } else {
        fun(data, file)
      }
    })
  return(dir)
}

# ==========================================================================
# biomart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

list_datasets <- function() {
  ensembl <- e(biomaRt::useEnsembl("genes"))
  e(biomaRt::listDatasets(ensembl))
}

new_biomart <- function(dataset = c("hsapiens_gene_ensembl", "sscrofa_gene_ensembl", "mmusculus_gene_ensembl"))
{
  if (missing(dataset)) {
    message("Missing `dataset`, try get from global options.")
    arg <- getOption("mart_dataset", "hsapiens_gene_ensembl")
    dataset <- match.arg(arg, dataset)
    message("Use ", dataset)
  } else {
    dataset <- match.arg(dataset)
  }
  predb <- getOption("biomart")[[ dataset ]]
  if (is.null(predb)) {
    ensembl <- e(biomaRt::useEnsembl(biomart = "ensembl", dataset = dataset))
    lst <- nl(dataset, list(ensembl))
    options(biomart = lst)
  } else {
    ensembl <- predb
  }
  ensembl
}

filter_biomart <- function(mart, attrs, filters = "", values = "", distinct = T) {
  anno <- e(biomaRt::getBM(attributes = attrs, filters = filters,
    values = values, mart = mart))
  anno <- relocate(anno, !!rlang::sym(filters))
  if (distinct)
    anno <- distinct(anno, !!rlang::sym(filters), .keep_all = T)
  tibble::as_tibble(anno)
}

list_attrs <- function(mart) {
  tibble::as_tibble(biomaRt::listAttributes(mart))
}

general_attrs <- function(pdb = F, ensembl_transcript_id = F) {
  attrs <- c("ensembl_gene_id",
    "entrezgene_id",
    "hgnc_symbol",
    "refseq_mrna",
    "chromosome_name",
    "start_position",
    "end_position",
    "description")
  if (pdb) {
    attrs <- c(attrs, "pdb")
  }
  if (ensembl_transcript_id) {
    attrs <- c(attrs, "ensembl_transcript_id")
  }
  attrs
}

hsym <- function() {
  "hgnc_symbol"
}

get_herb_data <- function(herb = "../HERB_herb_info.txt",
  component = "../HERB_ingredient_info.txt",
  target = "../HERB_target_info.txt")
{
  db <- lapply(namel(herb, component, target),
    function(file) {
      tibble::as_tibble(suppressWarnings(data.table::fread(file)))
    })
  names <- colnames(db$component)
  colnames(db$component) <- c(names[-1], "extra_id")
  db
}

ftibble <- function(files, ...) {
  if (length(files) > 1) {
    lapply(files,
      function(file){
        tibble::as_tibble(data.table::fread(file, ...))
      })
  } else {
    tibble::as_tibble(data.table::fread(files, ...))
  }
}

fxlsx <- function(file, ...) {
  as_tibble(openxlsx::read.xlsx(file, ...))
}

fxlsx2 <- function(file, ..., .id = "sheet") {
  sheets <- openxlsx::getSheetNames(file)
  lst <- lapply(1:length(sheets),
    function(n) {
      openxlsx::read.xlsx(file, sheet = n)
    })
  names(lst) <- sheets
  data <- data.table::rbindlist(lst, idcol = T, fill = T)
  as_tibble(dplyr::rename(data, !!!nl(.id, ".id")))
}

get_nci60_data <- function(comAct = "../comAct_nci60/DTP_NCI60_ZSCORE.xlsx",
  rna = "../rna_nci60/RNA__RNA_seq_composite_expression.xls")
{
  comAct <- readxl::read_xlsx(comAct, skip = 8)
  rna <- readxl::read_xls(rna, skip = 10)
  rna <- dplyr::mutate(rna, genes = gsub("\\.[0-9]*$", "", `Gene name d`))
  list(comAct = comAct, rna = rna)
}

# ==========================================================================
# cor test
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mul_corTest <- function(data.x, data.y, namecol.x = 1, namecol.y = 1, method = "pearson",
  p.cutoff = .05, cor.cutoff = NULL, saveAllRes = T)
{
  key.x <- unlist(data.x[, namecol.x])
  key.y <- unlist(data.y[, namecol.y])
  inter <- colnames(data.x)[colnames(data.x) %in% colnames(data.y)]
  data.x <- dplyr::select(data.x, dplyr::all_of(inter))
  data.y <- dplyr::select(data.y, dplyr::all_of(inter))
  n <- 1
  lst <- pbapply::pbapply(data.x, 1, simplify = F,
    function(x) {
      message("\nCalculate data.x of row ", n, " to data.y")
      n <<- n + 1
      data <- pbapply::pbapply(data.y, 1, simplify = F,
        function(y) {
          cor <- cor.test(as.numeric(x), as.numeric(y), method = method)
          data.frame(cor = cor$estimate, p.value = cor$p.value)
        })
      data <- tibble::as_tibble(data.table::rbindlist(data))
      data <- dplyr::mutate(data, name = !!key.y)
      odata <- data <- dplyr::relocate(data, name)
      if (saveAllRes) {
        data.table::fwrite(odata, "pearsonTest_allResults.csv")
      }
      if (!is.null(p.cutoff)) {
        data <- dplyr::filter(data, p.value < p.cutoff)
      }
      if (!is.null(cor.cutoff)) {
        data <- dplyr::filter(data, cor > cor.cutoff)
      }
      data <- dplyr::arrange(data, p.value)
      attr <- apply(data.y, 1, simplify = F,  
        function(y) {
          data.frame(x = as.numeric(x), y = as.numeric(y))
        })
      names(attr) <- key.y
      attr <- data.table::rbindlist(attr, idcol = T)
      attr <- dplyr::rename(attr, group = .id)
      attr <- dplyr::filter(attr, group %in% data$name)
      attr <- dplyr::arrange(attr, factor(group, levels = data$name))
      attr(data, "data") <- tibble::as_tibble(attr)
      data
    })
  names(lst) <- key.x
  lst
}

vis_relcurve <- function(data, anno, rev = F, lab.x = "Expression (FPKM)",
  lab.y = "IC50 (Z-score)")
{
  if (rev) {
    cols <- colnames(data)
    cols.xy <- which(cols %in% c("x", "y"))
    colnames(data)[ cols.xy ] <- rev(cols[ cols.xy ])
  }
  anno <- merge(anno, cal_annoCoord(data), by.x = "name", by.y = ".id")
  anno <- dplyr::rename(anno, group = name)
  lm.data <- cal_lm(data)
  p <- ggplot(data, aes(x = x, y = y)) +
    geom_point(shape = 21, color = "transparent", fill = "#4DBBD5FF") +
    geom_line(data = lm.data, aes(x = x, y = y)) +
    geom_text(data = anno,
      aes(x = x, y = y,
        label = paste0("Cor = ", round(cor, 2), "\n", "P-value = ", round(p.value, 6))),
      size = 3, hjust = 0, vjust = 1) +
    labs(x = lab.x, y = lab.y) +
    facet_wrap(~ factor(group, levels = unique(data$group)), scales = "free")
  p
}

cal_lm <- function(data, col = "group") {
  fun <- function(num) {
    zoRange(num, 1)
  }
  fun2 <- function(value, range) {
    ifelse(value > range[1] & value < range[2], T, F)
  }
  lst <- split(data, data[[ col ]])
  lst <- lapply(lst,
    function(data) {
      coef <- lm(data$y ~ data$x)$coefficients
      x2y <- function(val) coef[1] + coef[2] * val
      data <- dplyr::mutate(data, .y = x2y(x))
      data <- dplyr::filter(data, fun2(.y, fun(y)))
      dplyr::select(data, -y, y = .y)
    })
  data.table::rbindlist(lst)
}

cal_annoCoord <- function(data, col = "group", pos.x = .1, pos.y = 1.3) {
  fun <- function(num, shift) {
    range <- zoRange(num, 1)
    range[1] + (range[2] - range[1]) * shift
  }
  data <- split(data, data[[ col ]])
  anno <- sapply(names(data), simplify = F,
    function(name) {
      data <- data[[ name ]]
      data.frame(x = fun(data$x, pos.x), y = fun(data$y, pos.y))
    })
  data.table::rbindlist(anno, idcol = T)
}

# ==========================================================================
# network
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

new_stringdb <- function(
  score_threshold = 200,
  species = 9606,
  network_type = c("physical", "full"),
  input_directory = "../",
  version = "11.5")
{
  e(STRINGdb::STRINGdb$new(score_threshold = score_threshold,
    species = species, network_type = match.arg(network_type),
    input_directory = input_directory, version = version
  ))
}

create_interGraph <- function(sdb, data, col = "name", rm.na = T) {
  cli::cli_alert_info("sdb$map")
  mapped <- sdb$map(data.frame(data), col, removeUnmappedRows = rm.na)
  message()
  cli::cli_alert_info("sdb$get_subnetwork")
  graph <- sdb$get_subnetwork(mapped$STRING_id)
  list(mapped = mapped, graph = graph)
}

fast_layout.str <- function(res.str, sdb, layout = "fr", seed = 10) {
  ## get annotation
  target <- paste0(sdb$input_directory, "/", sdb$species, ".protein.info.v",
    sdb$version, ".txt")
  dir.create(dir <- paste0(sdb$input_directory, "/temp"), F)
  if (!file.exists(file <- paste0(dir, "/", get_filename(target)))) {
    file.copy(gfile <- paste0(target, ".gz"), dir)
    R.utils::gunzip(paste0(file, ".gz"))
  }
  anno <- data.table::fread(file)
  ## merge with annotation
  igraph <- add_attr.igraph(res.str$graph, anno, by.y = "#string_protein_id")
  set.seed(seed)
  graph <- fast_layout(igraph, layout)
  attr(graph, "igraph") <- igraph 
  graph
}

cal_pagerank <- function(igraph) {
  res <- igraph::page_rank(igraph)$vector
  data <- data.frame(name = names(res), weight = unname(res))
  igraph <- add_attr.igraph(igraph, data, by.y = "name")
  igraph
}

get_nodes <- function(igraph, from = "vertices") {
  tibble::as_tibble(igraph::as_data_frame(igraph, from))
} 

dedup.edges <- function(igraph){
  add_attr.igraph(igraph)
}

add_attr.igraph <- function(igraph, data, by.x = "name", by.y, dedup.edges = T)
{
  comps <- igraph::as_data_frame(igraph, "both")
  if (!missing(data)) {
    nodes <- merge(comps$vertices, data, by.x = by.x, by.y = by.y, all.x = T)
  } else {
    nodes <- comps$vertices
  }
  if (dedup.edges) {
    edges <- dplyr::distinct(comps$edges)
  } else {
    edges <- comps$edges
  }
  igraph <- igraph::graph_from_data_frame(edges, vertices = nodes)
  igraph
}

output_graph <- function(igraph, file, format = "graphml", toCyDir = T) {
  igraph::write_graph(igraph, file, format = format)
}

plot_network.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, sc = 5, ec = 5, 
  arr.len = 2, edge.color = 'grey70', edge.width = .4, label = F)
{
  if (label) {
    layer.nodes <- geom_node_label(aes(label = preferred_name), size = label.size)
  } else {
    layer.nodes <- geom_node_point(aes(x = x, y = y, color = centrality_degree))
  }
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      color = edge.color, width = edge.width) +
    layer.nodes +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  p
} 

plot_networkFill.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, node.size = 12, sc = 5, ec = 5, 
  arr.len = 2, edge.color = 'lightblue', edge.width = 1, lab.fill = "MCC score",
  label = "genes")
{
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      # arrow = arrow(length = unit(arr.len, 'mm')),
      color = edge.color, width = edge.width) +
    geom_node_point(aes(x = x, y = y, fill = ifelse(is.na(MCC_score), 0, MCC_score)),
      size = node.size, shape = 21, stroke = .3) +
    geom_node_text(aes(label = !!rlang::sym(label)), size = label.size) +
    scale_fill_gradient(low = "lightyellow", high = "red") +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    labs(fill = "MCC score") +
    theme_void() +
    theme(plot.margin = margin(r = .05, unit = "npc")) +
    geom_blank()
  p
} 

multi_enrichKEGG <- function(lst.entrez_id, organism = 'hsa')
{
  res <- pbapply::pblapply(lst.entrez_id,
    function(ids) {
      res.kegg <- clusterProfiler::enrichKEGG(ids, organism = organism)
      res.path <- tibble::as_tibble(res.kegg@result)
      res.path <- dplyr::mutate(res.path, geneID_list = lapply(strsplit(geneID, "/"), as.integer))
      res.path
    })
  res
}

multi_enrichGO <- function(lst.entrez_id, orgDb = 'org.Hs.eg.db', cl = NULL)
{
  res <- pbapply::pblapply(lst.entrez_id, cl = cl,
    function(ids) {
      onts <- c("BP", "CC", "MF")
      res <- sapply(onts, simplify = F,
        function(ont) {
          res.go <- try(clusterProfiler::enrichGO(ids, orgDb, ont = ont), T)
          if (inherits(res.go, "try-error")) {
            return("try-error of enrichment")
          }
          res.res <- try(res.go@result, T)
          if (inherits(res.res, "try-error")) {
            value <- "try-error of enrichment"
            attr(value, "data") <- res.res
            return(value)
          }
          res.path <- tibble::as_tibble(res.res)
          res.path
        })
    })
}

check_enrichGO <- function(res.go) {
  isthat <- lapply(res.go,
    function(res) {
      !vapply(res, FUN.VALUE = logical(1), is.character)
    })
  isthat
}

as_double.ratioCh <- function(ch) {
  values <- stringr::str_extract_all(ch, "[0-9]{1,}")
  vapply(values, FUN.VALUE = double(1),
    function(values) {
      values <- as.double(values)
      values[ 1 ] / values[ 2 ]
    })
}

vis_enrich.kegg <- function(lst, cutoff_p.adjust = .1, maxShow = 10) {
  res <- lapply(lst,
    function(data) {
      data <- dplyr::filter(data, p.adjust < cutoff_p.adjust)
      if (nrow(data) == 0) return()
      data <- dplyr::arrange(data, p.adjust)
      data <- head(data, n = maxShow)
      data <- dplyr::mutate(data, GeneRatio = as_double.ratioCh(GeneRatio))
      p <- ggplot(data) +
        geom_point(aes(x = reorder(Description, GeneRatio),
            y = GeneRatio, size = Count, fill = p.adjust),
          shape = 21, stroke = 0, color = "transparent") +
        scale_fill_gradient(high = "yellow", low = "red") +
        scale_size(range = c(4, 6)) +
        labs(x = "", y = "Gene Ratio") +
        guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
        coord_flip() +
        ylim(zoRange(data$GeneRatio, 1.3)) +
        theme_minimal()
      p
    })
}

vis_enrich.go <- function(lst, cutoff_p.adjust = .1, maxShow = 10) {
  fun <- function(data) {
    data <- lapply(data,
      function(data) {
        if (is.character(data)) return(NULL)
        data <- dplyr::filter(data, p.adjust < cutoff_p.adjust)
        data <- dplyr::arrange(data, p.adjust)
        data <- head(data, n = maxShow)
        data
      })
    data <- data.table::rbindlist(data, idcol = T)
    data <- dplyr::mutate(
      data, GeneRatio = as_double.ratioCh(GeneRatio),
      stringr::str_wrap(Description, width = 30)
    )
    p <- ggplot(data) +
      geom_point(aes(x = reorder(Description, GeneRatio),
          y = GeneRatio, size = Count, fill = p.adjust),
        shape = 21, stroke = 0, color = "transparent") +
      scale_fill_gradient(high = "yellow", low = "red") +
      scale_size(range = c(4, 6)) +
      guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
      coord_flip() +
      facet_grid(.id ~ ., scales = "free") +
      theme_minimal() +
      theme(axis.title.y = element_blank(),
        strip.background = element_rect(fill = "grey90", color = "grey70")) +
      geom_blank()
    p
  }
  res <- lapply(lst,
    function(x) {
      try(fun(x), silent = T)
    })
  res
}

.element_textbox <- 
  function(family = NULL, face = NULL, size = NULL,
    colour = "white", fill = "lightblue",
    box.colour = "white", linetype = 1, linewidth = NULL,
    hjust = NULL, vjust = NULL,
    halign = 0.5, valign = NULL, lineheight = NULL,
    margin = match.fun("margin")(3, 3, 3, 3),
    padding = match.fun("margin")(2, 0, 1, 0),
    width = grid::unit(1, "npc"),
    height = NULL, minwidth = NULL,
    maxwidth = NULL, minheight = NULL, maxheight = NULL,
    r = grid::unit(5, "pt"), orientation = NULL,
    debug = FALSE, inherit.blank = FALSE
    ){
    structure(as.list(environment()),
      class = c("element_textbox", "element_text", "element"))
  }

get_savedir <- function(target = NULL) {
  res <- getOption("savedir")
  if (is.null(res$figs)) {
    res$figs <- "figs"
  }
  if (is.null(res$tabs)) {
    res$tabs <- "tabs"
  }
  if (!is.null(target)) {
    res <- res[[ target ]]
    if (!dir.exists(res))
      dir.create(res, F)
    res
  } else {
    res
  }
}

fwrite2 <- function(data, name, ..., file = paste0(get_realname(name), ".csv"),
  mkdir = get_savedir("tabs"), fun = data.table::fwrite)
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  fun(data, file <- paste0(mkdir, "/", file))
  return(file)
}

write_tsv2 <- function(data, name, ..., file = paste0(get_realname(name), ".tsv"),
  mkdir = get_savedir("tabs"), fun = write_tsv)
{
  do.call(fwrite2, as.list(environment()))
}

write_xlsx2 <- function(data, name, ..., file = paste0(get_realname(name), ".xlsx"),
  mkdir = get_savedir("tabs"), fun = openxlsx::write.xlsx)
{
  do.call(fwrite2, as.list(environment()))
}

write_graphics <- function(data, name, ..., file = paste0(get_realname(name), ".pdf"), page = -1,
  mkdir = get_savedir("figs"))
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  file <- paste0(mkdir, "/", file)
  if (is(data, "wrap")) {
    pdf(file, width = data@width, height = data@height)
  } else {
    pdf(file)
  }
  show(data)
  dev.off()
  len <- qpdf::pdf_length(file)
  if (len > 1) {
    if (page == -1)
      page <- len
    newfile <- tempfile(fileext = ".pdf")
    qpdf::pdf_subset(file, page, newfile)
    file.copy(newfile, file, T)
  }
  return(file)
}

write_character <- function(x, name, ..., file = paste0(get_realname(name), ".txt"),
  mkdir = get_savedir("figs"))
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  writeLines(x, file <- paste0(mkdir, "/", file))
  return(file)
}

write_gg <- function(p, name, width = 7, height = 7, ...,
  file = paste0(get_realname(name), ".pdf"), mkdir = get_savedir("figs")) 
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  ggsave(file <- paste0(mkdir, "/", file), p, width = width, height = height)
  return(file)
}

write_grob <- function(grob, name, width = 7, height = 7, ...,
  file = paste0(get_realname(name), ".pdf"), mkdir = get_savedir("figs")) 
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  pdf(file <- paste0(mkdir, "/", file), width = width, height = height)
  draw(grob)
  dev.off()
  return(file)
}

sortDup_edges <- function(edges) {
  edges.sort <- apply(dplyr::select(edges, 1:2), 1,
    function(vec) {
      sort(vec)
    })
  edges.sort <- tibble::as_tibble(data.frame(t(edges.sort)))
  edges.sort <- dplyr::distinct(edges.sort)
  edges.sort
}

getBelong_edges <- function(edges) {
  nodes <- unique(unlist(c(edges[, 1], edges[, 2])))
  edges.rev <- edges
  edges.rev[, 1:2] <- edges.rev[, 2:1]
  links.db <- rbind(edges, edges.rev)
  lst.belong <- split(links.db, unlist(links.db[, 1]))
  lst.belong <- lapply(lst.belong,
    function(data) unlist(data[, 2], use.names = F))
  lst.belong
}

cal_mcc.str <- function(res.str, name = "name", rename = T){
  hubs_score <- cal_mcc(res.str$graph)
  hubs_score <- tbmerge(res.str$mapped, hubs_score, by.x = "STRING_id", by.y = "name", all.x = T)
  hubs_score <- dplyr::relocate(hubs_score, !!rlang::sym(name), MCC_score)
  hubs_score <- dplyr::arrange(hubs_score, dplyr::desc(MCC_score))
  if (rename)
    hubs_score <- dplyr::rename(hubs_score, genes = name)
  hubs_score
}

cal_mcc <- function(edges) 
{
  if (is(edges, "igraph")) {
    igraph <- edges
    edges <- igraph::as_data_frame(edges, "edges")
  } else if (is(edges, "data.frame")) {
    igraph <- igraph::graph_from_data_frame(edges, F)
  }
  nodes <- unique(unlist(c(edges[, 1], edges[ , 2])))
  maxCliques <- igraph::max_cliques(igraph)
  scores <- vapply(nodes, FUN.VALUE = double(1), USE.NAMES = F,
    function(node) {
      if.contains <- vapply(maxCliques, FUN.VALUE = logical(1), USE.NAMES = F,
        function(clique) {
          members <- attributes(clique)$names
          if (any(members == node)) T else F
        })
      in.cliques <- maxCliques[ if.contains ]
      scores <- vapply(in.cliques, FUN.VALUE = double(1),
        function(clique) {
          num <- length(attributes(clique)$names)
          factorial(num - 1)
        })
      sum(scores)
    })
  res <- data.frame(name = nodes, MCC_score = scores)
  res <- tibble::as_tibble(dplyr::arrange(res, dplyr::desc(MCC_score)))
  res
}

get_subgraph.mcc <- function(igraph, resMcc, top = 10) 
{
  tops <- dplyr::arrange(resMcc, dplyr::desc(MCC_score))
  tops <- head(tops$STRING_id, n = top)
  data <- igraph::as_data_frame(igraph, "both")
  nodes <- dplyr::filter(data$vertices, name %in% !!tops)
  nodes <- merge(nodes, resMcc, by.x = "name", by.y = "STRING_id", all.x = T)
  edges <- dplyr::filter(data$edges, (from %in% !!tops) & (to %in% !!tops))
  igraph <- igraph::graph_from_data_frame(edges, F, nodes)
  igraph <- dedup.edges(igraph)
  igraph
}

search_resFile <- function(file = "index.Rmd",
  pattern.field = "\\*\\*\\(对应文件.*?\\)\\*\\*",
  pattern.file = "(?<=`).*?(?=`)")
{
  md <- readLines(file)
  md <- paste0(md, collapse = "\n")
  matchs <- stringr::str_extract_all(md, pattern.field)[[1]]
  files <- unlist(stringr::str_extract_all(matchs, pattern.file))
  files <- files[ !grepl("^,", files) ]
  res <- lapply(files,
    function(file) {
      if (!file.exists(file))
        stop("file.exists(", file, ") == F")
      file
    })
  unlist(res)
}

package_results <- function(
  master_include = NULL,
  report = "output.pdf",
  main = unique(c(search_resFile("index.Rmd"), autoRegisters)), 
  headPattern = "[^r][^e][^s].*业务.*\\.zip", head = list.files(".", headPattern),
  masterSource = c("output.Rmd", "install.R", "read_me.txt", head),
  prefix = "results_", zip = paste0(prefix, head),
  clientZip = "client.zip", masterZip = "master.zip",
  clear = T, external_file = "order_material")
{
  if (!file.exists(external_file)) {
    stop("file.exists(external_file)")
  } else {
    if (length(exters <- list.files(external_file, full.names = T)) == 0) {
      stop("Guess you have forget to save the order material file to `external_file`.")
    } else {
      dir.create(extern0 <- "报单相关资料", F)
      file.copy(exters, extern0, T, T)
    }
  }
  files.client <- c(report, main, extern0)
  files.master <- c(files.client, head, masterSource)
  zip(clientZip, files.client)
  if (!is.null(masterZip)) {
    zip(masterZip, files.master)
    zip(zip, c(clientZip, masterZip), flags = "-ur9X")
    if (clear) {
      file.remove(c(clientZip, masterZip))
    }
  } else {
    message("Only client file were packaged.")
  }
}

colSum <- function(col) {
  length(unique(col))
}

show_multi <- function(layers, col, symbol = "Progein"){
  cat(crayon::silver("Data of", length(layers), "\n"))
  res <- mapply(layers, 1:length(layers), names(layers),
    FUN = function(com, seq, name){
      cat(crayon::silver("  +++ ", symbol, " ", seq, "+++\n"))
      cat("  ", crayon::yellow(name), "\n", rep(" ", 4), "Sum: ", nrow(com),
        "\n", sep = "")
      if (is.null(com)) {
        cat(paste0(rep(" ", 6), collapse = ""), crayon::silver("No data.\n"), sep = "")
        return()
      }
      args <- com[[ col ]]
      if (length(args) <= 4) {
        cat(paste0(paste0(rep(" ", 6), collapse = ""),
            col, ": ", paste0(args, collapse = ", ")), sep = "\n")
      } else {
        cat(paste0(paste0(rep(" ", 6), collapse = ""), col, ": ",
            paste0(c(args[1:4], "..."), collapse = ", ")), sep = "\n")
      }
      cat("\n")
    })
}

# map_from <- function(lst, db_names, db_values){}

.jour_bioinf <- function() {
  c("Nature Genetics", "Genome Biology", "Genome Research", "Nucleic Acids Res",
    "Briefings in Bioinformatics", "BMC genomics", "Nature Methods", "Nature Biotechnology")
}

esearch.mj <- function(key, jour = .jour_bioinf(), rbind = T)
{
  query <- paste(jour, "[JOUR] AND", key)
  sear <- pbapply::pblapply(query, esearch)
  names(sear) <- jour
  if (rbind) {
    sear <- tibble::as_tibble(data.table::rbindlist(sear, idcol = T, fill = T))
    if (nrow(sear) > 0)
      sear <- dplyr::arrange(sear, dplyr::desc(SortPubDate))
  }
  sear
}

esearch <- function(query = NULL, fetch.save = paste0(gsub(" ", "_", query), ".xml"),
  path = "search", tract.save = "res.tsv",
  fields = c("SortPubDate", "Title", "FullJournalName", "Name", "Id"))
{
  if (!dir.exists(path)) {
    dir.create(path)
  }
  if (!is.null(query)) {
    if (!file.exists(paste0(path, "/", fetch.save))) {
      cdRun("esearch -db pubmed -query \"", query, "\"",
        " | efetch -format docsum > ", fetch.save,
        path = path)
    }
  } else if (fetch.save == ".xml") {
    stop("fetch.save == \".xml\"")
  }
  cdRun(" cat ", fetch.save, " | xtract -pattern DocumentSummary ",
    " -sep \"|\" -element ", paste0(fields, collapse = " "),
    " > ", tract.save, path = path)
  file <- paste0(path, "/", tract.save)
  res <- ftibble(file, sep = "\t", quote = "", fill = T, header = F)
  if (nrow(res > 0)) {
    colnames(res) <- fields
    if (any("SortPubDate" == fields)) {
      res <- dplyr::mutate(res, SortPubDate = as.Date(SortPubDate))
      res <- dplyr::arrange(res, dplyr::desc(SortPubDate))
    }
    if (any("Id" == fields)) {
      res <- dplyr::mutate(res, Id = as.character(Id))
    }
  }
  res
}

split_lapply_rbind <- function(data, f, fun, ..., verbose = F) {
  data <- split(data, f)
  if (verbose)
    data <- pbapply::pblapply(data, fun, ...)
  else
    data <- lapply(data, fun, ...)
  data <- data.table::rbindlist(data)
  tibble::as_tibble(data)
}

.data_long <- setClass("data_long", 
  contains = c("data.frame"),
  representation = representation(),
  prototype = NULL)

setMethod("show", signature = c(object = "data_long"),
  function(object){
    suppressWarnings(print(tibble::as_tibble(object)))
  })

handling_na.long <- function(data) {
  .check_columns(data, c("sample", "group"))
  metadata <- dplyr::select(data, sample, group)
  data.wide <- dplyr::select(data, -group)
  data.wide <- tidyr::gather(data.wide, "var", "value", -sample)
  data.wide <- tidyr::spread(data.wide, sample, value)
  data.wide <- handling_na(data.wide, "var", metadata)
  data.long <- tidyr::gather(data.wide, "sample", "value", -var)
  data.long <- dplyr::mutate(data.long, group = metadata$group)
  dplyr::relocate(data.long, sample, group)
}

variance_analysis <- function(data.long){
  .check_columns(data.long, c("sample", "group"))
  ## shapiro
  p.shapiro <- lapply(split(data.long, ~ group),
    function(data){
      shapiro.test(data$value)$p.value
    })
  p.shapiro <- as_df.lst(p.shapiro, "name", "p.value")
  p.shapiro$type <- "Shapiro"
  ## bartlett
  p.bartlett <- bartlett.test(data.long$value, data.long$group)$p.value
  p.bartlett <- data.frame(type = "Bartlett", name = "All", p.value = p.bartlett)
  ## aov
  aov <- aov(value ~ group, data = data.long)
  p.aov <- unlist(summary(aov))[[ "Pr(>F)1" ]]
  p.aov <- data.frame(type = "Variance.", name = "All", p.value = p.aov)
  ## lsd or hsd
  res <- lapply(c("hsd", "lsd"),
    function(method){
      res <- DescTools::PostHocTest(aov, method = method)$group
      data.frame(type = switch(method, hsd = "Tukey HSD", lsd = "Fisher LSD"),
        name = rownames(res), p.value = data.frame(res)$pval)
    })
  res <- data.table::rbindlist(res)
  all.p <- rbind(p.shapiro, p.bartlett, p.aov, res)
  all.p
}

.split_data.long <- function(data.long){
  .check_columns(data.long, c("sample", "group"))
  metadata <- dplyr::select(data.long, group, sample)
  metadata <- tibble::as_tibble(metadata)
  data <- dplyr::select(data.long, -group, -sample)
  rownames(data) <- metadata$sample
  namel(data, metadata)
}

# ==========================================================================
# pca and oplsda
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.andata <- setClass("andata", 
  contains = character(),
  representation = 
    representation(
      data = "ANY",
      anno = "ANY",
      metadata = "ANY",
      palette = "ANY",
      raw = "ANY",
      extra = "ANY"),
    prototype = NULL
)

setGeneric("pal", 
  function(x) standardGeneric("pal"))

setMethod("pal", signature = c(x = "andata"),
  function(x){ 
    pal <- x@palette
    if (is.null(pal)) {
      group <- unique(x@metadata$group)
      pal <- nl(group, color_set()[1:length(group)], F)
    }
    pal
  })

.andata_pca <- setClass("andata_pca", contains = c("andata"))

pca_data.long <- function(x, fun_scale = function(x) scale(x)) {
  lst <- .split_data.long(x)
  data <- t(lst$data)
  metadata <- lst$metadata
  if (!is.null(fun_scale)) {
    data <- fun_scale(data)
  }
  data <- prcomp(data, retx = F)
  anno <- round(summary(data)$importance[2, ], 3)
  data <- tibble::as_tibble(data$rotation)
  data <- dplyr::mutate(data, sample = metadata$sample,
    group = metadata$group)
  data <- dplyr::select(data, sample, group, PC1, PC2)
  .andata_pca(data = data, anno = anno, metadata = metadata)
}

.andata_opls <- setClass("andata_opls", 
  contains = c("andata"),
  representation = representation(vip = "ANY"), prototype = NULL)

opls_data.long <- function(data.long, combns, inter.fig = F) {
  lst <- .split_data.long(data.long)
  data <- lst$data
  metadata <- lst$metadata
  if (!is.list(combns)) {
    stop("is.list(combns) == F")
  }
  res <- lapply(combns,
    function(combn){
      whi <- which(metadata$group %in% combn)
      metadata <- metadata[whi, ]
      data <- data[whi, ]
      opls <- ropls::opls(
        data, as.character(metadata$group),
        predI = 1, orthoI = NA, fig.pdfC = F
      )
      data <- cbind(opls@scoreMN[, 1], opls@orthoScoreMN[, 1])
      colnames(data) <- c("h1", "o1")
      data <- tibble::as_tibble(data)
      anno <- c(
        x = paste0("T score[1](", opls@modelDF[1, "R2X"] * 100, "%)"),
        y = paste0("Orthogonal T score[1](", opls@modelDF[2, "R2X"] * 100, "%)")
      )
      vip <- data.frame(opls@vipVn)
      vip <- cbind(rownames(vip), vip)
      colnames(vip) <- c("var", "vip")
      vip <- tibble::as_tibble(vip)
      .andata_opls(data = data, anno = anno, metadata = metadata, vip = vip)
    })
  names(res) <- paste0("combn", 1:length(combns))
  res
}

setGeneric("plot_andata", 
  function(x, ...) standardGeneric("plot_andata"))

setMethod("plot_andata", signature = c(x = "andata_pca"),
  function(x){
    p <- ggplot(x@data, aes(x = PC1, y = PC2, fill = group)) +
      geom_point(size = 3, shape = 21, stroke = 0, color = "transparent") +
      stat_ellipse(aes(color = group), level  =  0.95) +
      labs(x = paste0("PC1 (", x@anno[1] * 100, "%)"),
        y = paste0("PC2 (", x@anno[2] * 100, "%)"), color = "Group", fill = "Group") +
      scale_fill_manual(values = pal(x)) +
      scale_color_manual(values = pal(x)) +
      theme()
    p
  })

setMethod("plot_andata", signature = c(x = "andata_opls"),
  function(x){
    p <- ggplot(x@data, aes(x = h1, y = o1)) +
      geom_point(size = 3, shape = 21,
        aes(fill = x@metadata$group), stroke = 0, color = "transparent") +
      stat_ellipse(aes(color = x@metadata$group), level  =  0.95) +
      scale_fill_manual(values = pal(x)) +
      scale_color_manual(values = pal(x)) +
      labs(x = x@anno[1], y = x@anno[2], color = "Group", fill = "Group") +
      theme()
    p
  })

get_metadata.geo <- function(lst,
  select = rlang::quos(rownames, title),
  pattern = c("diagnosis", "Sex", "^age", "^time point", "data_processing"),
  abbrev = c("data_processing"))
{
  res <- lapply(lst,
    function(eset){
      as_tibble(eset@phenoData@data)
    })
  if (!is.null(select)) {
    main <- lapply(res,
      function(data){
        cols <- colnames(data)
        extra <- unlist(.find_and_sort_strings(cols, pattern))
        dplyr::select(data, !!!select, dplyr::all_of(extra))
      })
    if (!is.null(abbrev)) {
      abbrev <- lapply(res,
        function(data){
          cols <- colnames(data)
          extra <- unlist(.find_and_sort_strings(cols, abbrev))
          dplyr::distinct(data, dplyr::pick(extra))
        })
    } else abbrev <- NULL
    res <- namel(main, abbrev, res)
  }
  lst_clear0(res)
}

show_lst.ch <- function(lst, width = 60) {
  sapply(names(lst),
    function(name){
      message("+++ ", name, " +++\n")
      textSh(lst[[ name ]], wrap_width = width, pre_wrap = T)
    })
  message()
}

get_prod.geo <- function(lst) {
  res <- as.list(dplyr::distinct(do.call(rbind, lst$abbrev)))
  .lich(res)
}

.lich <- setClass("lich", 
  contains = c("list"),
  representation = 
    representation(),
  prototype = NULL)

new_lich <- function(lst) {
  lst <- lapply(lst, paste0, collapse = ", ")
  .lich(lst)
}

setMethod("show", signature = c(object = "lich"),
  function(object){
    show_lst.ch(object)
  })

# ==========================================================================
# limma and edgeR
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

prepare_expr_data <- function(metadata, counts, genes, message = T) 
{
  ## sort and make name for metadata and counts
  checkDup <- function(x) {
    if (any(duplicated(x))) {
      stop("The value in ID column are duplicated.")
    }
  }
  lapply(list(counts[[ 1 ]], genes[[ 1 ]], metadata[[ 1 ]]), checkDup)
  colnames(metadata) %<>% make.names()
  metadata <- dplyr::rename(metadata, sample = 1)
  counts <- dplyr::select(counts, 1, dplyr::all_of(metadata$sample))
  metadata$sample %<>% make.names()
  colnames(counts) %<>% make.names()
  ## sort genes
  data_id <- do.call(data.frame, nl(colnames(genes)[1], list(counts[[1]])))
  genes <- dplyr::distinct(genes, !!rlang::sym(colnames(genes)[1]), .keep_all = T)
  genes <- tbmerge(
    data_id, genes,
    by.x = colnames(data_id)[1], by.y = colnames(genes)[1],
    sort = F, all.x = T
  )
  if (ncol(genes) > 1) {
    if (message) {
      message("## The missing genes in `genes`")
      check.na <- dplyr::filter(
        genes, is.na(!!rlang::sym(colnames(genes)[2]))
      )
      print(check.na)
    }
  }
  counts <- dplyr::select(counts, -1)
  colnames(counts) <- metadata$sample
  namel(counts, metadata, genes)
}

new_dge <- function(metadata, counts, genes, message = T)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  e(edgeR::DGEList(lst$counts, samples = lst$metadata, genes = lst$genes))
}

new_elist <- function(metadata, counts, genes, message = T)
{
  lst <- do.call(prepare_expr_data, as.list(environment()))
  elist(list(E = tibble::as_tibble(lst$counts), targets = lst$metadata, genes = lst$genes))
}

# #' @importClassesFrom edgeR DGEList
# #' @importClassesFrom limma EList
.eset <- c("DGEList", "EList")
lapply(.eset, setClass, where = topenv())
setClassUnion("eset", .eset)

elist <- setClass("elist", 
  contains = c("EList", "list"),
  representation = representation(),
    prototype = NULL)

setValidity("elist", 
  function(object){
    if (nrow(object$E) != nrow(object$genes))
      F
    else T
  })

.data_long.eset <- setClass("data_long.eset", 
  contains = c("data_long"),
  representation = representation(),
  prototype = NULL)

setGeneric("as_data_long", 
  function(x, y, ...) standardGeneric("as_data_long"))

setMethod("as_data_long", signature = c(x = "DGEList"),
  function(x){
    data <- tibble::as_tibble(x$counts)
    data$gene <- x$genes[[ 1 ]]
    data <- tidyr::gather(data, sample, value, -gene)
    if (!all(unique(data$sample) %in% x$samples$sample)) {
      stop("!all(unique(data$sample) %in% x$samples$sample) == T")
    }
    data <- tbmerge(data, dplyr::select(x$samples, sample, group),
      by = "sample", all.x = T)
    .data_long.eset(data)
  })

setMethod("as_data_long", signature = c(x = "EList"),
  function(x){
    data <- tibble::as_tibble(x$E)
    data$gene <- x$genes[[ 1 ]]
    data <- tidyr::gather(data, sample, value, -gene)
    data <- tbmerge(data, x$targets, by = "sample", all.x = T)
    .data_long.eset(data)
  })

setGeneric("pca_data.long", 
  function(x, fun_scale) standardGeneric("pca_data.long"))

setMethod("pca_data.long", signature = c(x = "data_long.eset"),
  function(x){
    x <- dplyr::select(tibble::as_tibble(data.frame(x)), sample, group, gene, value)
    x <- tidyr::spread(x, gene, value)
    callNextMethod(x, NULL)
  })

filter_low.dge <- function(dge, group., min.count = 10, prior.count = 2) {
  dge$samples$group = group.
  mean <- mean(dge$samples$lib.size) * 1e-6
  median <- median(dge$samples$lib.size) * 1e-6
  cutoff <- log2(min.count / median + 2 / mean)
  ## raw
  raw_dge <- dge
  raw_dge$counts <- edgeR::cpm(raw_dge, log = T, prior.count = prior.count)
  ## filter
  keep.exprs <- e(edgeR::filterByExpr(dge, group = group., min.count = min.count))
  pro_dge <- dge <- e(edgeR::`[.DGEList`(dge, keep.exprs, , keep.lib.sizes = F))
  pro_dge$counts <- e(edgeR::cpm(pro_dge, log = T))
  ## plot
  data <- list(Raw = as_data_long(raw_dge), Filtered = as_data_long(pro_dge))
  data <- data.table::rbindlist(data, idcol = T)
  data <- dplyr::select(data, .id, sample, value)
  p <- ggplot(data) +
    geom_density(aes(x = value, color = sample), alpha = .7) +
    geom_vline(xintercept = cutoff, linetype = "dashed") +
    facet_wrap(~ factor(.id, c("Raw", "Filtered"))) +
    labs(x = "Log2-cpm", y = "Density")
  if (length(unique(data$sample)) > 50) {
    p <- p + theme(legend.position = "none")
  }
  attr(dge, "p") <- p
  dge
}

mx <- function(...){
  design <- model.matrix(...)
  colnames(design) %<>% gsub(".*?group\\.?", "", .)
  colnames(design) %<>% gsub(".*?batch\\.?", "batch.", .)
  design
}

norm_genes.dge <- function(dge, design, prior.count = 2, fun = limma::voom, ..., vis = T){
  ## raw
  raw_dge <- dge
  raw_dge$counts <- edgeR::cpm(raw_dge, log = T, prior.count = prior.count)
  ## pro
  dge <- e(edgeR::calcNormFactors(dge, method = "TMM"))
  pro_dge <- dge <- fun(dge, design, ...)
  ## data long
  if (vis) {
    cli::cli_alert_info("as_data_long")
    data <- list(Raw = as_data_long(raw_dge), Normalized = as_data_long(pro_dge))
    data <- data.table::rbindlist(data, idcol = T, fill = T)
    data <- dplyr::select(data, .id, sample, value)
    if (length(unique(data$sample)) < 50) {
      p <- ggplot(data) +
        geom_boxplot(aes(x = sample, y = value),
          outlier.color = "grey60", outlier.size = .5) +
        coord_flip() +
        facet_wrap(~ factor(.id, c("Raw", "Normalized"))) +
        labs(x = "Sample", y = "Log2-cpm")
    } else {
      p <- plot_median_expr_line(data)
    }
  } else {
    p <- NULL
  }
  attr(dge, "p") <- p
  dge
}

diff_test <- function(x, design, contr = NULL, block = NULL){
  if (!is.null(block)){
    dupcor <- e(limma::duplicateCorrelation(x, design, block = block))
    cor <- dupcor$consensus.correlation
    message("## Within-donor correlation:", cor)
  } else {
    cor <- NULL
  }
  fit <- e(limma::lmFit(x, design, block = block, correlation = cor))
  if (!is.null(contr)) {
    fit <- e(limma::contrasts.fit(fit, contrasts = contr))
  }
  fit <- e(limma::eBayes(fit))
  fit
}

extract_tops <- function(x, use = "adj.P.Val", use.cut = 0.05, cut.fc = 0.3){
  res <- e(lapply(1:ncol(x$contrasts),
    function(coef){
      res <- limma::topTable(x, coef = coef, number = Inf)
      if (!is.null(use.cut) & !is.null(cut.fc)) {
        res <- dplyr::filter(res, !!rlang::sym(use) < use.cut, abs(logFC) > cut.fc)
      }
      as_tibble(res) 
    }))
  names(res) <- colnames(x$contrasts)
  res
}

limma_downstream <- function(dge.list, group., design, contr.matrix,
    min.count = 10, voom = T, cut.q = 0.05, cut.fc = 0.3,
    get_ebayes = F, get_normed.exprs = F, block = NULL)
  {
    # if (NA %in% dge.list$samples[["lib.size"]]){
    # dge.list$samples[["lib.size"]] <- apply(dge.list$counts, 2, sum, na.rm = T)
    # }
    keep.exprs <- e(edgeR::filterByExpr(dge.list, group = group., min.count = min.count))
    dge.list <- edgeR::`[.DGEList`(dge.list, keep.exprs, , keep.lib.sizes = F)
    if (voom){
      dge.list <- e(edgeR::calcNormFactors(dge.list, method = "TMM"))
      dge.list <- e(limma::voom(dge.list, design))
    }else{
      genes <- dge.list$genes
      targets <- dge.list$samples
      dge.list <- scale(dge.list$counts, scale = F, center = T)
    }
    if (get_normed.exprs)
      return(dge.list)
    if (!is.null(block)){
      dupcor <- e(limma::duplicateCorrelation(dge.list, design, block = block))
      cor <- dupcor$consensus.correlation
      cat("## Within-donor correlation:", cor, "\n")
    }else{
      cor <- NULL
    }
    fit <- e(limma::lmFit(dge.list, design, block = block, correlation = cor))
    if (!voom){
      fit$genes <- genes
      fit$targets <- targets
    }
    fit.cont <- e(limma::contrasts.fit(fit, contrasts = contr.matrix))
    ebayes <- e(limma::eBayes(fit.cont))
    if (get_ebayes)
      return(ebayes)
    res <- e(lapply(1:ncol(contr.matrix),
        function(coef){
          results <- limma::topTable(ebayes, coef = coef, number = Inf) %>% 
            dplyr::filter(adj.P.Val < cut.q, abs(logFC) > cut.fc) %>% 
            dplyr::as_tibble() 
          return(results)
        }))
    names(res) <- colnames(contr.matrix)
    return(res)
  }

fuzzy <- function(str) {
  make.names(gsub(" ", "", str))
}

# ==========================================================================
# heatmap
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

get_fun <- function(name, envir = topenv()) {
  get(name, envir = envir, mode = "function")
}

.heatdata <- setClass("heatdata", 
  contains = c("can_not_be_draw"),
  representation = 
    representation(
      raw = "ANY",
      data_long = "ANY",
      main = "ANY", aesn = "ANY", para = "ANY",
      aesh = "ANY",
      ymeta = "ANY", y_aesn = "ANY", y_pal = "ANY",
      xmeta = "ANY", x_aesn = "ANY", x_pal = "ANY",
      gg_main = "ANY",
      gg_xtree = "ANY", gg_ytree = "ANY",
      gg_xgroup = "ANY", gg_ygroup = "ANY",
      fun_plot = "ANY"),
    prototype = prototype(
      aesh = "fill",
      aesn = list(
        x = "cname", y = "rname", fill = "value",
        lab_x = "Columns", lab_y = "Rows", lab_fill = "Level"),
      x_aesn = list(x = "cname", y = "Group", fill = "group",
        lab_x = "", lab_y = "", lab_fill = "Column group"),
      x_pal = color_set(),
      y_aesn = list(x = "Group", y = "rname", fill = "group",
        lab_x = "", lab_y = "", lab_fill = "Row group"),
      para = list(
        clust_row = T, clust_col = T, method = 'average'),
      fun_plot = list(
        xtree = get_fun("plot_xtree"),
        ytree = get_fun("plot_ytree"),
        main = get_fun("tile_heatmap"),
        xgroup = get_fun("add_xgroup.tile.heatmap"),
        ygroup = get_fun("add_ygroup.tile.heatmap")
      )
      ))

setMethod("show", signature = c(object = "heatdata"),
  function(object){
    if (!is.null(object@gg_main))
      print(draw(object))
    else
      message("A 'heatdata' object")
  })

heatdata_gene <- setClass("heatdata_gene", 
  contains = c("heatdata"),
  representation = representation(),
  prototype = prototype(
    aesn = list(
      x = "sample", y = "gene", fill = "value",
      lab_x = "Samples", lab_y = "Genes", lab_fill = "Gene level"),
    x_aesn = list(x = "sample", y = "group", fill = "group",
      lab_x = "", lab_y = "", lab_fill = "Group"),
    x_pal = color_set(),
    y_aesn = list(x = "module", y = "gene", fill = "module",
      lab_x = "", lab_y = "", lab_fill = "Module"),
    y_pal = wgcna_colors()
    ))

.heatdata_cor <- setClass("heatdata_cor", 
  contains = c("heatdata"),
  representation = representation(),
  prototype = prototype(
    aesh = "color",
    aesn = list(
      x = "col_var", y = "row_var", color = "cor", size = "-log2(P.value)",
      shape = "significant", lab_x = "Columns", lab_y = "Rows",
      lab_color = "Correlation", lab_size = "-log2(P.value)",
      lab_shape = "Significant"),
    para = list(
      clust_row = T, clust_col = T, method = 'average'),
    x_aesn = list(x = "col_var", y = "group", color = "group",
      lab_x = "", lab_y = "", lab_color = "Column Group"),
    x_pal = color_set(),
    y_aesn = list(x = "Group", y = "row_var", color = "group",
      lab_x = "", lab_y = "", lab_color = "Row group"),
    y_pal = wgcna_colors(),
    fun_plot = list(
      xtree = get_fun("plot_xtree"),
      ytree = get_fun("plot_ytree"),
      main = get_fun("dot_heatmap"),
      xgroup = get_fun("add_xgroup.dot.heatmap"),
      ygroup = get_fun("add_ygroup.dot.heatmap"))
    ))

.heatdata_gene_cor <- setClass("heatdata_gene_cor", 
  contains = c("heatdata_cor", "heatdata_gene"),
  representation = representation(),
  prototype = prototype(
    aesn = list(
      x = "trait", y = "module", color = "cor", size = "-log2(P.value)",
      shape = "significant", lab_x = "Traits", lab_y = "",
      lab_color = "Correlation", lab_size = "-log2(P.value)",
      lab_shape = "Significant"),
    x_aesn = list(x = "sample", y = "group", color = "group",
      lab_x = "", lab_y = "", lab_color = "Group"),
    y_aesn = list(x = "Module", y = "module", color = "module",
      lab_x = "", lab_y = "", lab_color = "Module")
    ))

setGeneric("naviRaw", 
  function(x) standardGeneric("naviRaw"))

setMethod("naviRaw", signature = c(x = "heatdata_gene"),
  function(x){
    data_long <- as_data_long(x@raw)
    x@data_long <- tibble::as_tibble(data_long)
    x
  })

setGeneric("new_heatdata", 
  function(x, y, ...) standardGeneric("new_heatdata"))

setMethod("new_heatdata", signature = c(x = "EList"),
  function(x){
    x <- heatdata_gene(raw = x)
    naviRaw(x)
  })

setGeneric("standby", 
  function(x) standardGeneric("standby"))

setMethod("standby", signature = c(x = "heatdata"),
  function(x){
    cols <- unlist(x@aesn[ names(x@aesn) %in% c("x", "y", x@aesh) ])
    .check_columns(x@data_long, cols, "x@data_long")
    main <- dplyr::select(x@data_long, dplyr::all_of(unname(cols)))
    main <- tidyr::spread(main, cols[[ "x" ]], cols[[ x@aesh ]])
    main <- data.frame(main, check.names = F)
    rownames(main) <- main[[ cols[[ "y" ]] ]]
    main <- dplyr::select(main, -!!rlang::sym(cols[[ "y" ]]))
    x@main <- main
    x
  })

setGeneric("set_xmeta", 
  function(x, metadata) standardGeneric("set_xmeta"))

setMethod("set_xmeta", signature = c(x = "heatdata_gene"),
  function(x){
    x@xmeta <- dplyr::distinct(x@data_long, sample, group)
    x
  })

setGeneric("callheatmap", 
  function(x, y, ...) standardGeneric("callheatmap"))

setGeneric("corheatmap", 
  function(x, y, ...) standardGeneric("corheatmap"))

setMethod("corheatmap", signature = c(x = "df", y = "df"),
  function(x, y, row_var = "row_var", col_var = "col_var"){
    callheatmap(new_heatdata(cal_corp(x, y, row_var, col_var)))
  })

setMethod("callheatmap", signature = c(x = "heatdata"),
  function(x){
    y.reformat <- x.reformat <- 0L
    if (!is.null(x@ymeta)) {
      x@data_long <- dplyr::filter(
        x@data_long, 
        !!rlang::sym(x@aesn[[ "y" ]]) %in% x@ymeta[[ x@aesn[[ "y" ]] ]]
      )
      y.reformat <- 1L
    }
    if (!is.null(x@xmeta)) {
      x@data_long <- dplyr::filter(
        x@data_long, 
        !!rlang::sym(x@aesn[[ "x" ]]) %in% x@xmeta[[ x@aesn[[ "x" ]] ]]
      )
      x.reformat <- 1L
    }
    if (any(c(y.reformat, x.reformat)) | is.null(x@main)) {
      x <- standby(x)
    }
    x@gg_main <- do.call(x@fun_plot[[ "main" ]], c(list(x@data_long), x@aesn))
    if (y.reformat) {
      if (any(colnames(x@ymeta) == x@y_aesn[[ x@aesh ]])) {
        args <- c(list(data = x@ymeta, p = NULL, pal = x@y_pal), x@y_aesn)
        x@gg_ygroup <- do.call(x@fun_plot[[ "ygroup" ]], args)
      }
    }
    x@gg_xtree <- x@fun_plot[[ "xtree" ]](x@main, x@para[[ "method" ]])
    x@gg_ytree <- x@fun_plot[[ "ytree" ]](x@main, x@para[[ "method" ]])
    if (x.reformat) {
      if (any(colnames(x@xmeta) == x@x_aesn[[ x@aesh ]])) {
        args <- c(list(data = x@xmeta, p = NULL, pal = x@x_pal), x@x_aesn)
        x@gg_xgroup <- do.call(x@fun_plot[[ "xgroup" ]], args)
      }
    }
    return(x)
  })

draw_sampletree <- function(x) {
  aplot::insert_top(x@gg_xgroup, x@gg_xtree, height = 8) 
}

draw_genetree <- function(x) {
  aplot::insert_left(x@gg_ygroup, x@gg_ytree, width = 8)
}

setMethod("draw", signature = c(x = "heatdata"),
  function(x){
    p <- x@gg_main
    if (is.null(x@gg_xgroup)) {
      p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    }
    if (!is.null(x@gg_ygroup)) {
      p <- aplot::insert_left(p, x@gg_ygroup, width = 0.02) 
    }
    if (x@para[[ "clust_col" ]] & !is.null(x@gg_xtree)) {
      p <- aplot::insert_top(p, x@gg_xtree, height = 0.2)
    }
    if (x@para[[ "clust_row" ]] & !is.null(x@gg_ytree)) {
      fun <- function(n) length(unique(x@data_long[[n]]))
      width <- if (fun(2) / fun(1) > 4) .1 else .2
      p <- aplot::insert_left(p, x@gg_ytree, width = width)
    }
    if (!is.null(x@gg_xgroup)) {
      p <- aplot::insert_bottom(p, x@gg_xgroup, height = 0.05) 
    }
    p
  })

as_df.distframe <- function(data, threshold = NULL) {
  row <- 0
  names <- rownames(data)
  data <- apply(data, 1,
    function(value){
      row <<- row + 1
      data.frame(from = names[row], to = names(value), value = unname(value))
    }, simplify = F)
  data <- do.call(rbind, data)
  data <- tibble::as_tibble(data)
  if (!is.null(threshold)) {
    data <- dplyr::filter(data, value > threshold)
  }
  data
}

as_edges.distframe <- function(data, threshold = NULL) 
{
  data <- as_df.distframe(data, threshold)
  meta <- apply(data[, 1:2], 1, sort)
  meta <- data.frame(t(meta))
  colnames(meta) <- c("from", "to")
  data <- cbind(meta, data[, -(1:2)])
  data <- dplyr::distinct(data, from, to, .keep_all = T)
  data <- dplyr::filter(data, from != to)
  tibble::as_tibble(data)
}

# ==========================================================================
# wgcna
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.wgcData <- setClass("wgcData", 
  contains = c("data.frame"),
  representation = representation(),
    prototype = NULL)

setMethod("show", signature = c(object = "wgcData"),
  function(object){
    print(tibble::as_tibble(data.frame(object)))
    message(crayon::green("Rownames (Sample names):"))
    textSh(crayon::cyan(paste0(rownames(object), collapse = ", ")),
      pre_trunc = T, pre_wrap = T)
  })

setGeneric("as_wgcData", 
  function(x) standardGeneric("as_wgcData"))

setMethod("as_wgcData", signature = c(x = "EList"),
  function(x){
    data <- data.frame(t(x$E))
    colnames(data) <- x$genes[[ 1 ]]
    .wgcData(data)
  })

.wgcTrait <- setClass("wgcTrait", 
  contains = c("wgcData"),
  representation = representation(),
    prototype = NULL)

setGeneric("as_wgcTrait", 
  function(x) standardGeneric("as_wgcTrait"))

setMethod("as_wgcTrait", signature = c(x = "elist"),
  function(x){
    data <- dplyr::select_if(x$targets, is.numeric)
    data <- data.frame(data)
    rownames(data) <- x$targets$sample
    .wgcTrait(data)
  })

cut_tree <- function(tree, height, size) {
  clust <- e(WGCNA::cutreeStatic(tree, height, size))
  clust == 1
}

setGeneric("exclude", 
  function(x, y) standardGeneric("exclude"))

setMethod("exclude", signature = c(x = "wgcData"),
  function(x, y){
    .wgcData(x[y, ])
  })

setGeneric("draw_sampletree", 
  function(x) standardGeneric("draw_sampletree"))

setMethod("draw_sampletree", signature = c(x = "wgcData"),
  function(x){
    x <- hclust(dist(x), method = "average")
    plot(x, main = "Sample clustering", sub = "", xlab = "", cex = .7)
    return(x)
  })

cal_sft <- function(data, powers = c(c(1:10), seq(12, 20, by = 2))) 
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == F")
  }
  sft <- e(WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5))
  sft
}

plot_sft <- function(sft) 
{
  p1 <- ggplot(sft$fitIndices, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
    geom_line(color = "darkred", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2") +
    theme_classic()
  p2 <- ggplot(sft$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_line(color = "darkgreen", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Mean Connectivity") +
    theme_classic()
  require(patchwork)
  p1 + p2
}

wrap <- function(data, width = 10, height = 8) {
  if (is(data, "wrap")) {
    data@width <- width
    data@height <- height
    data
  } else {
    .wrap(data = data, width = width, height = height)
  }
}

.wrap <- setClass("wrap", 
  contains = c("can_not_be_draw"),
  representation = representation(data = "ANY", width = "ANY", height = "ANY"),
  prototype = NULL)

setGeneric("zoom", 
  function(x, s.w, s.h, ...) standardGeneric("zoom"))

setMethod("zoom", signature = c(x = "wrap"),
  function(x, s.w, s.h){
    x@width <- x@width * s.w
    x@height <- x@height * s.h
    return(x)
  })

z7 <- function(x, s.w = .7, s.h = .7) {
  zoom(x, s.w, s.h)
}

setMethod("show", 
  signature = c(object = "wrap"),
  function(object){
    setdev(width = object@width, height = object@height)
    if (is(object@data, "grob.obj")) {
      grid.draw(object@data)
    } else {
      print(object@data)
    }
  })

.wgcNet <- setClass("wgcNet", 
  contains = c("can_not_be_draw", "list"),
  representation = representation(),
  prototype = NULL)

.wgcEigen <- setClass("wgcEigen", 
  contains = c("wgcData"),
  representation = representation(colors = "df", members = "list"),
    prototype = NULL)

.corp <- setClass("corp", 
  contains = c("data_long"),
  representation = representation(),
  prototype = NULL)

setGeneric("cal_corp", 
  function(x, y, ...) standardGeneric("cal_corp"))

setMethod("cal_corp", signature = c(x = "df", y = "df"),
  function(x, y, row_var = "row_var", col_var = "col_var", trans = F){
    x <- data.frame(x)
    y <- data.frame(y)
    if (is.character(x[[1]])) {
      rownames(x) <- x[[1]]
      x <- x[, -1]
      message("Set rolnames of x as the first columns.")
    }
    if (is.character(y[[1]])) {
      rownames(y) <- y[[1]]
      y <- y[, -1]
      message("Set rolnames of y as the first columns.")
    }
    if (trans) {
      x <- t(x)
      y <- t(y)
    }
    cor <- agricolae::correlation(x, y)
    data <- as_data_long(cor$correlation, cor$pvalue, row_var, col_var, "cor", "pvalue")
    .corp(add_anno(.corp(data)))
  })

setMethod("new_heatdata", signature = c(x = "df"),
  function(x){
    if (is.character(x[[1]])) {
      x <- data.frame(x)
      rownames(x) <- x[[1]]
      x <- x[, -1]
    }
    x <- as_data_long(x)
    new_heatdata(x)
  })

setMethod("new_heatdata", signature = c(x = "data_long"),
  function(x){
    object <- .heatdata()
    object@data_long <- tibble::as_tibble(x)
    object@aesn$x <- object@x_aesn$x <- colnames(x)[2]
    object@aesn$y <- object@y_aesn$y <- colnames(x)[1]
    object@aesn$fill <- colnames(x)[3]
    object@aesn$lab_x <- Hmisc::capitalize(colnames(x)[2])
    object@aesn$lab_y <- Hmisc::capitalize(colnames(x)[1])
    object
  })

setMethod("new_heatdata", signature = c(x = "corp"),
  function(x){
    object <- .heatdata_cor()
    object@data_long <- tibble::as_tibble(x)
    object@aesn$x <- object@x_aesn$x <- colnames(x)[2]
    object@aesn$y <- object@y_aesn$y <- colnames(x)[1]
    object@aesn$lab_x <- Hmisc::capitalize(colnames(x)[2])
    object@aesn$lab_y <- Hmisc::capitalize(colnames(x)[1])
    object
  })

setMethod("as_data_long", signature = c(x = "df"),
  function(x, row_var = "rname", col_var = "cname", x_value = "value"){
    x <- dplyr::mutate(tibble::as_tibble(x), !!!nl(row_var, list(rownames(x))))
    x <- tidyr::gather(x, !!col_var, !!x_value, -!!rlang::sym(row_var))
    .data_long(x)
  })

setMethod("as_data_long", signature = c(x = "df", y = "df"),
  function(x, y, row_var = "rname", col_var = "cname", 
    x_value = "x_value", y_value = "y_value")
  {
    x <- dplyr::mutate(tibble::as_tibble(x), !!!nl(row_var, list(rownames(x))))
    x <- tidyr::gather(x, !!col_var, !!x_value, -!!rlang::sym(row_var))
    y <- tidyr::gather(tibble::as_tibble(y), !!col_var, !!y_value)
    x <- dplyr::mutate(x, !!!nl(y_value, list(y[[ y_value ]])))
    .data_long(x)
  })

setGeneric("add_anno", 
  function(x) standardGeneric("add_anno"))

setMethod("add_anno", signature = c(x = "corp"),
  function(x){
    min <- min(x$pvalue[x$pvalue != 0])
    dplyr::mutate(tibble::as_tibble(data.frame(x)),
      `-log2(P.value)` = -log2(ifelse(pvalue == 0, min / 10, pvalue)),
      significant = ifelse(pvalue > .05, "> 0.05",
        ifelse(pvalue > .001, "< 0.05", "< 0.001")),
      sign = ifelse(pvalue > .05, "-",
        ifelse(pvalue > .001, "*", "**"))
    )
  })

setMethod("new_heatdata", signature = c(x = "wgcEigen", y = "wgcTrait"),
  function(x, y){
    cor <- e(WGCNA::cor(x, y, use = "p"))
    pvalue <- e(WGCNA::corPvalueStudent(cor, nrow(x)))
    data <- as_data_long(cor, pvalue, "module", "trait", "cor", "pvalue")
    data <- add_anno(.corp(data))
    .heatdata_gene_cor(
      data_long = data,
      ymeta = dplyr::select(x@colors, module),
      y_pal = nl(x@colors$module, x@colors$color, F)
    )
  })

setMethod("draw", signature = c(x = "heatdata_gene_cor"),
  function(x){
    if (!is.null(x@gg_xgroup)) {
      x@gg_xgroup <- x@gg_xgroup +
        guides(color = "none")
    }
    if (!is.null(x@gg_ygroup)) {
      x@gg_ygroup <- x@gg_ygroup +
        guides(color = "none")
    }
    callNextMethod(x)
  })

get_eigens <- function(net) {
  eigens <- e(WGCNA::orderMEs(net$MEs))
  colors <- e(WGCNA::labels2colors(colorIndex <- unique(net$colors)))
  color_data <- tibble::tibble(module = paste0("ME", colorIndex), color = colors)
  members <- split(names(net$colors), net$colors)
  names(members) <- paste0("ME", names(members))
  .wgcEigen(eigens, colors = color_data, members = members)
}

setMethod("show", signature = c(object = "wgcNet"),
  function(object){
    setdev(width = 12, height = 9)
    mergedColors = e(WGCNA::labels2colors(object$colors))
    e(WGCNA::plotDendroAndColors(object$dendrograms[[1]], mergedColors[object$blockGenes[[1]]],
        "Module colors", dendroLabels = FALSE, hang = 0.01,
        addGuide = TRUE, guideHang = 0.05))
  })

cal_module <- function(data, power, cut_hight = .25, min_size = 30, save_tom = "tom",
  maxBlockSize = 25000, ...)
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == F")
  }
  require(WGCNA)
  net <- e(WGCNA::blockwiseModules(
      data, power = power,
      TOMType = "unsigned", minModuleSize = min_size,
      reassignThreshold = 0, mergeCutHeight = cut_hight,
      numericLabels = TRUE, pamRespectsDendro = FALSE, loadTOM = T,
      saveTOMs = TRUE, saveTOMFileBase = save_tom,
      maxBlockSize = maxBlockSize, verbose = 3, ...
      ))
  .wgcNet(net)
}

setGeneric("clip_data", 
  function(x, by) standardGeneric("clip_data"))

setMethod("clip_data", signature = c(x = "elist", by = "wgcData"),
  function(x, by){
    ## filter sample in Counts data
    validObject(x)
    fun <- function(data, by) {
      data <- data[, colnames(data) %in% rownames(by)]
      data
    }
    x$E <- fun(x$E, by)
    ## filter sample in metadata
    fun <- function(data, by) {
      data <- data[data[[ 1 ]] %in% rownames(by), ]
      data
    }
    x$targets <- fun(x$targets, by)
    ## filter genes in genes and counts data
    fun <- function(data, by) {
      data[[1]] %in% colnames(by)
    }
    logi <- fun(x$genes, by)
    x$genes <- x$genes[logi, ]
    x$E <- x$E[logi, ]
    validObject(x)
    x
  })

# ==========================================================================
# Fast display the content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

send_eval <- function(to,
  subject = "需要评估绩效的业务",
  content = "Hello, 慧姐\n\n这是这个月需要评估绩效的业务，已附在表格中。\n\nBest wish!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0("~/outline/lixiao/", "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "need_eval.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_registers <- function(to,
  subject = "业务表格更新",
  content = "Hello, 慧姐\n\n这是每月末需提交的更新的业务登记表。\n\nBest wish!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0("~/outline/lixiao/", "summary"),
  atts = paste0(path_summary, "/", "生信组表格登记-黄礼闯.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_prin <- function(to,
  subject = "月底交个人行为准则考核表",
  content = "Hello, 慧姐\n\n个人行为准则考核表已交。\n\nBest wish!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0("~/outline/lixiao/", "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "2023年行为准则考核表.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_summary <- function(to,
  subject = "月底交付考核",
  content = "Hello, 慧姐\n\n绩效表已提交。\n\nBest wish!",
  time = Sys.time(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path_summary = paste0("~/outline/lixiao/", "summary_", paste0(year, "-", month)),
  atts = paste0(path_summary, "/", "assess_绩效+软性考核表.xlsx"))
{
  send_that(to, subject, content, atts)
}

send_that <- function(to, subject, content, atts = NULL, from = "huanglichuang@wie-biotech.com") {
  isthat <- usethis::ui_yeah("Are you sure sending the mail?")
  if (!isthat) {
    return(message("Sending cancelled."))
  }
  if (is.character(content)) {
    file <- tempfile("mail_", fileext = ".txt")
    writeLines(content, file)
  }
  if (!is.null(atts)) {
    code_atts <- paste0(" -a ", paste0(atts, collapse = " "), " -- ")
  } else {
    code_atts <- ""
  }
  cdRun("mutt ", to,
    " -e 'my_hdr From:", from, "'",
    " -s '", subject, "'",
    code_atts,
    " < ", file
  )
  message("Done")
}

update_registers <- function(orders = get_orders(),
  path = "~/outline/lixiao/summary",
  name = "黄礼闯",
  register = paste0(path, "/", "生信组表格登记-", name, ".xlsx"),
  templ_dir = "~/outline/lixiao/performance_templ/")
{
  target <- paste0(templ_dir, "/", "生信组表格登记-XX.xlsx")
  if (!dir.exists(path)) {
    dir.create(path)
  }
  if (!file.exists(register)) {
    file.copy(target, register)
  }
  wb <- openxlsx2::wb_load(register)
  pos <- list(2, 1)
  ## base
  data <- dplyr::filter(orders, type != "备单业务")
  data <- dplyr::mutate(data, seq = 1:nrow(data), note = "")
  data <- dplyr::select(data,
    date, seq, info, id, score, member, receive_date, status, title, note
  )
  wb <- openxlsx2::wb_add_data(wb, 1, data, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos))
  ## extra
  data <- dplyr::filter(orders, type == "备单业务")
  data <- dplyr::mutate(data, seq = 1:nrow(data), note = "")
  data <- dplyr::select(data,
    date, seq, info, id, score, member, receive_date, status, title, note
  )
  wb <- openxlsx2::wb_add_data(wb, 2, data, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos))
  openxlsx2::wb_save(wb, register)
  browseURL(normalizePath(register))
}

summary_month <- function(
  time = Sys.time(),
  orders = get_orders(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  path = "~/outline/lixiao",
  templ_dir = "~/outline/lixiao/performance_templ/",
  rm = F)
{
  dir <- paste0(path, "/", "summary_", paste0(year, "-", month))
  targets <- c(ass = "assess_绩效+软性考核表.xlsx", prin = "2023年行为准则考核表.xlsx")
  if (rm) {
    unlink(list.files(dir, full.names = T, all.files = T, recursive = T), T, T)
    file.copy(paste0(templ_dir, "/", targets), dir)
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
    file.copy(paste0(templ_dir, "/", targets), dir)
  }
  targets[] <- paste0(dir, "/", targets)
  wb <- openxlsx2::wb_load(targets[[ "ass" ]])
  ## prepare orders (this month)
  orders <- dplyr::filter(orders, lubridate::year(belong) == !!year,
    lubridate::month(belong) == !!month)
  ## modify the assess table
  pos.date <- list(3, 1)
  date <- paste0(year, "年 ", month, "月 01日", "至",
    year, "年 ", month, "月 ", lubridate::days_in_month(time), "日")
  wb <- openxlsx2::wb_add_data(wb, 1, date,
    dims = do.call(openxlsx2::wb_dims, pos.date))
  if (T) {
    pos.data_ass <- list(7, 1)
    data_ass <- dplyr::mutate(orders,
      seq = 1:nrow(orders), num = 1, title.en = "", note = "", coef = round(coef, 3))
    data_ass <- dplyr::select(data_ass,
      member, seq, id, type, score, num, title, title.en, status, note, coef
    )
    fun <- function(wb, data) {
      openxlsx2::wb_add_data(wb, 1, data, col_names = F,
        dims = do.call(openxlsx2::wb_dims, pos.data_ass), na.strings = "")
    }
    if (any(is.na(data_ass$coef))) {
      wb_eval <- fun(wb, dplyr::filter(data_ass, is.na(coef)))
      openxlsx2::wb_save(wb_eval, paste0(dir, "/need_eval.xlsx"))
    }
    wb <- fun(wb, data_ass)
  }
  if (T) {
    pos.data_sum <- list(29, 6)
    data_ass <- dplyr::mutate(data_ass,
      type = ifelse(type %in% c("固定业务", "备单业务"), "base", "other"))
    data_sum <- split_lapply_rbind(data_ass, ~type,
      function(data) {
        sum <- data.frame(
          num = nrow(data),
          score = "",
          coef = sum(data$coef)
        )
        dplyr::mutate(sum,
          coef_fomu = paste0(paste0(data$coef, collapse = "+"), "=", coef)
        )
      })
    wb <- openxlsx2::wb_add_data(wb, 1, data_sum, col_names = F,
      dims = do.call(openxlsx2::wb_dims, pos.data_sum))
  }
  if (T) {
    pos.sum <- list(29, 10)
    sum <- sum(data_sum$coef)
    wb <- openxlsx2::wb_add_data(wb, 1, sum,
      dims = do.call(openxlsx2::wb_dims, pos.sum))
  }
  openxlsx2::wb_save(wb, targets[[ "ass" ]])
  ########################################
  ########################################
  wb <- openxlsx2::wb_load(targets[[ "prin" ]])
  info <- data.frame(depart = "部门：技术部",
    x1 = "", x2 = "",
    job = "岗位：生物信息工程师",
    x3 = "", name = "姓名：黄礼闯", x4 = "",
    month = paste0("月度：", year, "年", month, "月"))
  pos.info <- list(2, 1)
  wb <- openxlsx2::wb_add_data(wb, 1, info, col_names = F,
    dims = do.call(openxlsx2::wb_dims, pos.info))
  pos.sign <- list(22, 1)
  sign <- paste0("本人确认：            日 期 ：", year, "年 ", month, "月  日")
  wb <- openxlsx2::wb_add_data(wb, 1, sign,
    dims = do.call(openxlsx2::wb_dims, pos.sign))
  openxlsx2::wb_save(wb, targets[[ "prin" ]])
  ## check
  browseURL(normalizePath(targets[[ "ass" ]]))
}

cf <- function(remuneration, base_wage = 6000) {
  remuneration / base_wage
}

odate <- function(month, year = format(Sys.time(), "%Y")) {
  as.Date(paste0(year, "-", month, "-01"))
}

remu <- function(coef, date, base_wage) {
  use <- vapply(date, FUN.VALUE = double(1),
    function(date) {
      which <- tail(which(as.Date(names(base_wage)) <= as.Date(date)), n = 1)
      base_wage[[which]]
    })
  coef * use
}

get_orders <- function(
  dir = "~/outline/lixiao/", pattern = ".items.rds",
  base_wage = list("2023-07-01" = 3000, "2023-08-01" = 4500, "2023-09-01" = 6000))
{
  files <- list.files(dir, pattern, T, T, T)
  lst <- lapply(files,
    function(file) {
      lst <- readRDS(file)
      maybeMulti <- c("coef", "belong", "id", "type")
      lst.m <- lst[ names(lst) %in% maybeMulti ]
      lst <- lst[ !names(lst) %in% maybeMulti ]
      data <- do.call(tibble::tibble, lst)
      if (length(lst.m$coef) > 1) {
        data <- do.call(rbind, rep(list(data), length(lst.m$coef)))
      }
      data <- do.call(dplyr::mutate, c(list(data), lapply(lst.m, unlist)))
      data <- dplyr::mutate(data,
        coef = as.double(coef),
        date = as.Date(date),
        receive_date = as.Date(receive_date)
      )
      data
    })
  data <- as_tibble(data.table::rbindlist(lst, fill = T))
  data <- dplyr::arrange(data, belong, receive_date)
  data <- dplyr::mutate_if(data, is.character,
    function(x) ifelse(is.na(x), "", x))
  data <- dplyr::mutate(data, remuneration = remu(coef, belong, base_wage))
  data <- dplyr::relocate(data, belong, id, title, remuneration)
  data
}

od_get_title <- function() {
  getOption("title")
}

items <- function(
  type = "固定业务",
  title = od_get_title(),
  status = "完成",
  coef = .25,
  date = "2023-07-12",
  info = od_get_info(),
  id = od_get_id(),
  receive_date = od_get_date(),
  score = od_get_score(),
  member = "黄礼闯",
  save = ".items.rds",
  belong = as.Date(receive_date))
{
  if (is.null(id)) {
    stop("The `id` can not be a NULL !!!")
  }
  if (identical(id, "")) {
    stop("The `id` can not be a empty character !!!")
  }
  items <- as.list(environment())
  saveRDS(items, save)
  items
}

od_get_id <- function(...) {
  it <- od_get(..., key = "id")
}

od_get_score <- function(...) {
  it <- od_get(..., key = "score")
}

od_get_info <- function(...) {
  it <- od_get(..., key = "info")
}

od_get_date <- function(file = "./mailparsed/date.md") {
  line <- readLines(file, n = 1)
  as.Date(line, "%A, %d %B %Y")
}

odb <- function(...) {
  n <- 0L
  res <- lapply(list(...),
    function(key) {
      res <- odk(key, if (!n) T else F)
      if (is.null(res)) {
        stop("No keywords of ", key, " found.")
      }
      res
    })
  paste0(unlist(res), collapse = "")
}

odk <- function(key, fresh = F, file = "./mailparsed/part_1.md")
{
  kds <- getOption("od_keywords")
  if (is.null(kds) | fresh) {
    lines <- readLines(file)
    pattern <- paste0("[a-z]+\\[\\[.*?\\]\\]")
    res <- stringr::str_extract_all(paste0(lines, collapse = " "), pattern)
    res <- unlist(res)
    res <- lapply(res,
      function(x) {
        name <- gs(x, "^([a-z]+).*", "\\1")
        value <- gs(x, ".*\\[\\[(.*?)\\]\\].*", "\\1")
        list(name, value)
      })
    kds <- lapply(res, function(x) x[[2]])
    names(kds) <- lapply(res, function(x) x[[1]])
    options(od_keywords = kds)
  }
  kds[[ key ]]
}

od_get <- function(file = "./mailparsed/part_1.md", key = "id",
  pattern = paste0("(?<=", key, "\\{\\{).*?(?=\\}\\})"))
{
  if (file.exists(file)) {
    lines <- readLines(file)
    if (is.null(pattern)) {
      lines
    } else {
      res <- stringr::str_extract(paste0(lines, collapse = " "), pattern)
      if (is.na(res)) {
        ""
      } else {
        res
      }
    }
  } else {
    return()
  }
}

deparse_mail <- function(dir = "mail",
  savedir = "mailparsed", attsdir = "order_material",
  force = F)
{
  if (!dir.exists(dir)) {
    message("No mail directory found.")
    return()
  }
  if (!dir.exists(savedir)) {
    dir.create(savedir)
  }
  if (!force & file.exists(sig <- paste0(savedir, "/", ".parsed.md"))) {
    return()
  }
  testfile <- list.files(dir, full.names = T, recursive = T)[1]
  fewLines <- readLines(testfile, n = 100)
  if (!any(grpl(fewLines, "Content-Transfer-Encoding: base64"))) {
    message("Maybe this mail is not the raw file.")
    return()
  }
  ## import python package
  bt <- e(reticulate::import_builtins())
  m <- e(reticulate::import("mailbox"))
  ## load mail
  mdir <- m$Maildir(dir)
  obj <- mdir$get_message(mdir$keys()[1]) 
  ## payload
  contents <- obj$get_payload()
  isMulti <- vapply(contents, function(x) x$is_multipart(), logical(1))
  ## all attachments
  atts <- contents[ !isMulti ]
  if (!dir.exists(attsdir)) {
    dir.create(attsdir)
  }
  if (length(atts) == 0) {
    writeLines("", paste0(attsdir, "empty.txt"))
  } else {
    n <- 0L
    lapply(atts,
      function(att) {
        bins <- att$get_payload(decode = T)
        filename <- att$get_filename()
        if (length(filename) == 0) {
          return()
        }
        if (grpl(filename, "^=\\?UTF-8")) {
          path <- paste0(savedir, "/", filename)
        } else {
          n <<- n + 1L
          suffix <- stringr::str_extract(filename, "\\.[a-zA-Z0-9]*$")
          if (is.na(suffix)) {
            suffix <- ""
          }
          path <- paste0(attsdir, "/", "file_", n, suffix)
        }
        fp <- bt$open(path, "wb")
        fp$write(bins)
        fp$close()
      })
  }
  if (!length(list.files(attsdir)))
    writeLines("", paste0(attsdir, "empty.txt"))
  ## multipart
  main <- contents[ isMulti ]
  if (length(main) == 0) {
    # return()
    main <- contents
  } else {
    main <- main[[1]]
    main <- main$get_payload()
  }
  n <- 0L
  files_multipart <- lapply(main,
    function(part) {
      n <<- n + 1L
      if (part$get_content_type() == "multipart/alternative") {
        part <- part$get_payload()[[1]]
      }
      bin <- part$get_payload(decode = T)
      if (part$get_content_type() == "text/plain") {
        filename <- paste0("part_", n, ".md")
      } else if (part$get_content_type() == "text/html") {
        filename <- paste0("part_", n, ".html")
      } else {
        return()
      }
      fp <- bt$open(file <- paste0(savedir, "/", filename), "wb")
      fp$write(bin)
      fp$close()
      return(file)
    })
  files_multipart <- lst_clear0(files_multipart)
  writeLines("", sig)
  ## basic information
  date <- lapply(obj$items(),
    function(lst) {
      if (lst[[ 1 ]] == "Date") {
        return(lst[[ 2 ]])
      }
    })
  date <- unlist(date)
  writeLines(date, paste0(savedir, "/", "date.md"))
}

auto_material <- function(class = "job_geo", envir = .GlobalEnv) {
  names <- ls(envir = envir, all.names = all.names)
  info <- lapply(names,
    function(name) {
      obj <- get(name, envir = envir)
      if (is(obj, class)) {
        list(gse = object(obj),
          design = obj@params$about[[1]]@experimentData@other$overall_design)
      } else NULL
    })
  info <- lst_clear0(info)
  info <- lapply(info,
    function(x) {
      c(paste0("- **", x$gse, "**: ", stringr::str_trunc(x$design, 200)), "")
    }
  )
  info <- unlist(info)
  info <- c("All used GEO expression data and their design: ", "", info)
  writeLines(info)
}

auto_method <- function(class = "job", envir = .GlobalEnv) {
  names <- ls(envir = envir, all.names = all.names)
  info <- lapply(names,
    function(name) {
      obj <- get(name, envir = envir)
      if (is(obj, class)) {
        res <- try(list(cite = obj@cite, method = obj@method), T)
        if (inherits(res, "try-error")) {
          NULL
        } else {
          res
        }
      } else NULL
    })
  info <- lst_clear0(info)
  methods <- lapply(info,
    function(lst) {
      if (!is.null(lst$method)) {
        if (!identical(lst$method, character(0))) {
          meth <- paste("-", lst$method)
          if (!is.null(lst$cite)) {
            meth <- paste(meth, lst$cite)
          }
          paste0(meth, ".")
        }
      }
    })
  methods <- unique(unlist(methods))
  methods <- c("Mainly used method:", "", methods,
  "- Other R packages (eg., `dplyr` and `ggplot2`) used for statistic analysis or data visualization.")
  writeLines(methods)
}

set_cover <- function(title, author = "LiChuang Huang", date = Sys.Date(),
  coverpage = "../cover_page.pdf", institution = "@立效研究院")
{
  options(title = title)
  content <- strwrap(paste0("\\begin{titlepage}
      \\newgeometry{top=7.5cm}
      \\ThisCenterWallPaper{1.12}{", coverpage, "}
      \\begin{center}
      \\textbf{\\Huge ", title, "}
      \\vspace{4em}
      \\begin{textblock}{10}(3,5.9)
      \\huge \\textbf{\\textcolor{white}{", date, "}}
      \\end{textblock}
      \\begin{textblock}{10}(3,7.3)
      \\Large \\textcolor{black}{", author, "}
      \\end{textblock}
      \\begin{textblock}{10}(3,11.3)
      \\Large \\textcolor{black}{", institution, "}
      \\end{textblock}
      \\end{center}
      \\end{titlepage}
      \\restoregeometry
      "
      ), 50)
  writeLines(content)
}

set_index <- function() {
  if (knitr::is_latex_output()) {
    cat("\\pagenumbering{roman}\n\n")
    cat("\\tableofcontents\n\n")
    cat("\\listoffigures\n\n")
    cat("\\listoftables\n\n")
    cat("\\newpage\n\n")
    cat("\\pagenumbering{arabic}\n\n")
  }
}

autor_preset <- function(...) {
  knitr::opts_chunk$set(
    echo = F, eval = F, message = F, warning = F,
    fig.cap = character(0),
    out.width = "\\linewidth", ...)
  fun_fig.cap <- function(options) {
    options$fig.cap <- Hmisc::capitalize(gsub("-", " ", options$label))
    options
  }
  knitr::opts_hooks$set(fig.cap = fun_fig.cap)
}

# ==========================================================================
# autor: save object, summarise object, show object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## orinal function, for save file, and return file name
autosv <- function(x, name, ...) {
  if (!exists("autoRegisters")) {
    autoRegisters <- character(0)
  }
  if (!any(name == names(autoRegisters))) {
    if (!is(x, "files"))
      file <- select_savefun(x)(x, name, ...)
    else if (file.exists(x)) {
      file <- as.character(x)
      if (is(x, "file_fig")) {
        file.copy(file, dir <- get_savedir("figs"), T)
        file <- paste0(dir, "/", get_filename(file))
      }
    } else {
      stop("file.exists(x) == F")
    }
    autoRegisters <<- c(autoRegisters, nl(name, file, F))
  } else {
    file <- autoRegisters[[ name ]]
    .message_info("autor", "file.exists:", name, sig = "[INFO]")
  }
  return(file)
}

setGeneric("autor", 
  function(x, name, ...) standardGeneric("autor"))

## autor get chunk label for name.
## Be carefull, if the method not passed into this sub method,
## the reference of table or figure may be error, unless the `name`
## is identical with the chunk label.
## **Note**
## the `name` decided the auto reference, savename, and the values
## in `autoRegisters`.
## However, the 'fig.cap' was decided by chunk label, unless
## the `autor_preset` was not performed.
## The best way to use `autor`, don't manualy set the `name`.
## As there was no way to modify the chunk label by `name` in
## codes.
setMethod("autor", signature = c(x = "ANY", name = "missing"),
  function(x, ...){
    name <- knitr::opts_current$get("label")
    if (!is.null(name))
      autor(x, name, ...)
    else {
      message("Not in knitr circumstance, show object only.")
      show(x)
    }
  })

setMethod("autor", signature = c(x = "list", name = "character"),
  function(x, name, ...){
    file <- autosv(x, name, ...)
    autor(file, name)
  })

setMethod("autor", signature = c(x = "can_not_be_draw", name = "character"),
  function(x, name, ...){
    file <- autosv(x, name, ...)
    autor(file, name, ...)
    if (!is.null(lich <- attr(x, "lich"))) {
      abstract(lich, name)
      file <- autosv(lich, name <- paste0(name, "-content"))
      locate_file(name, "上述信息框内容已保存至")
    }
  })

setClassUnion("can_be_draw", c("gg.obj", "heatdata", "grob.obj"))
## autor for ggplot
setMethod("autor", signature = c(x = "can_be_draw", name = "character"),
  function(x, name, ...){
    file <- autosv(x, name, ...)
    autor(file, name, ...)
  })

## autor for data.frame
setMethod("autor", signature = c(x = "df", name = "character"),
  function(x, name, ..., asis = T){
    x <- dplyr::select_if(x,
      function(x) is.character(x) | is.numeric(x) | is.logical(x) | is.factor(x))
    file <- autosv(x, name, ...)
    if (asis)
      abstract(x, name = name, ...)
    include(x, name, ...)
  })

## autor for figures of file
setMethod("autor", signature = c(x = "fig", name = "character"),
  function(x, name, ..., asis = T){
    file <- autosv(x, name, ...)
    if (asis)
      abstract(x, name, ...)
    include(x, name, ...)
  })

setMethod("autor", signature = c(x = "files", name = "character"),
  function(x, name, ..., asis = T){
    file <- autosv(x, name, ...)
    if (asis)
      abstract(x, name, ...)
  })

setMethod("autor", signature = c(x = "lich", name = "character"),
  function(x, name, ...){
    abstract(x, name, ...)
  })

setMethod("autor", signature = c(x = "character", name = "character"),
  function(x, name, ...){
    if (length(x) > 1)
      stop("length(x) == 1")
    if (!file.exists(x))
      stop("file.exists(x) == F")
    fig.type <- c(".jpg", ".png", ".pdf")
    file.type <- stringr::str_extract(x, "\\.[a-zA-Z]+$")
    if (!is.na(file.type)) {
      if (any(file.type == fig.type))
        return(autor(fig(x), name, ...))
    }
    autor(files(x), name, ...)
  })

# ==========================================================================
# show object in report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setGeneric("include", 
  function(x, name, ...) standardGeneric("include"))

## include for fig
setMethod("include", signature = c(x = "fig"),
  function(x, name, ...){
    if (knitr::is_latex_output()) {
      cat("\n\\def\\@captype{figure}\n")
      cat("\\begin{center}\n",
        "\\includegraphics[width = 0.9\\linewidth]{", as.character(x), "}\n",
        "\\caption{", knitr::opts_current$get("fig.cap"),
        "}\\label{fig:", name, "}\n",
        "\\end{center}\n", sep = "")
    } else {
      inclu.fig(as.character(x), saveDir = "report_picture")
    }
  })

setMethod("include", signature = c(x = "df"),
  function(x, name, ...){
    x <- tibble::as_tibble(x)
    if (knitr::is_latex_output()) {
      x <- trunc_table(x)
      knitr::kable(x, "markdown", caption = as_caption(name))
    } else {
      print(x)
    }
  })

as_caption <- function(str) {
  Hmisc::capitalize(gsub("-", " ", str))
}

trunc_table <- function(x) {
  x <- dplyr::mutate_all(x, as.character)
  x <- dplyr::mutate_all(x, function(str) stringr::str_trunc(str, 8))
  colnames(x) %<>% stringr::str_trunc(8)
  if (nrow(x) > 15) {
    x <- head(x, n = 15)
    blank <- head(x, n = 0)
    blank[1, ] <- "..."
    x <- dplyr::bind_rows(x, blank)
  }
  if (ncol(x) > 10) {
    x <- x[, 1:10]
  }
  col <- vapply(colnames(x), nchar, integer(1))
  col <- which(cumsum(col) > 80)
  if (length(col) > 0) {
    col <- col[1]
    x <- x[, 1:col]
    x$... <- "..."
  }
  x
}

asis <- function(object) {
  knitr::asis_output(object)
}

# ==========================================================================
# select save function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setGeneric("select_savefun", 
  function(x, ...) standardGeneric("select_savefun"))

setMethod("select_savefun", signature = c(x = "list"),
  function(x){
    fun <- function(x, class) {
      vapply(x, function(obj) is(obj, class), logical(1))
    }
    if (all(fun(x, "easywrite"))) {
      function(lst, name)
        get_fun("writeDatas")(lst, get_realname(name))
    } else if (all(fun(x, "gg.obj"))) {
      function(lst, name, ...)
        get_fun("writePlots")(lst, get_realname(name))
    } else if (all(fun(x, "wrap"))) {
      function(lst, name, ...)
        get_fun("writeWraps")(lst, get_realname(name))
    } else {
      stop("None function found for save")
    }
  })

setMethod("select_savefun", signature = c(x = "can_not_be_draw"),
  function(x){
    get_fun("write_graphics")
  })

## select_savefun for ggplot
setMethod("select_savefun", signature = c(x = "gg.obj"),
  function(x){
    get_fun("write_gg")
  })

setMethod("select_savefun", signature = c(x = "heatdata"),
  function(x){
    function(x, ...) {
      write_gg(draw(x), ...)
    }
  })

setMethod("select_savefun", signature = c(x = "character"),
  function(x){
    get_fun("write_character")
  })

setMethod("select_savefun", signature = c(x = "grob.obj"),
  function(x){
    get_fun("write_grob")
  })

setMethod("select_savefun", signature = c(x = "df"),
  function(x){
    if (!is(x, "data.frame")) {
      x <- tibble::as_tibble(x)
    }
    data <- dplyr::select_if(x, is.character)
    check <- apply(data, 2,
      function(ch) {
        any(grepl(",", ch))
      })
    if (any(check)) {
      ## around 4.8 Mb
      if (object.size(x) < 5e6)
        get_fun("write_xlsx2")
      else
        get_fun("write_tsv2")
    } else get_fun("fwrite2")
  })

# ==========================================================================
# summarise the object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

setGeneric("abstract", 
  function(x, name, latex, ...) {
    standardGeneric("abstract")
  })

setMethod("abstract", signature = c(x = "ANY", name = "character", latex = "missing"),
  function(x, name, ...){
    restype <- knitr::opts_current$get("results")
    if (is.null(restype)) {
      warning("is.null(knitr::opts_current$get(\"results\")) == T")
    } else if (restype != "asis") {
      warning("restype != \"asis\"")
    }
    if (knitr::is_latex_output())
      latex <- T
    else
      latex <- NULL
    reCallMethod("abstract", namel(x, name, latex), ...)
  })

## abstract for table
setMethod("abstract", signature = c(x = "df", name = "character", latex = "logical"),
  function(x, name, latex, ..., key = 1, abs = NULL, summary = T, sum.ex = NULL){
    x <- tibble::as_tibble(x)
    cat("Table \\@ref(tab:", name, ")", " (下方表格) ",
      "为表格", gsub("-", " ", name), "概览。\n", sep = "")
    if (!is.null(abs))
      cat(abs, "\n")
    locate_file(name)
    if (!is.null(summary)) {
      cat(text_roundrect(fix.tex(sumTbl(x, key, sum.ex))))
    }
  })

## abstract for figure of file
setMethod("abstract", signature = c(x = "fig", name = "character", latex = "logical"),
  function(x, name, latex, ..., abs = NULL){
    cat("Figure \\@ref(fig:", name, ")", " (下方图) ",
      "为图", gsub("-", " ", name), "概览。\n", sep = "")
    if (!is.null(abs))
      cat(abs, "\n")
    locate_file(name)
  })

setMethod("abstract", signature = c(x = "lich", name = "character", latex = "logical"),
  function(x, name, latex, ..., abs = NULL){
    str <- sapply(names(x),
      function(name){
        ch <- c("\n\\textbf{", name, ":}\n\n\\vspace{0.5em}\n")
        ch <- c(ch, strwrap(stringr::str_trunc(x[[ name ]], 1500), indent = 4, width = 60))
        ch <- c(ch, "\n\\vspace{2em}\n")
        ch
      })
    cat(text_roundrect(paste0(fix.tex(unlist(str)), collapse = "\n")))
  })

setMethod("abstract", signature = c(x = "files", name = "character", latex = "logical"),
  function(x, name, latex, ..., abs = NULL, sum.ex = NULL){
    if (dir.exists(x)) {
      cat(abs, "\n")
      cat("`", as_caption(name), "' 数据已全部提供。", "\n", sep = "")
      locate_file(name)
      cat(text_roundrect(sumDir(autoRegisters[[ name ]], sum.ex)))
    } else {
      cat(abs, "\n")
      cat("`", as_caption(name), "' 数据已提供。", "\n", sep = "")
      locate_file(name)
    }
  })

sumDir <- function(dir, sum.ex = NULL) {
  files <- list.files(dir)
  num <- length(files)
  if (length(files) > 5)
    files <- c(head(files, n = 5), "...")
  files <- fix.tex(files)
  paste0("注：文件夹", fix.tex(dir), "共包含", num, "个文件。\n",
    sum.ex, "\n",
    paste0(
      "\\begin{enumerate}",
      "\\tightlist\n",
      paste0(paste0("\\item ", files), collapse = "\n"),
      "\n\\end{enumerate}",
      collapse = "\n")
  )
}

fix.tex <- function(str) {
  gsub("_", "\\\\_", str)
}

sumTbl <- function(x, key, sum.ex) {
  paste0("注：表格共有", nrow(x), "行", ncol(x), "列，",
    "以下预览的表格可能省略部分数据；",
    "表格含有", colSum(x[[ key ]]),
    "个唯一`", colnames(x[, key]), "'。\n", sum.ex)
}

locate_file <- function(name, des = "对应文件为") {
  if (!exists('autoRegisters'))
    stop("!exists('autoRegisters')")
  if (!file.exists(autoRegisters[[ name ]]))
    stop("file.exists(autoRegisters[[ name ]] == F)")
  cat("\n**(", des, " `", autoRegisters[[ name ]], "`)**", "\n", sep = "")
}

text_roundrect <- function(str, collapse = "\n") {
  paste0("\\begin{center}",
    "\\begin{tcolorbox}[colback=gray!10, colframe=gray!50, width=0.9\\linewidth, arc=1mm, boxrule=0.5pt]",
    str, "\\end{tcolorbox}\n\\end{center}", collapse = collapse
  )
}

# ==========================================================================
# wrapper for dplyr tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# mutate <- dplyr::mutate
# filter <- dplyr::filter
# arrange <- dplyr::arrange
# distinct <- dplyr::distinct
# select <- dplyr::select
# rename <- dplyr::rename
# relocate <- dplyr::relocate
# slice <- dplyr::slice
# slice_max <- dplyr::slice_max
# slice_min <- dplyr::slice_min
# group_by <- dplyr::group_by

lapply(c("mutate", "filter", "arrange", "distinct",
    "select", "rename", "relocate", "slice", "slice_max",
    "slice_min", "group_by"),
  function(name) {
    setGeneric(name, function(DF_object, ...) DF_object)
    setMethod(name, signature = c(DF_object = "df"),
      function(DF_object, ..., fun_name = name){
        fun <- get_fun(fun_name, asNamespace("dplyr"))
        if (!is(DF_object, "tbl_df")) {
          DF_object <- tibble::as_tibble(DF_object)
        }
        object <- fun(DF_object, ...)
        object
      })
  })

setGeneric("as_tibble", 
  function(x) standardGeneric("as_tibble"))

setMethod("as_tibble", signature = c(x = "df"),
  function(x){
    rownames <- rownames(x)
    x <- tibble::as_tibble(x)
    if (!identical(rownames, as.character(1:nrow(x)))) {
      x <- dplyr::mutate(x, rownames = !!rownames)
      x <- dplyr::relocate(x, rownames)
    }
    return(x)
  })

# ==========================================================================
# for fast combine object
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

xf <- function(..., DATA_object = list(), Env = parent.frame(1)) {
  DATA_object <- get_from_env(weight <- list(...), DATA_object, Env)
  frame_col(weight, DATA_object)
}

yf <- function(..., DATA_object = list(), Env = parent.frame(1)) {
  DATA_object <- get_from_env(weight <- list(...), DATA_object, Env)
  frame_row(weight, DATA_object)
}

get_from_env <- function (weight, data = list(), env = parent.frame(1)){
  names <- names(weight)
  sapply(names, simplify = F,
    function(name) {
      if (is.null(data[[ name ]]))
        get(name, envir = env)
      else
        data[[ name ]]
    })
}

# ==========================================================================
# ROC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

plot_roc <- function(roc) {
  plot(1- x$specificities, x$sensitivities)
}

# ==========================================================================
# upset plot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.upset <- setClass("upset_data", 
  contains = c("can_not_be_draw", "tbl_df"),
  representation = representation(params = "list"),
  prototype = NULL)

setMethod("show", signature = c(object = "upset_data"),
  function(object){
    data <- suppressWarnings(data.frame(as_tibble(object)))
    if (!is.null(object@params$trunc))
      colnames(data) %<>% stringr::str_trunc(object@params$width, object@params$trunc)
    data <- dplyr::select(data, -members)
    maxnum <- max(apply(data, 2, sum))
    upset <- UpSetR::upset(data, sets = colnames(data), nintersects = NA,
      sets.bar.color = "lightblue", order.by = "freq",
      set_size.show = T, set_size.scale_max = 1.5 * maxnum,
      text.scale = c(1.4, 1, 1.4, 1.1, 1, 1.1)
    )
    show(wrap(upset, ncol(data) * 1.4, ncol(data) * 1.3))
  })

new_upset2 <- function(..., lst = NULL, trunc = "left", width = 30, convert = T) {
  p <- new_upset(..., lst = lst, trunc = trunc, width = width, convert = convert)
  if (is.null(lst)) {
    ins <- ins(...)
  } else {
    ins <- ins(lst = lst)
  }
  namel(p, ins)
}

new_upset <- function(..., lst = NULL, trunc = "left", width = 30, convert = T) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  lst <- lapply(lst, unique)
  members <- unique(unlist(lst, use.names = F))
  data <- data.frame(members = members)
  lst <- lapply(lst,
    function(set) {
      ifelse(data$members %in% set, 1L, 0L)
    })
  data <- do.call(mutate, c(list(data), lst))
  data <- .upset(data, params = namel(trunc, width))
  if (convert) {
    show(data)
    recordPlot()
  } else {
    data
  }
}

new_venn <- function(..., lst = NULL, wrap = T, fun_pre = rm.no) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  lst <- lapply(lst, fun_pre)
  p <- ggVennDiagram::ggVennDiagram(lst) +
    scale_fill_gradient(low = "grey95", high = sample(color_set(), 1)) +
    rstyle("theme") +
    theme(axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()) +
    geom_blank()
  if (wrap) {
    p <- wrap(p, 4, 2.5)
  }
  attr(p, "ins") <- ins <- intersect(lst[[1]], lst[[2]])
  attr(p, "lich") <- new_lich(list(Intersection = ins))
  p
}

setdev <- function(width, height) {
  name <- names(dev.cur())
  if (name == "null device")
    dev.new(width = width, height = height)
}

new_allu <- function(data, col.fill = 1, axes = 1:2, label.auto = F, label.freq = NULL, label.factor = 1) {
  require(ggalluvial)
  fill <- colnames(data)[col.fill]
  data <- dplyr::mutate(data, fill = !!rlang::sym(fill))
  data <- to_lodes_form(data, key = "Types", axes = axes)
  if (label.auto) {
    freq <- table(data$stratum)
    if (is.null(label.freq)) {
      label.notshow <- names(freq)[ as.integer(freq) <= fivenum(as.integer(freq))[4] * label.factor ]
    } else {
      label.notshow <- names(freq)[ as.integer(freq) <= label.freq ]
    }
    fun <- function(x) {
      ifelse(x %in% label.notshow, "", as.character(x))
    }
    data <- dplyr::mutate(data, label = fun(stratum))
  } else {
    data <- dplyr::mutate(data, label = stratum)
  }
  aes <- aes(x = Types, y = 1, label = label, stratum = stratum, alluvium = alluvium)
  p.alluvial <- ggplot(data, aes) +
    geom_alluvium(aes(fill = fill)) +
    geom_stratum(fill = "lightyellow") +
    geom_text(stat = "stratum") +
    labs(fill = "", y = "") +
    scale_fill_manual(values = color_set()) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
      axis.title = element_blank(),
      legend.position = "none") +
    geom_blank()
  p.alluvial
}

new_col <- function(..., lst = NULL, fun = function(x) x[ !is.na(x) & x != ""]) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  lst <- vapply(lst, function(x) length(fun(unique(x))), double(1))
  data <- data.frame(var = names(lst), value = unname(lst))
  p <- ggplot(data, aes(x = reorder(var, value), y = value, fill = value)) +
    geom_col(width = .5) +
    geom_text(aes(x = var, y = value + max(value) * .01, label = value), hjust = 0, size = 3) +
    ylim(c(0, max(data$value) * 1.2)) +
    coord_flip() +
    labs(x = "", y = "") +
    theme(legend.position = "")
  wrap(p, 7, nrow(data) * .5 + .5)
}

new_pie <- function(x, title = NULL, use.ggplot = T, fun_text = ggplot2::geom_text) {
  x <- split(x, x)
  x <- vapply(x, length, integer(1))
  if (use.ggplot) {
    data <- data.frame(var = names(x), value = unname(x))
    data <- dplyr::arrange(data, dplyr::desc(var))
    data <- dplyr::mutate(data,
      label = paste0("(", round(value / sum(value) * 100, 1), "%)"),
      label = paste0(var, " ", label), lab.x = .2,
      lab.y = value / 2 + c(0, cumsum(value)[ -length(value) ])
    )
    palette <- if (length(x) > 10) color_set() 
      else ggsci::pal_npg()(10)
    p <- ggplot(data, aes(x = 0L, y = value, fill = var)) +
      geom_bar(stat = 'identity', position = 'stack', width = 1) +
      fun_text(aes(x = lab.x, y = lab.y, label = label)) +
      scale_fill_manual(values = palette) +
      labs(x = '', y = '', title = '') +
      coord_polar(theta = 'y') +
      theme_minimal() +
      theme(legend.position = "none",
        axis.text = element_blank(),
        plot.margin = margin(-.1, -.1, -.1, -.1, "npc"),
        panel.grid = element_blank()) +
      geom_blank()
    wrap(as_grob(p), 5, 4)
  } else {
    grob <- ggplotify::base2grob(expression({
      par(mar = rep(1, 4))
      pie(x, main = title)
    }))
    wrap(grob, 5, 4)
  }
}

plot_median_expr_line <- function(data) {
  lst <- pbapply::pblapply(split(as_tibble(data), ~ .id),
    function(data) {
      data <- lapply(split(data, ~ sample),
        function(x) {
          fn <- fivenum(x$value)
          names(fn) <- n(v, 5)
          do.call(data.frame, as.list(fn))
        })
      data <- data.table::rbindlist(data, idcol = T)
      rename(data, sample = .id)
    })
  seq <- 1:nrow(lst[[1]])
  e(spiralize::spiral_initialize(range(seq)))
  e(spiralize::spiral_track(range(lst[[1]]$v3)))
  e(spiralize::spiral_lines(seq, max(lst[[1]]$v3), type = "h", gp = gpar(col = "grey70")))
  spiralize::spiral_lines(seq, lst[[1]]$v3, gp = gpar(col = "red"))
  spiralize::spiral_lines(seq, lst[[2]]$v3, gp = gpar(col = "blue"))
  grid.text("Median Expression line", y = .97, gp = gpar(cex = 2))
  lgd <- e(ComplexHeatmap::packLegend(
    ComplexHeatmap::Legend(title = "From", type = "lines", legend_gp = gpar(col = c("blue", "red"), lwd = 2),
      at = c("Raw", "Normalized"))
  ))
  ComplexHeatmap::draw(lgd)
  p <- recordPlot()
  dev.off()
  attr(p, "data") <- lst
  p
}

fix.html.str <- function(x) {
  gs(x, " ", "%20")
}

rm.no <- function(x) {
  unique(x[ !is.na(x) & x != "" ])
}

get_fe_data <- function(use.symbol = T, for_gsea = F,
  path = "../ferroptosis_2023-10-24.rds")
{
  if (F) {
    fe_db <- list(marker = "~/Downloads/ferroptosis_marker.csv",
      driver = "~/Downloads/ferroptosis_driver.csv",
      suppressor = "~/Downloads/ferroptosis_suppressor.csv",
      unclassifier = "~/Downloads/ferroptosis_unclassified.csv",
      inducer = "~/Downloads/ferroptosis_inducer.csv",
      inhibitor = "~/Downloads/ferroptosis_inhibitor.csv",
      disease = "~/Downloads/ferroptosis_disease.csv"
    )
    fe_db <- lapply(fe_db, ftibble)
    saveRDS(fe_db, "../ferroptosis_2023-10-24.rds")
  }
  data <- readRDS(path)
  if (use.symbol) {
    data <- lapply(data,
      function(x) {
        if (any(colnames(x) == "symbol")) {
          return(x)
        }
      })
    data <- lst_clear0(data)
  }
  if (for_gsea) {
    gsea <- data.table::rbindlist(data, idcol = T, fill = T)
    gsea <- dplyr::mutate(gsea, term = paste0("Ferroptosis_", .id))
    gsea <- as_tibble(dplyr::relocate(gsea, term, symbol))
    return(gsea)
  }
  data
}

grp <- function(x, pattern, ...) {
  grep(pattern, x, ...)
}

grpl <- function(x, pattern, ...) {
  grepl(pattern, x, ...)
}

grpf <- function(x, pattern, ...) {
  x[grepl(pattern, x, ...)]
}

## igraph tips
# grid star circle

get_layout <- function(edges = NULL, layout = "grid", nodes = NULL, ...) {
  data <- as_tibble(data.frame(fast_layout(edges, layout, nodes = nodes, ...)))
  data
}
