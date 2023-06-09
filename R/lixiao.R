# ==========================================================================
# work and function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
  if (is.null(sinkFile)) {
    tryCatch(system(paste0(...)), finally = setwd(owd))
  } else {
    tryCatch(
      capture.output(system(paste0(...)),
        file = sinkFile, split = T),
      finally = setwd(owd))
  }
}

vina_limit <- function(lig, recep, timeLimit = 120, ...) {
  try(vina(lig, recep, ..., timeLimit = timeLimit), T)
}

vina <- function(lig, recep, dir = "vina_space",
  exhaustiveness = 32, scoreing = "ad4", stout = "/tmp/res.log", timeLimit = 60)
{
  if (!file.exists(dir)) {
    dir.create(dir, F)
  } 
  subdir <- paste0(reals <- get_realname(c(lig, recep)), collapse = "_into_")
  wd <- paste0(dir, "/", subdir)
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
      " --scoring ", scoreing,
      " --exhaustiveness ", exhaustiveness,
      " --out ", subdir, "_out.pdbqt",
      " >> ", stout), T)
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

nl <- function(names, values, as.list = T, ...) {
  .as_dic(values, names, as.list = as.list, ...)
}

smiles_as_sdfs.obabel <- function(smiles) {
  lst.sdf <- pbapply::pbsapply(smiles, simplify = F,
    function(smi) {
      ChemmineOB::convertFormat("SMI", "SDF", source = test)
    })
  lst.sdf
}

tbmerge <- function(...) {
  tibble::as_tibble(merge(...))
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

get_pdb <- function(ids, cl = 3, mkdir.pdb = "protein_pdb") {
  dir.create(mkdir.pdb, F)
  ids <- tolower(ids)
  ids <- grouping_vec2list(ids, round(length(ids) / cl), T)
  pbapply::pblapply(ids, cl = cl,
    function(ids) {
      tmp <- tempfile(fileext = ".txt")
      cat(ids, sep = ", ", file = tmp)
      system(paste0("get_pdb.sh -f ", tmp, " -p -o ", mkdir.pdb))
    })
  lapply(list.files(mkdir.pdb, ".*\\.gz$", full.names = T), R.utils::gunzip)
  files <- list.files(mkdir.pdb, ".*\\.pdb$", full.names = T)
  files <- nl(gsub(".*?([^/]{1,})\\.pdb$", "\\1", files), files)
  files
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

start_drive <- function(command = "java -jar ~/operation/selenium.jar",
  port = 4444, extra = NULL, ...)
{
  system(paste(command, "-port", port, extra), wait = F)
  new_link(port = port)
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
  checks <- list.files("~/Downloads", "\\.xlsx$", full.names = T)
  if (length(checks) > 0)
    file.remove(checks)
  link$open()
  jobn <- 0
  lapply(urls,
    function(url) {
      link$navigate(url)
      jobn <<- jobn + 1
      if (jobn > 1) {
        link$refresh()
      }
      res <- try(silent = T, {
        ele <- link$findElement(
          using = "xpath",
          value = paste0("//main//div//h4[text()='", heading, "']/../div//button")
        )
        ele$sendKeysToElement(list("Download", key = "enter"))
      })
      if (inherits(res, "try-error")) {
        writeLines("", paste0("~/Downloads/empty_", jobn, ".xlsx"))
      }
    })
}

moveToDir_herbs <- function(ids,
  file.pattern = "\\.xlsx$", index.pfun = file_seq.by_time,
  from = "~/Downloads", to = "herbs_ingredient", suffix = ".xlsx", .id = "herb_id")
{
  args <- as.list(environment())
  files <- do.call(moveToDir, args)
  data <- lapply(files,
    function(file) {
      res <- try(openxlsx::read.xlsx(file), T)
      if (inherits(res, "try-error")) {
        tibble::tibble()
      } else res
    })
  names(data) <- ids
  data <- tibble::as_tibble(data.table::rbindlist(data, idcol = T, fill = T))
  cols <- colnames(data)
  data <- dplyr::relocate(data, .id)
  colnames(data)[ which(cols == ".id") ] <- .id
  data
}

moveToDir <- function(ids,
  file.pattern, index.pfun = file_seq.by_time,
  from, to, suffix)
{
  files <- list.files(from, file.pattern, full.names = T)
  index <- index.pfun(files)
  files <- files[order(index)][1:length(ids)]
  files <- rev(files)
  dir.create(to, F)
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

new_link <- function(port = 4444L, browser = "firefox", addr = "localhost")
{
  RSelenium::remoteDriver(
    remoteServerAddr = addr, port = port,
    browserName = browser
  )
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

writePlots <- function(lst, dir, width = 7, height = 7, postfix = ".pdf") {
  fun <- function(p, file) ggsave(file, p, width = width, height = height)
  writeDatas(lst, dir, fun, postfix = postfix)
}

writeDatas <- function(lst, dir, fun = data.table::fwrite, postfix = ".csv") 
{
  if (is.null(names(lst))) {
    stop("is.null(names(lst)) == T")
  }
  dir.create(dir, F)
  n <- 1
  lapply(names(lst),
    function(name) {
      data <- lst[[name]]
      file <- paste0(dir, "/", n, "_", name, postfix)
      n <<- n + 1
      if (is.null(data)) {
        writeLines("", file)
      } else {
        fun(data, file)
      }
    })
}

new_biomart <- function(dataset = c("hsapiens_gene_ensembl"))
{
  ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = match.arg(dataset))
  ensembl
}

filter_biomart <- function(mart, attrs, filters = "", values = "") {
  anno <- biomaRt::getBM(attributes = attrs, filters = filters,
    values = values, mart = mart)
  tibble::as_tibble(anno)
}

list_attrs <- function(mart) {
  tibble::as_tibble(biomaRt::listAttributes(mart))
}

general_attrs <- function() {
  attrs <- c("ensembl_gene_id",
    "entrezgene_id",
    "hgnc_symbol",
    "pdb",
    "chromosome_name",
    "start_position",
    "end_position",
    "description")
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

get_nci60_data <- function(comAct = "../comAct_nci60/DTP_NCI60_ZSCORE.xlsx",
  rna = "../rna_nci60/RNA__RNA_seq_composite_expression.xls")
{
  comAct <- readxl::read_xlsx(comAct, skip = 8)
  rna <- readxl::read_xls(rna, skip = 10)
  rna <- dplyr::mutate(rna, genes = gsub("\\.[0-9]*$", "", `Gene name d`))
  list(comAct = comAct, rna = rna)
}

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

new_stringdb <- function(
  score_threshold = 200,
  species = 9606,
  network_type = c("physical", "full"),
  input_directory = "../")
{
  STRINGdb::STRINGdb$new(score_threshold = score_threshold,
    species = species, network_type = match.arg(network_type),
    input_directory = input_directory
  )
}

create_interGraph <- function(sdb, data, col = "name", rm.na = T) {
  mapped <- sdb$map(data.frame(data), col, removeUnmappedRows = rm.na)
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

add_attr.igraph <- function(igraph, data, by.x = "name", by.y, dedup.edges = T)
{
  comps <- igraph::as_data_frame(igraph, "both")
  nodes <- merge(comps$vertices, data, by.x = by.x, by.y = by.y, all.x = T)
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
  arr.len = 2, edge.color = 'lightblue', edge.width = 1)
{
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      arrow = arrow(length = unit(arr.len, 'mm')),
      color = edge.color, width = edge.width) +
    # geom_node_point(aes(x = x, y = y), shape = 21, stroke = .3, fill = "grey70") +
    geom_node_label(aes(label = preferred_name), size = label.size) +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    theme_void()
  p
} 

plot_networkFill.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, node.size = 12, sc = 5, ec = 5, 
  arr.len = 2, edge.color = 'lightblue', edge.width = 1, lab.fill = "MCC score")
{
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      # arrow = arrow(length = unit(arr.len, 'mm')),
      color = edge.color, width = edge.width) +
    geom_node_point(aes(x = x, y = y, fill = MCC_score),
      size = node.size, shape = 21, stroke = .3) +
    geom_node_text(aes(label = preferred_name), size = label.size) +
    scale_fill_gradient(low = "lightyellow", high = "red") +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    theme_void() +
    theme(plot.margin = margin(r = .05, unit = "npc")) +
    geom_blank()
  p
} 

multi_enrichKEGG <- function(lst.entrez_id) 
{
  res <- pbapply::pblapply(lst.entrez_id,
    function(ids) {
      res.kegg <- clusterProfiler::enrichKEGG(ids)
      res.path <- tibble::as_tibble(res.kegg@result)
      res.path
    })
  res
}

multi_enrichGO <- function(lst.entrez_id, orgDb = 'org.Hs.eg.db')
{
  res <- pbapply::pblapply(lst.entrez_id,
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
        guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
        coord_flip() +
        ylim(zoRange(data$GeneRatio, 1.3)) +
        theme_minimal() +
        theme(axis.title.x = element_blank())
      p
    })
}

vis_enrich.go <- function(lst, cutoff_p.adjust = .1, maxShow = 10) {
  res <- lapply(lst,
    function(data) {
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
    })
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

fwrite2 <- function(data, file, mkdir = "tabs") {
  if (!file.exists(mkdir))
    dir.create(mkdir)
  data.table::fwrite(data, paste0(mkdir, "/", file))
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
  tops <- head(tops$name, n = 10)
  data <- igraph::as_data_frame(igraph, "both")
  nodes <- dplyr::filter(data$vertices, name %in% !!tops)
  nodes <- merge(nodes, resMcc, by = "name", all.x = T)
  edges <- dplyr::filter(data$edges, (from %in% !!tops) & (to %in% !!tops))
  igraph <- igraph::graph_from_data_frame(edges, F, nodes)
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
  main = search_resFile("index.Rmd"), 
  headPattern = "[^r][^e][^s].*业务.*\\.zip", head = list.files(".", headPattern),
  masterSource = c("output.Rmd", "install.R", "read_me.txt", head),
  prefix = "results_", zip = paste0(prefix, head),
  clientZip = "client.zip", masterZip = "master.zip",
  clear = T)
{
  file.copy(report, "readMe_introduction.pdf", T)
  report <- "readMe_introduction.pdf"
  files.client <- c(report, main)
  files.master <- c(files.client, head, masterSource)
  zip(clientZip, files.client)
  zip(masterZip, files.master)
  zip(zip, c(clientZip, masterZip), flags = "-ur9X")
  if (clear) {
    file.remove(c(clientZip, masterZip))
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

esearch <- function(query = NULL, fetch.save = paste0(gsub(" ", "_", query), ".xml"),
  path = "~/operation", tract.save = "res.tsv",
  fields = c("SortPubDate", "Title", "FullJournalName", "Name", "Id"))
{
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
  colnames(res) <- fields
  if (any("SortPubDate" == fields)) {
    res <- dplyr::mutate(res, SortPubDate = as.Date(SortPubDate))
    res <- dplyr::arrange(res, dplyr::desc(SortPubDate))
  }
  if (any("Id" == fields)) {
    res <- dplyr::mutate(res, Id = as.character(Id))
  }
  res
}

split_lapply_rbind <- function(data, f, fun, ...) {
  data <- split(data, f)
  data <- lapply(data, fun, ...)
  data <- data.table::rbindlist(data)
  tibble::as_tibble(data)
}

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

setMethod("pal", 
  signature = c(x = "andata"),
  function(x){ 
    pal <- x@palette
    if (is.null(pal)) {
      group <- unique(x@metadata$group)
      pal <- nl(group, MCnebula2:::.get_color_set()[1:length(group)], F)
    }
    pal
  })

.andata_pca <- setClass("andata_pca", contains = c("andata"))

pca_data.long <- function(data.long, fun_scale = function(x) scale(x)) {
  lst <- .split_data.long(data.long)
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

setMethod("plot_andata", 
  signature = c(x = "andata_pca"),
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

setMethod("plot_andata", 
  signature = c(x = "andata_opls"),
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
  select = rlang::quos(title, ),
  pattern = c("diagnosis", "Sex", "^age", "^time point", "data_processing"),
  abbrev = c("data_processing"))
{
  res <- lapply(lst,
    function(eset){
      tibble::as_tibble(eset@phenoData@data)
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
      MCnebula2:::textSh(lst[[ name ]], wrap_width = width, pre_wrap = T)
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

setMethod("show", 
  signature = c(object = "lich"),
  function(object){
    show_lst.ch(object)
  })

