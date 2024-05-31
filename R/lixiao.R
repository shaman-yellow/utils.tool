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

setMethod("$", signature = c(x = "fig"),
  function(x, name){
    attr(x, name)
  })

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

get_c2_data <- function(pattern = NULL,
  mode = c("hsa", "mmu"))
{
  mode <- match.arg(mode)
  path <- if (mode == "hsa") {
    .prefix("human_c2_v5p2.rdata", "db")
  } else if (mode == "mmu"){
    .prefix("mouse_c2_v5p2.rdata", "db")
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
  job <- .job(method = "Database of `MSigDB` (c2, curated gene sets) was used for signiture screening", cite = "")
  .add_internal_job(job)
  return(db)
}

writeWraps <- function(lst, dir, width = 7, height = 7, ..., postfix = ".pdf") 
{
  fun <- function(p, file) write_graphics(p, file, mkdir = get_path(file))
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

frbind <- function(lst, ...) {
  dplyr::as_tibble(data.table::rbindlist(lst, ...))
}

ftibble <- function(files, ...) {
  fun <- function(file, ...) {
    data <- data.table::fread(file, ...)
    if (any(dup <- duplicated(colnames(data)))) {
      message("Duplicated colnames found, adding suffix.")
      thedup <- colnames(data)[dup]
      colnames(data)[dup] <- paste0(thedup, ".dup.", 1:length(thedup))
    }
    tibble::as_tibble(data)
  }
  if (length(files) > 1 || is(files, "list")) {
    lapply(files, fun, ...)
  } else {
    fun(files, ...)
  }
}

fxlsx <- function(file, ...) {
  as_tibble(data.frame(openxlsx::read.xlsx(file, ...)))
}

fxlsx2 <- function(file, n = NULL, .id = "sheet", bind = T, pattern = NULL, ...) {
  sheets <- openxlsx::getSheetNames(file)
  if (is.null(n)) {
    if (is.null(pattern)) {
      n <- 1:length(sheets)
    } else {
      n <- grp(sheets, pattern)
    }
  }
  lst <- lapply(n,
    function(n) {
      openxlsx::read.xlsx(file, sheet = n, ...)
    })
  names(lst) <- sheets[n]
  if (bind) {
    data <- data.table::rbindlist(lst, idcol = .id, fill = T)
    as_tibble(data)
  } else {
    lapply(lst, as_tibble)
  }
}

get_nci60_data <- function(comAct = .prefix("comAct_nci60/DTP_NCI60_ZSCORE.xlsx", "db"),
  rna = .prefix("rna_nci60/RNA__RNA_seq_composite_expression.xls", "db"))
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

set_showtext <- function() {
  options(SHOWTEXT = T)
}

write_graphics <- function(data, name, ..., file = paste0(get_realname(name), ".pdf"), page = -1,
  mkdir = get_savedir("figs"))
{
  if (!file.exists(mkdir))
    dir.create(mkdir)
  file <- paste0(mkdir, "/", file)
  if (is(data, "wrap")) {
    if (!is.null(fam <- getOption("font_family"))) {
      pdf(file, width = data@width, height = data@height, family = fam)
    } else {
      pdf(file, width = data@width, height = data@height)
    }
  } else {
    pdf(file)
  }
  showtext <- getOption("SHOWTEXT", F)
  if (showtext) {
    showtext::showtext_begin()
  }
  show(data)
  if (showtext) {
    showtext::showtext_end()
  }
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

split_lapply_rbind <- function(data, f, fun, ..., verbose = F, args = list(use.names = T)) {
  data <- split(data, f)
  if (verbose)
    data <- pbapply::pblapply(data, fun, ...)
  else
    data <- lapply(data, fun, ...)
  data <- do.call(data.table::rbindlist, c(list(data), args))
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

show_lst.ch <- function(lst, width = 60) {
  sapply(names(lst),
    function(name){
      message("+++ ", name, " +++\n")
      textSh(lst[[ name ]], wrap_width = width, pre_wrap = T)
    })
  message()
}

.lich <- setClass("lich", 
  contains = c("list"),
  representation = 
    representation(),
  prototype = NULL)

new_lich <- function(lst) {
  if (!is(lst, "list")) {
    message("Not 'list', try format as 'list'.")
    lst <- list(Content = lst)
  }
  if (length(lst) == 1) {
    if (is.null(names(lst))) {
      names(lst) <- "Content"
    }
  }
  if (is.null(names(lst))) {
    stop("The `lst` without names.")
  }
  lst <- lapply(lst, paste0, collapse = ", ")
  .lich(lst)
}

new_hp.cor <- function(data, ..., sig = T, fontsize = 6, names = c("From", "To")) {
  if (length(unique(data[[ names[2] ]])) == 1) {
    cluster_columns <- F
  } else {
    cluster_columns <- T
  }
  if (length(unique(data[[ names[1] ]])) == 1) {
    cluster_rows <- F
  } else {
    cluster_rows <- T
  }
  groups <- dplyr::groups(data)
  if (length(groups)) {
    data <- dplyr::ungroup(data)
    groups <- lapply(groups,
      function(x) {
        if (length(unique(data[[ x ]])) > 1) x
      })
    groups <- lst_clear0(groups)
    if (length(groups)) {
      data <- dplyr::group_by(data, !!!rlang::syms(groups))
    }
  }
  p.hp <- tidyHeatmap::heatmap(data, !!rlang::sym(names[1]), !!rlang::sym(names[2]), cor,
    cluster_columns = cluster_columns, cluster_rows = cluster_rows,
    column_names_gp = gpar(fontsize = fontsize),
    row_names_gp = gpar(fontsize = fontsize),
    ...
  )
  if (sig) {
    p.hp <- tidyHeatmap::layer_star(p.hp, pvalue < .05)
  }
  palette <- fun_color()
  p.hp@input$col <- palette
  if (!identical(p.hp@input$col, palette)) {
    stop("!identical(p.hp@input$col, palette)")
  }
  x.num <- length(unique(data[[ names[1] ]]))
  y.num <- length(unique(data[[ names[2] ]]))
  p.hp <- wrap(p.hp,
    if (y.num > 100) 18 else 1.5 + .12 * y.num,
    if (x.num > 50) 12 else 1.5 + .12 * x.num)
  p.hp <- .set_lab(p.hp, "correlation heatmap")
  p.hp
}

show_col.thp <- function(x) {
  scales::show_col(attr(x, "colors"))
}

setMethod("show", signature = c(object = "lich"),
  function(object){
    show_lst.ch(object)
  })

# ==========================================================================
# limma and edgeR
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# #' @importClassesFrom edgeR DGEList
# #' @importClassesFrom limma EList
.eset <- c("DGEList", "EList")
setFakeClasses(.eset)
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

mx <- function(...){
  design <- model.matrix(...)
  colnames(design) %<>% gsub(".*?group\\.?", "", .)
  colnames(design) %<>% gsub(".*?batch\\.?", "batch.", .)
  design
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
  str <- unique(str)
  str <- tolower(make.names(str[ !is.na(str) ]))
  gs(gs(str, "_|^x\\.", "."), "[.]+", ".")
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
  function(x, HLs = NULL){
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
    x@gg_main <- do.call(x@fun_plot[[ "main" ]], c(list(x@data_long), x@aesn, list(HLs = HLs)))
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

setMethod("[[", signature = c(x = "wrap"),
  function(x, i, ...){
    attr(x, i)
  })

setMethod("[[<-", signature = c(x = "wrap"),
  function(x, i, ..., value){
    attr(x, i) <- value
    return(x)
  })

setMethod("$", signature = c(x = "wrap"),
  function(x, name){
    x[[ name ]]
  })

setMethod("$<-", signature = c(x = "wrap"),
  function(x, name, value){
    x[[ name ]] <- value
    return(x)
  })

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
  function(x, y, row_var = "row_var", col_var = "col_var", trans = F, fast = T)
  {
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
    if (fast) {
      cor <- agricolae::correlation(x, y)
      data <- as_data_long(cor$correlation, cor$pvalue, row_var, col_var, "cor", "pvalue")
      .corp(add_anno(.corp(data)))
    } else {
      alls <- pbapply::pbapply(x, 2, simplify = F,
        function(i) {
          group <- apply(y, 2, simplify = F,
            function(j) {
              fit <- lm(j ~ i)
              tibble::tibble(
                cor = fit$coefficients[[2]],
                pvalue = signif(summary(fit)$coefficients[2, 4], 5),
                model = list(fit$model)
              )
            })
          names(group) <- colnames(y)
          frbind(group, idcol = col_var)
        })
      names(alls) <- colnames(x)
      alls <- frbind(alls, idcol = row_var)
      add_anno(.corp(alls))
    }
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
  function(x, ...) standardGeneric("add_anno"))

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

.di <- setClass("di", 
  contains = c("list"),
  representation = representation(),
  prototype = NULL)

d <- function(ref, mode = c("all", "chinese", "english")) {
  dic <- getOption("dic")
  if (is.null(dic)) {
    stop("No `dic` in options found. Please use function `dic` to set that.")
  }
  mode <- match.arg(mode)
  names(dic) <- tolower(names(dic))
  ref <- match.arg(ref, names(dic))
  res <- dic[[ ref ]]
  if (mode == "all") {
    paste0(res$ch, " (", res$en, ", ", res$abs, ") ")
  } else if (mode == "chinese") {
    res$ch
  } else if (mode == "english") {
    res$en
  }
}

dic <- function(...) {
  dic <- lapply(list(...),
    function(x) {
      if (is(x, "di")) {
        x
      } else if (is.character(x)) {
        if (grpl(x, "^[a-zA-Z]")) {
          di(en = x)
        } else {
          di(ch = x)
        }
      }
    })
  names(dic) <- lapply(dic, function(x) x$abs)
  options(dic = dic)
  mess <- vapply(dic, FUN.VALUE = character(1),
    function(x) {
      paste0("# ", x$abs, ": ", x$en, " ", x$ch)
    })
  writeLines(mess)
}

di <- function(ch = NULL, en = NULL, abs = NULL) {
  if (is.null(ch) & is.null(en)) {
    stop("`ch` and `en` can not be both 'NULL'")
  }
  if (is.null(ch)) {
    ch <- trans.google(en, to = "zh-CN")
  } else if (is.null(en)) {
    en <- trans.google(ch, to = "en")
  }
  if (is.null(abs)) {
    abs <- paste0(stringr::str_extract_all(en, "(?<= |^).")[[1]], collapse = "")
    abs <- toupper(abs)
  }
  .di(as.list(environment()))
}

trans <- function(str) {
  fun <- function(str) {
    if (grpl(str, "^[A-Za-z]")) {
      trans.google(str, to = "zh-CN")
    } else {
      trans.google(str, to = "en")
    }
  }
  if (length(str) > 1) {
    ref <- unique(str)
    ref <- sapply(ref, simplify = F,
      function(x) {
        fun(x)
    })
    print(ref)
    unlist(ref[ match(str, names(ref)) ])
  } else {
    fun(str)
  }
}

trans.google <- function(str, from = "auto", to = "zh-CN") {
  e(gtranslate::translate(str, from = from, to = to))
}

trans.youdao <- function(str) {
  res <- try(e(ecce::translate(str)), silent = T)
  if (inherits(res, "try-error")) {
    stop("The `str` (", str, ") could not be translated (error).")
  }
  res
}

.set_API_youdao <- function() {
  if (Sys.getenv("app_key") == "") {
    Sys.setenv(
      app_key = "5206fcce5ba590d5",
      app_secret = "gORc2vFSztfXn1JmJHA5VNHdgsrzzBmi"
    )
  }
}

workflow_publish <- function(file = "index.Rmd", output = "output.Rmd", title = "",
  fun = write_articlePdf)
{
  browseURL(path <- fun(file, output, title))
  file.copy(path, paste0(getOption("title"), ".pdf"), T)
}

order_publish <- function(file = "index.Rmd", output = "output.Rmd", title = "",
  fun = write_articlePdf)
{
  browseURL(fun(file, output, title))
}

auto_material <- function(class = "job_PUBLISH", envir = .GlobalEnv) {
  names <- .get_job_list(envir)
  info <- lapply(names,
    function(name) {
      obj <- .obtain_job(name, envir, class)
      if (is(obj, "job_geo")) {
        if (obj@step >= 1) {
          x <- list(gse = object(obj),
            design = obj@params$about[[1]]@experimentData@other$overall_design)
          list(type = "geo",
            content = c(paste0("- **", x$gse, "**: ", stringr::str_trunc(x$design, 200)), ""))
        }
      } else if (is(obj, "job_publish")) {
        if (length(obj@cite)) {
          x <- list(cite = obj@cite, method = obj@method)
          list(type = "publish",
            content = c(paste0("- ", x$method, " ", gs(x$cite, "\\[@(.*)\\]", "\\1"), x$cite, "."))
          )
        }
      }
    })
  info <- lst_clear0(info)
  showThat <- function(name, des) {
    info <- lapply(info, function(x) if (x$type == name) x$content)
    info <- unique(unlist(info))
    if (length(info)) {
      info <- c(des, "", info)
      writeLines(info)
    }
  }
  showThat("geo", "All used GEO expression data and their design:")
  showThat("publish", "Other data obtained from published article (e.g., supplementary tables):")
}

auto_method <- function(rm = NULL, class = "job", envir = .GlobalEnv, exclude = "job_publish")
{
  names <- .get_job_list(envir)
  info <- lapply(names,
    function(name) {
      obj <- .obtain_job(name, envir, class)
      if (!is.null(rm)) {
        if (any(class(obj) == rm)) {
          return(NULL)
        }
      }
      if (is(obj, exclude)) {
        NULL
      } else if (is(obj, class)) {
        if (any(obj@cite == rm)) {
          NULL
        } else {
          res <- try(list(cite = obj@cite, method = obj@method), T)
          if (inherits(res, "try-error")) {
            NULL
          } else {
            res
          }
        }
      } else NULL
    })
  info <- lst_clear0(info)
  methods <- lapply(info,
    function(lst) {
      if (!is.null(lst$method)) {
        if (!identical(lst$method, character(0))) {
          meth <- paste("-", lst$method)
          if (!identical(lst$cite, character(0))) {
            meth <- paste0(meth, lst$cite)
          }
          paste0(meth, ".")
        }
      }
    })
  methods <- unique(unlist(methods))
  methods <- c("Mainly used method:", "", methods,
    paste0("- ",
      R.version$version.string, "; ",
      "Other R packages (eg., `dplyr` and `ggplot2`) used for statistic analysis or data visualization.")
  )
  writeLines(methods)
}

.get_job_list <- function(
  envir = .GlobalEnv, extra = getOption("internal_job"),
  type = c("name", "class"))
{
  type <- match.arg(type)
  names <- ls(envir = envir, all.names = F)
  names <- lapply(names,
    function(x) {
      if (is(get(x, envir = envir), "job")) {
        if (type == "name") {
          x
        } else {
          class(get(x, envir = envir))
        }
      } else NULL
    })
  names <- unlist(names)
  if (!is.null(extra)) {
    if (type == "name") {
      if (identical(class(extra), "list")) {
        names <- c(as.list(names), extra)
      } else {
        names <- c(as.list(names), list(extra))
      }
    } else {
      if (!identical(class(extra), "list")) {
        extra <- list(extra)
      }
      names <- c(as.list(names), unlist(lapply(extra, class)))
      return(unique(unlist(names)))
    }
  }
  names
}

get_job_source <- function(jobs = .get_job_list(type = "class"), path = "~/utils.tool/") {
  files <- list.files(path, "workflow.*\\.R$", full.names = T, recursive = T)
  names <- get_filename(files)
  jobs <- gs(jobs, "^job_", "")
  isThat <- vapply(names, FUN.VALUE = logical(1),
    function(name) {
      any(grpl(name, paste0(paste0("_", jobs, "\\."), collapse = "|")))
    })
  files[ isThat ]
}

.obtain_job <- function(name, envir, class = "job") {
  if (!is(name, class) & is.character(name)) {
    obj <- get(name, envir = envir)
  } else {
    obj <- name
  }
  obj
}

.add_internal_job <- function(job, clear = F, limit = 20) {
  if (clear) {
    job@params <- list()
    job@plots <- list()
    job@tables <- list()
    job@object <- NULL
  }
  size <- as.double(gs(obj.size(job), "[a-zA-Z]", ""))
  if (size > limit) {
    warning("Too large `job` (", size, ") add into 'internal_job' options (limit: ", limit, ").")
  }
  injobs <- getOption("internal_job", list())
  if (!is(job, "job_publish")) {
    hasMethods <- unlist(lapply(injobs, function(x) if (is(x, "job")) x@method else NULL))
    if (!any(job@method == hasMethods)) {
      injobs <- c(injobs, nl(class(job), list(job)))
      options(internal_job = injobs)
    }
  } else {
    hasCites <- unlist(lapply(injobs, function(x) if (is(x, "job")) x@cite else NULL))
    if (!any(job@cite == hasCites)) {
      injobs <- c(injobs, nl(class(job), list(job)))
      options(internal_job = injobs)
    }
  }
}

set_cover <- function(title, author = "LiChuang Huang", date = Sys.Date(),
  coverpage = .prefix("cover_page.pdf"), institution = "@立效研究院")
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

set_index <- function(fig = T, tab = T) {
  if (knitr::is_latex_output()) {
    cat("\\pagenumbering{roman}\n\n")
    cat("\\tableofcontents\n\n")
    if (fig)
      cat("\\listoffigures\n\n")
    if (tab)
      cat("\\listoftables\n\n")
    cat("\\newpage\n\n")
    cat("\\pagenumbering{arabic}\n\n")
  }
}

autor_preset <- function(echo = F, eval = F, ...) {
  knitr::opts_chunk$set(
    echo = echo, eval = eval, message = F, warning = F,
    fig.cap = character(0), collapse = F,
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

autorm <- function(names) {
  autoRegisters <<- autoRegisters[ -which(names(autoRegisters) %in% names) ]
}

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
        .dir <- get_savedir("figs")
        nfile <- paste0(.dir, "/", name, strx(file, "\\.[^.]+$"))
        if (get_path(file) != .dir) {
          file.copy(file, nfile, T)
        }
        file <- nfile
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
    if (is(x, "rms") || is(x, "validate")) {
      if (knitr::is_latex_output()) {
        options(prType = "latex")
      } else {
        options(prType = "plain")
      }
      set_showtext()
      return(x)
    }
    name <- knitr::opts_current$get("label")
    if (!is.null(name)) {
      autor(x, name, ...)
    } else {
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
  function(x, name, ..., asis = getOption("autor_asis", T)){
    if (knitr::is_latex_output()) {
      cat("\\begin{center}\\vspace{1.5cm}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{89}\\end{center}")
    }
    if (any(vapply(x, class, character(1)) == "list")) {
      x <- dplyr::mutate(x,
        dplyr::across(dplyr::where(is.list),
          function(x) {
            vapply(x, paste0, collapse = " | ", FUN.VALUE = character(1))
          }))
    }
    x <- dplyr::select_if(x,
      function(x) is.character(x) | is.numeric(x) | is.logical(x) | is.factor(x))
    file <- autosv(x, name, ...)
    if (asis) {
      abstract(x, name = name, ...)
      if (!is.null(lich <- attr(x, "lich"))) {
        abstract(lich, name = name)
      }
    }
    include(x, name, ...)
    if (knitr::is_latex_output()) {
      cat("\n\n\\begin{center}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{89}\\vspace{1.5cm}\\end{center}")
    }
  })

## autor for figures of file
setMethod("autor", signature = c(x = "fig", name = "character"),
  function(x, name, ..., asis = getOption("autor_asis", T)){
    if (knitr::is_latex_output()) {
      cat("\\begin{center}\\vspace{1.5cm}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{88}\\end{center}")
    }
    file <- autosv(x, name, ...)
    if (asis) {
      abstract(x, name = name, ...)
    }
    include(x, name, ...)
    if (!is.null(lich <- attr(x, "lich"))) {
      if (asis) {
        abstract(lich, name = name)
      }
    }
    if (knitr::is_latex_output()) {
      cat("\n\n\\begin{center}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{88}\\vspace{1.5cm}\\end{center}")
    }
  })

setMethod("autor", signature = c(x = "files", name = "character"),
  function(x, name, ..., asis = getOption("autor_asis", T)){
    file <- autosv(x, name, ...)
    if (knitr::is_latex_output()) {
      cat("\n\n\\begin{center}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{85}\\vspace{1.5cm}\\end{center}")
    }
    if (asis)
      abstract(x, name, ...)
    if (knitr::is_latex_output()) {
      cat("\n\n\\begin{center}\\pgfornament[anchor=center,ydelta=0pt,width=9cm]{85}\\vspace{1.5cm}\\end{center}")
    }
  })

setMethod("autor", signature = c(x = "lich", name = "character"),
  function(x, name, ...){
    abstract(x, name, ...)
    file <- autosv(x, name <- paste0(name, "-content"))
    locate_file(name, "上述信息框内容已保存至")
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
      print(knitr::kable(x, "markdown", caption = as_caption(name)))
    } else {
      print(x)
    }
  })

as_caption <- function(str) {
  Hmisc::capitalize(gsub("-", " ", str))
}

trunc_table <- function(x) {
  if (ncol(x) > 10) {
    x <- x[, 1:10]
  }
  width <- if (ncol(x) > 7) {
    9
  } else if (ncol(x) > 5) {
    13
  } else if (ncol(x) > 3) {
    20
  } else {
    30
  }
  x <- dplyr::mutate_all(x, as.character)
  x <- dplyr::mutate_all(x, function(str) stringr::str_trunc(str, width))
  colnames(x) %<>% stringr::str_trunc(width)
  if (nrow(x) > 15) {
    x <- head(x, n = 15)
    blank <- head(x, n = 0)
    blank[1, ] <- "..."
    x <- dplyr::bind_rows(x, blank)
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
    if (knitr::is_latex_output()) {
      latex <- T
    } else {
      latex <- NULL
    }
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
      cat(text_roundrect(fix.tex(
            c(sumTbl(x, key, sum.ex),
              .enumerate_items(.get_des(colnames(x)))
            )
            )))
    }
  })

.enumerate_items <- function(ch) {
  if (length(ch)) {
    paste0(
      "\\begin{enumerate}",
      "\\tightlist\n",
      paste0(paste0("\\item ", names(ch), ": ", unname(ch)), collapse = "\n"),
      "\n\\end{enumerate}",
      collapse = "\n")
  }
}

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
    if (length(x) > 5) {
      x <- head(x, n = 5)
      x <- c(x, list("(Others)" = "..."))
    }
    str <- sapply(names(x),
      function(name){
        text <- x[[ name ]]
        text <- gs(text, "\\\\href\\{", "\\\\url{")
        ch <- c("\n\\textbf{", name, ":}\n\n\\vspace{0.5em}\n")
        ch <- c(ch, strwrap(stringr::str_trunc(text, 300), indent = 4, width = 60))
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

sumTbl <- function(x, key, sum.ex = NULL, mustSum = getOption("abstract.mustSum")) {
  if (!is.null(mustSum)) {
    exkey <- which( colnames(x) %in% mustSum )
    key <- unique(sort(c(key, exkey)))
  }
  sums <- paste0("含有", apply(x[, key], 2, function(x) length(unique(x))),
    "个唯一`", colnames(x[, key]), collapse = "；")
  paste0("注：表格共有", nrow(x), "行", ncol(x), "列，",
    "以下预览的表格可能省略部分数据；", sums, "'。\n", sum.ex)
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

setGeneric("mutate", 
  function(x, ...) standardGeneric("mutate"))

setMethod("mutate", signature = c(x = "df"),
  function(x, ...){
    dplyr::mutate(x, ...)
  })

setGeneric("filter", 
  function(x, ...) standardGeneric("filter"))

setMethod("filter", signature = c(x = "df"),
  function(x, ...){
    dplyr::filter(x, ...)
  })

setGeneric("arrange", 
  function(x, ...) standardGeneric("arrange"))

setMethod("arrange", signature = c(x = "df"),
  function(x, ...){
    dplyr::arrange(x, ...)
  })

setGeneric("distinct", 
  function(x, ...) standardGeneric("distinct"))

setMethod("distinct", signature = c(x = "df"),
  function(x, ...){
    dplyr::distinct(x, ...)
  })

setGeneric("select", 
  function(x, ...) standardGeneric("select"))

setMethod("select", signature = c(x = "df"),
  function(x, ...){
    dplyr::select(x, ...)
  })

setGeneric("rename", 
  function(x, ...) standardGeneric("rename"))

setMethod("rename", signature = c(x = "df"),
  function(x, ...){
    dplyr::rename(x, ...)
  })

setGeneric("relocate", 
  function(x, ...) standardGeneric("relocate"))

setMethod("relocate", signature = c(x = "df"),
  function(x, ...){
    dplyr::relocate(x, ...)
  })

setGeneric("slice", 
  function(x, ...) standardGeneric("slice"))

setMethod("slice", signature = c(x = "df"),
  function(x, ...){
    dplyr::slice(x, ...)
  })

# lapply(c("mutate", "filter", "arrange", "distinct",
#     "select", "rename", "relocate", "slice", "slice_max",
#     "slice_min", "group_by"),
#   function(name) {
#     setGeneric(name, function(DF_object, ...) DF_object)
#     setMethod(name, signature = c(DF_object = "df"),
#       function(DF_object, ..., fun_name = name){
#         fun <- get_fun(fun_name, asNamespace("dplyr"))
#         if (!is(DF_object, "tbl_df")) {
#           DF_object <- tibble::as_tibble(DF_object)
#         }
#         object <- fun(DF_object, ...)
#         object
#       })
#   })

setGeneric("as_tibble", 
  function(x, ...) standardGeneric("as_tibble"))

setMethod("as_tibble", signature = c(x = "df"),
  function(x, ...){
    rownames <- rownames(x)
    x <- tibble::as_tibble(x, ...)
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

new_upset <- function(..., lst = NULL, trunc = "left", width = 30, convert = T, ins = NULL) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  raw.lst <- lst <- lapply(lst, unique)
  members <- unique(unlist(lst, use.names = F))
  data <- data.frame(members = members)
  lst <- lapply(lst,
    function(set) {
      ifelse(data$members %in% set, 1L, 0L)
    })
  data <- do.call(dplyr::mutate, c(list(data), lst))
  data <- .upset(as_tibble(data), params = namel(trunc, width))
  if (convert) {
    show(data)
    p <- wrap(recordPlot())
    p$ins <- ins(lst = raw.lst)
    lich <- list(All_intersection = p$ins)
    if (!is.null(ins)) {
      exLich <- list("Other intersection" = paste0(
          vapply(ins, FUN.VALUE = character(1),
            function(x) {
              if (length(x) != 2) {
                stop("Each `ins` should be length 2 of integer.")
              }
              anno <- paste0(paste0("\"", names(raw.lst)[x], "\""), collapse = " WITH ")
              paste0(anno, ": ", paste0(ins(lst = raw.lst[x]), collapse = ", "))
            }),
          collapse = "\n\\newline\n"))
      lich <- c(lich, exLich)
    }
    p$lich <- new_lich(lich)
    p$raw <- raw.lst
    p
  } else {
    data
  }
}

new_venn <- function(..., lst = NULL, wrap = T, fun_pre = rm.no) {
  if (is.null(lst)) {
    lst <- list(...)
  }
  lst <- lapply(lst, function(x) as.character(fun_pre(x)))
  p <- ggVennDiagram::ggVennDiagram(lst, label_percent_digit = 1) +
    scale_fill_gradient(low = "grey95", high = sample(color_set(), 1)) +
    theme_void() +
    theme(axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()) +
    geom_blank()
  if (wrap) {
    p <- wrap(p, 4, 2.5)
  }
  attr(p, "ins") <- ins <- ins(lst = lst)
  attr(p, "lich") <- new_lich(list(All_intersection = ins))
  lab(p) <- paste0("Intersection of ", paste0(names(lst), collapse = " with "))
  p
}

setdev <- function(width, height) {
  name <- names(dev.cur())
  if (name == "null device")
    dev.new(width = width, height = height)
}

## logistic
new_lrm <- function(data, formula, rev.level = F, lang = c("cn", "en"), B = 500, ...)
{
  fun_escape_bug_of_rms <- function() {
    wh <- which(vapply(1:ncol(data), function(n) is.factor(data[[ n ]]), FUN.VALUE = logical(1)))
    if (any(grpl(colnames(data)[ wh ], "\\s"))) {
      stop("`Due to the bug of `rms`, the names of columns of which is factor, can not contains blank")
    }
  }
  fun_escape_bug_of_rms()
  if (is.character(formula)) {
    message("Detected `formula` input with 'character', use as 'y'.")
    if (length(formula) > 1) {
      stop("length(formula) > 1")
    }
    if (!formula %in% colnames(data)) {
      stop("The 'y' not found in colnames of data.")
    }
    formula <- paste0(formula, " ~ ",
      paste0(
        paste0("`", colnames(data)[ colnames(data) != formula ], "`"),
        collapse = " + "
      )
    )
    message("Guess Formula: ", formula)
    formula <- as.formula(formula)
  }
  isColChar <- apply(data, 2, is.character)
  if (any(isColChar)) {
    message("Convert all character columns as factor.")
    data <- dplyr::mutate(data, dplyr::across(dplyr::where(is.character), as.factor))
  }
  lang <- match.arg(lang)
  y <- as.character(formula[[ 2 ]])
  outcome <- data[[ y ]]
  if (rev.level) {
    data[[ y ]] <- factor(outcome, levels = rev(levels(outcome)))
  }
  message("Check levels: ", paste0(levels <- levels(data[[ y ]]), collapse = ", "))
  set_rms_datadist(data)
  fit <- e(rms::lrm(formula, data = data, x = T, y = T))
  # 95\\% CL
  # exp(confint.default(lrm.eff$fit))
  if (T) {
    message("Use boot to calculate average C-index ...")
    ## use boot to calculate average C-index and 95\\% CI
    fun_c_index <- function(data, indices) {
      data <- data[ indices, ]
      fit <- rms::lrm(formula, data = data)
      fit$stats[ c("C", "P") ]
    }
    boot <- boot::boot(data, fun_c_index, R = B)
    boot.ci <- boot::boot.ci(boot, .95, type = "basic")
    fun <- function() {
      means <- apply(boot$t, 2, mean)
      ci <- tail(boot.ci$basic[1, ], n = 2)
      new_lich(list(`Re-sample` = B, `C-index` = means[1],
          `P-value` = means[2], "95\\% CI" = paste0(ci, collapse = " ~ "))
      )
    }
    lich <- fun()
    lich <- .set_lab(lich, "bootstrap others")
    boots <- namel(boot, boot.ci, lich)
  }
  if (F) {
    old <- getOption("prType")
    options(prType = "html")
    html <- paste0(print(lrm.supp$fit))
    coefs <- get_table.html(html)
    options(prType = old)
  } else {
    coefs <- NULL
  }
  cal <- rms::calibrate(fit, method = "boot", B = B)
  if (lang == "cn") {
    xlab <- paste0("预测", levels[2], "概率")
    ylab <- paste0("实际", levels[2], "概率")
  } else {
    xlab <- "Predicted Probability"
    ylab <- "observed Probability"
  }
  p.cal <- as_grob(
    expression(plot(cal, xlim = c(0, 1), ylim = c(0, 1),
        xlab = xlab, ylab = ylab, subtitle = F)), environment()
  )
  if (lang == "cn") {
    message("Try convert English legend as Chinese.")
    labels <- p.cal$children[[15]]$label
    if (identical(labels, c("Apparent", "Bias-corrected", "Ideal"))) {
      p.cal$children[[15]]$label <- c("表观状态", "误差纠正", "理想状态")
    }
  }
  p.cal <- wrap(p.cal, 7, 7)
  p.cal <- .set_lab(p.cal, "bootstrap calibration")
  roc <- new_roc(data[[ y ]], predict(fit), lang = lang, ...)
  .add_internal_job(.job(method = "R package `rms` used for Logistic regression and nomogram visualization"))
  namel(fit, coefs, data, levels, cal, p.cal, roc, lang, boots)
}

set_rms_datadist <- function(data) {
  .RMS_datadist <- e(rms::datadist(data))
  assign(".RMS_datadist", .RMS_datadist, envir = .GlobalEnv)
  message("A global variable defined herein: `.RMS_datadist`.")
  options(datadist = ".RMS_datadist")
}

new_nomo <- function(lrm, fun_label = lrm$levels[2], lang = lrm$lang,
  lp = F, fun_at = seq(.1, .9, by = .1))
{
  if (!is(lrm, "list")) {
    stop("the `lrm` should be object 'list' return by `new_lrm`")
  }
  lang <- match.arg(lang, c("cn", "en"))
  set_rms_datadist(lrm$data)
  fun_label <- if (lang == "en") {
    paste0("Risk of ", fun_label)
  } else {
    paste0(fun_label, "风险")
  }
  nomo <- e(rms::nomogram(lrm$fit, fun = stats::plogis,
      funlabel = fun_label, lp = lp, fun.at = fun_at))
  if (lang == "en") {
    plot(nomo, lplabel = "Linear Predictor",
      points.label = 'Points', total.points.label = 'Total Points'
    )
  } else {
    plot(nomo, lplabel = "线性预测",
      points.label = '分数', total.points.label = '总分'
    )
  }
  p.nomo <- wrap(recordPlot(), 10, .6 * length(lrm$fit$coefficients) + .3)
  p.nomo <- .set_lab(p.nomo, "nomogram plot")
  namel(p.nomo, nomo)
}

new_roc <- function(y, x, ..., plot.thres = NULL, lang = c("en", "cn"), cn.mode = c("zhen", "1-"))
{
  lang <- match.arg(lang)
  roc <- pROC::roc(y, x, ..., ci = T)
  if (T) {
    # try get p-value
    # https://stackoverflow.com/questions/61997453/how-to-get-p-value-after-roc-analysis-with-proc-package
    fun_p <- function() {
      v <- pROC::var(roc)
      b <- roc$auc - .5
      se <- sqrt(v)
      z <- (b / se)
      2 * pt(-abs(z), df = Inf)
    }
    p.value <- fun_p()
  }
  thres <- pROC::coords(roc, "best")
  if (lang == "cn") {
    cn.mode <- match.arg(cn.mode)
    if (cn.mode == "zhen") {
      xlab <- "假阳性率"
      ylab <- "真阳性率"
    } else {
      xlab <- "1-特异性"
      ylab <- "敏感性"
    }
  } else {
    xlab <- "Specificity"
    ylab <- "Sensitivity"
  }
  p.roc <- as_grob(
    expression(pROC::plot.roc(roc, print.thres = plot.thres, print.auc = T,
        print.auc.x = .4, print.auc.y = .05,
        xlab = xlab, ylab = ylab)), environment()
  )
  p.roc <- wrap(p.roc, 7, 7)
  p.roc <- .set_lab(p.roc, "ROC")
  .add_internal_job(.job(method = "R package `pROC` used for building ROC curve"))
  lich <- new_lich(
    list(AUC = as.double(roc$auc),
      "95\\% CI" = roc$ci[-2],
      `P-value` = p.value
    )
  )
  lich <- .set_lab(lich, "ROC others")
  data <- tibble::tibble(Sensitivities = roc$sensitivities, Specificities = roc$specificities)
  namel(p.roc, thres, roc, p.value, lich, data)
}

new_allu <- function(data, col.fill = 1, axes = 1:2,
  label.auto = F, label.freq = NULL, label.factor = 1, shiny = F, trunc = F)
{
  require(ggalluvial)
  fill <- colnames(data)[col.fill]
  data <- dplyr::mutate(data, fill = !!rlang::sym(fill))
  if (trunc) {
    data <- dplyr::mutate_all(data, function(x) stringr::str_trunc(x, 20))
  }
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
  scale_fill <- if (is.numeric(data[[ fill ]])) {
    scale_fill_gradientn(colors = color_set2())
  } else {
    scale_fill_manual(values = color_set(T))
  }
  geom_stratum <- if (shiny) {
    geom_stratum(aes(fill = fill))
  } else {
    geom_stratum(fill = "lightyellow")
  }
  p.alluvial <- ggplot(data, aes) +
    geom_alluvium(aes(fill = fill)) +
    geom_stratum +
    geom_text(stat = "stratum") +
    labs(fill = "", y = "") +
    scale_fill +
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
  pal <- color_set2()
  p <- ggplot(data, aes(x = reorder(var, value), y = value, fill = value)) +
    geom_col(width = .5) +
    geom_text(aes(x = var, y = value + max(value) * .01, label = value), hjust = 0, size = 3) +
    ylim(c(0, max(data$value) * 1.2)) +
    coord_flip() +
    labs(x = "", y = "") +
    scale_fill_gradient(low = pal[2], high = pal[1]) +
    rstyle("theme") +
    theme(legend.position = "")
  wrap(p, 7, nrow(data) * .5 + .5)
}

new_pie <- function(x, title = NULL, use.ggplot = T, overlap = 30,
  fun_text = function(...) ggrepel::geom_label_repel(..., max.overlaps = overlap))
{
  x <- split(x, x)
  x <- vapply(x, length, integer(1))
  if (use.ggplot) {
    data <- data.frame(var = names(x), value = unname(x))
    data <- dplyr::mutate(data,
      var = factor(var, levels = var[ order(value) ])
    )
    data <- dplyr::arrange(data, var)
    data <- dplyr::mutate(data,
      label = paste0("(", round(value / sum(value) * 100, 1), "%, ", value, ")"),
      label = paste0(var, " ", label), lab.x = .2,
      lab.y = sum(value) - (value / 2 + c(0, cumsum(value)[ -length(value) ]))
    )
    palette <- if (length(x) > 40)
      color_set(T)
    else if (length(x) > 10)
      color_set()
    else
      ggsci::pal_npg()(10)
    p <- ggplot(data, aes(x = 0L, y = value, fill = var)) +
      geom_bar(stat = 'identity', position = 'stack', width = 1) +
      fun_text(aes(x = lab.x, y = lab.y, label = label)) +
      scale_fill_manual(values = palette) +
      labs(x = '', y = '', title = '') +
      coord_polar(theta = 'y') +
      theme_minimal() +
      ggtitle(title) +
      theme(legend.position = "none",
        plot.title = element_text(size = 15, hjust = .5, vjust = -.5),
        axis.text = element_blank(),
        plot.margin = margin(if (is.null(title)) -.1 else 0, -.1, -.1, -.1, "npc"),
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
  path = .prefix("ferroptosis_2023-10-24.rds", "db"), add_internal_job = T)
{
  if (F) {
    # <http://www.zhounan.org/ferrdb/current/>
    fe_db <- list(marker = "~/Downloads/ferroptosis_marker.csv",
      driver = "~/Downloads/ferroptosis_driver.csv",
      suppressor = "~/Downloads/ferroptosis_suppressor.csv",
      unclassifier = "~/Downloads/ferroptosis_unclassified.csv",
      inducer = "~/Downloads/ferroptosis_inducer.csv",
      inhibitor = "~/Downloads/ferroptosis_inhibitor.csv",
      disease = "~/Downloads/ferroptosis_disease.csv"
    )
    fe_db <- lapply(fe_db, ftibble)
    saveRDS(fe_db, .prefix("ferroptosis_2023-10-24.rds", "db"))
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
  if (add_internal_job) {
    job <- .job(method = "Database of `FerrDb V2` used for obtaining ferroptosis regulators",
      cite = "[@FerrdbV2UpdaZhou2023]")
    .add_internal_job(job)
  }
  data <- .set_lab(data, "Ferroptosis regulators", names(data))
  lab(data) <- "Ferroptosis regulators"
  data
}

grp <- function(x, pattern, ignore.case = F, ...) {
  grep(pattern, x, ignore.case = ignore.case, ...)
}

grpl <- function(x, pattern, ignore.case = F, ...) {
  grepl(pattern, x, ignore.case = ignore.case, ...)
}

grpf <- function(x, pattern, ignore.case = F, ...) {
  x[grepl(pattern, x, ignore.case = ignore.case, ...)]
}

## igraph tips
# grid star circle

get_layout <- function(edges = NULL, layout = "grid", nodes = NULL, ...) {
  data <- as_tibble(data.frame(fast_layout(edges, layout, nodes = nodes, ...)))
  data
}

ink2d <- function(inchikey) {
  stringr::str_extract(inchikey, "^[A-Z]{14}")
}

strx <- function(...) {
  stringr::str_extract(...)
}

search.scopus <- function(data, try_format = T, sleep = 3, group.sleep = 5, n = 10, db = "scopusV2.rds", port = 7777)
{
  if (try_format) {
    .check_columns(data, c("name", "inst"))
    data <- dplyr::mutate(data,
      last.name = strx(name, "^[^,]+"),
      first.name = gs(name, "[^,]*, (.*?)", "\\1"),
      first.name = gs(first.name, "([A-Z])", " \\1"),
      first.name = gs(first.name, "^\\s", ""),
      .id = paste0(last.name, ", ", first.name, "#", inst)
    )
  }
  .check_columns(data, c("last.name", "first.name", ".id"))
  #######################
  #######################
  fun_input <- function(xpath, x) {
    ele <- link$findElement("xpath", xpath)
    ele$clearElement()
    Sys.sleep(.1)
    ele$sendKeysToElement(list(x))
    Sys.sleep(.1)
  }
  fun_search <- function(xpath) {
    ele <- link$findElement("xpath", xpath)
    ele$sendKeysToElement(list(key = "enter"))
  }
  fun_format <- function(x) {
    x <- as_tibble(x[[1]])
    colnames(x) <- gs(make.names(gs(colnames(x), "\n", "")), "^([^.]+\\.[^.]+).*", "\\1")
    x <- dplyr::select(x, Author, Documents, h.index, Affiliation, City, Country.Territory)
    x <- dplyr::filter(x, !is.na(h.index))
    x <- dplyr::mutate_all(x, function(x) gs(strx(x, "^[^\n]+"), "^\\s*", ""))
    x <- dplyr::mutate(x, Documents = as.double(Documents), h.index = as.double(h.index))
    x
  }
  db <- new_db(db, ".id")
  db <- not(db, data$.id)
  query <- dplyr::filter(data, .id %in% db@query)
  query <- dplyr::distinct(query, .id, .keep_all = T)
  #######################
  #######################
  if (nrow(query)) {
    link <- start_drive(port = port)
    Sys.sleep(3)
    link$open()
    if (is.numeric(n)) {
      group <- grouping_vec2list(1:nrow(query), n)
      group <- rep(1:length(group), lengths(group))
      lst <- split(query, group)
    } else {
      lst <- list(query)
    }
    N <- 0L
    for (query in lst) {
      N <- N + 1L
      cli::cli_h1(paste0("Group: ", N, " (", length(lst), ")"))
      Sys.sleep(sleep)
      res <- apply(query, 1, simplify = F,
        function(x) {
          x <- as.list(x)
          link$navigate("https://www.scopus.com/freelookup/form/author.uri?zone=TopNavBar&origin=NO%20ORIGIN%20DEFINED")
          fun_input("//div//input[@id='lastname']", x$last.name)
          fun_input("//div//input[@id='firstname']", x$first.name)
          fun_input("//div//input[@id='institute']", x$inst)
          Sys.sleep(1)
          fun_search("//div//button[@id='authorSubmitBtn']")
          Sys.sleep(3)
          html <- link$getPageSource()[[1]]
          table <- get_table.html(html)
          if (length(table)) {
            table <- fun_format(table)
          } else {
            table <- data.frame()
          }
          Sys.sleep(sleep)
          table
        })
      .names <- names(res) <- query$.id
      res <- frbind(res, fill = T, idcol = T)
      db <- upd(db, res, .names)
    }
    link$close()
    end_drive()
  }
  res <- dplyr::filter(db@db, .id %in% !!data$.id)
  # keep the first search result
  res <- dplyr::distinct(res, .id, .keep_all = T)
  res <- res[match(data$.id, res$.id), ]
  res <- dplyr::rename(res, `query (Name + Affiliation)` = .id)
  ################################
  ################################
  return(res)
}

read_from_compoundDiscovery <- function(file_xlsx, exdir = "compound_discovery") {
  db <- openxlsx::read.xlsx(file_xlsx, colNames = F)
  fun_format <- function(x) {
    x <- dplyr::filter(x, !is.na(X2))
    x <- dplyr::select_if(x, function(x) !all(is.na(x)))
    x <- dplyr::select(x, name = X2, rt.min = X5, formula = X6, mw = X8)
    x <- dplyr::filter(x, name != "Name")
    pos.na <- which(is.na(x$rt.min))
    pos.napre <- pos.na - 1
    str.na <- x$name[ pos.na ]
    x$name[ pos.napre ] <- paste0(x$name[ pos.napre ], x$name[ pos.na ])
    x <- dplyr::filter(x, !is.na(rt.min))
    x <- dplyr::mutate(x, en.name = gs(name, "([^\u4e00-\u9fa5]*).*", "\\1"),
      cn.name = gs(name, ".*?([\u4e00-\u9fa5].*)$", "\\1"),
      tail = strx(en.name, " [^ ]+-$"),
      tail = ifelse(is.na(tail), "", tail),
      en.name = substr(en.name, 1, nchar(en.name) - nchar(tail)),
      en.name = gs(en.name, "\\s*$", ""),
      cn.name = gs(paste0(tail, cn.name), "^\\s*", ""),
      formula = gs(formula, "\\s", ""),
      rt.min = round(as.double(rt.min), 2),
      mw = round(as.double(mw), 4)
    )
    x <- dplyr::select(x, -tail, -name)
    dplyr::relocate(x, en.name, cn.name)
  }
  table <- as_tibble(fun_format(db))
  ## peak area image
  dir.create(exdir, F)
  unzip(file_xlsx, exdir = exdir)
  fun <- function(dir) {
    files <- list.files(dir, ".png$", full.names = T)
    scale <- sapply(files, simplify = F,
      function(x) {
        try(bitmap_info(x))
      })
    scale <- frbind(scale, idcol = "file")
    scale <- dplyr::filter(scale, width > height * 15)
    scale <- dplyr::mutate(scale, order = as.integer(gs(file, ".*?([0-9]+).png", "\\1")))
    scale <- dplyr::arrange(scale, order)
  }
  scales <- fun(paste0(exdir, "/xl/media/"))
  if (nrow(table) == nrow(scales)) {
    message("Images of Peak area matched number.")
    table$file_area <- scales$file
  }
  lst <- namel(table, scales)
}

try_get_area.compoundDiscovery <- function(lstcd, res_ocr) {
  res_ocr <- lapply(res_ocr,
    function(x) {
      obj <- x$pages[[1]]$blocks[[1]]$lines[[1]]$words[[1]]
      data.frame(value = obj$value, confidence = obj$confidence)
    })
  res_ocr <- frbind(res_ocr)
  res_ocr <- dplyr::mutate(res_ocr,
    format = gs(value, "^([0-9])\\.?", "\\1."),
    format = gs(format, "B$", "8"),
    format = as.double(format))
  if (all(!is.na(res_ocr$format)) & identical(sort(res_ocr$format, decreasing = T), res_ocr$format)) {
    message("Got the expected value.")
  } else {
    Terror <<- res_ocr
    stop("Not the expected value")
  }
  lstcd$table$peak_area <- res_ocr$format
  lab(lstcd$table) <- "Identified compounds records in table CompoundDiscovery"
  lstcd
}


code <- function(lines, lang = "Bash", color = "red") {
  lines <- unlist(strsplit(lines, "\n"))
  lines <- gs(lines, "^\\s+", "")
  begin <- c(
    paste0("\\begin{tcolorbox}[colback = gray!10,",
      " colframe = ", color, "!50, width = 16cm,",
      " arc = 1mm, auto outer arc, title = {", lang, " input}]"),
    "\\begin{verbatim}")
  end <- c("\\end{verbatim}",
    "\\end{tcolorbox}")
  writeLines(c(begin, lines, end))
}

