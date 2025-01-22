# ==========================================================================
# workflow of ogwas
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ogwas <- setClass("job_ogwas", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@MungesumstatsMurphy2021; @TheMrcIeuOpeElswor2020; @TheMrBasePlaHemani2018]",
    method = "R package `MungeSumstats` used for downloading or formatting GWAS summary data (from Open GWAS)",
    tag = "gwas",
    analysis = "MungeSumstats 获取 GWAS 数据"
    ))

job_ogwas <- function(traits, api_token = NULL, dir = .prefix("ogwas", "db"), 
  force = FALSE)
{
  if (missing(traits)) {
    stop('missing(traits), please give the pattern to match.')
  }
  if (is.null(api_token)) {
    if (is.null(api_token <- getOption("gwas_token", NULL))) {
      stop('is.null(getOption("gwas_token", NULL))')
    }
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  file_ogwasInfo <- file.path(dir, "ogwasinfo.rds")
  methodAdd_onExit("x", "以 R 包 `MungeSumstats` ({packageVersion('MungeSumstats')}) {cite_show('MungesumstatsMurphy2021')} 和 R 包 `ieugwasr` 获取 Open GWAS 的可用数据。")
  if (file.exists(file_ogwasInfo) && !force) {
    ogwasInfo <- readRDS(file_ogwasInfo)
  } else {
    ogwasInfo <- try(
      e(ieugwasr::gwasinfo(opengwas_jwt = api_token)), TRUE
    )
    if (inherits(ogwasInfo, "try-error")) {
      stop(
        'inherits(ogwasInfo, "try-error"), maybe the token Expired?\n',
        'Move to <https://api.opengwas.io/profile/>.'
      )
    }
    if (!is(ogwasInfo, "data.frame")) {
      stop('!is(ogwasInfo, "data.frame"), not valid.')
    }
    cli::cli_alert_info("Sys.setenv")
    saveRDS(ogwasInfo, file_ogwasInfo)
  }
  Sys.setenv(OPENGWAS_JWT = api_token)
  ogwasInfo <- as_tibble(data.frame(ogwasInfo))
  x <- .job_ogwas()
  x$ogwasInfo <- ogwasInfo
  if (length(traits) > 1) {
    traits <- paste0(traits, collapse = "|")
  }
  x <- snapAdd(x, "获取 Open GWAS 的可用数据，匹配 {traits} (trait)。")
  res <- dplyr::filter(ogwasInfo, grpl(trait, !!traits, TRUE))
  if (!nrow(res)) {
    stop('!nrow(res), found nothing.')
  }
  lab(res) <- "Traits found in Open GWAS"
  res <- setLegend(res, "在 Open GWAS 中匹配到的可用数据集 (GWAS统计数据)。")
  object(x) <- res
  return(x)
}

setMethod("step0", signature = c(x = "job_ogwas"),
  function(x){
    step_message("Prepare your data with function `job_ogwas`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_ogwas"),
  function(x, ids, which = NULL, ref_genome = "GRCH38",
    vcf_dir = .prefix("ogwas_vcf", "db"),
    save_dir = .prefix("ogwas_data", "db"), try_catalog = TRUE)
  {
    step_message("Import gwas summary data")
    if (!dir.exists(vcf_dir)) {
      dir.create(vcf_dir)
    }
    if (!dir.create(save_dir)) {
      dir.create(save_dir)
    }
    if (!is.null(which) && missing(ids)) {
      ids <- object(x)$id[which]
    }
    .suggest_bio_package("SNPlocs.Hsapiens.dbSNP155.GRCh38")
    .suggest_bio_package("BSgenome.Hsapiens.NCBI.GRCh38")
    .suggest_bio_package("GenomicFiles")
    x$db <- e(MungeSumstats::import_sumstats(ids = ids, ref_genome = ref_genome,
        vcf_dir = vcf_dir, save_dir = save_dir))
    if (!inherits(x$db[[1]], "simpleError")) {
      x <- methodAdd(x, "以 `MungeSumstats::import_sumstats` 导入 GWAS summary 数据。")
    } else {
      message('!inherits(x$db[[1]], "simpleError"), can not import via `MungeSumstats::import_sumstats`.')
      if (try_catalog && length(ids) == 1 && grpl(ids, "GCST")) {
        id <- strx(ids, "GCST[0-9]+")
        cli::cli_alert_info("Try download from GWAS Catalog.")
        x$db <- .download_full_summary_from_catalog(id)
        x <- methodAdd(
          x, "从 GWAS Catalog (<https://www.ebi.ac.uk/>) 下载 {id} 的
          Full Summary Statistic 数据 (<{attr(x$db, 'url')}>)。"
        )
      }
    }
    object(x) <- dplyr::filter(object(x), id %in% !!ids)
    return(x)
  })

.download_full_summary_from_catalog <- function(id,
  db_dir = .prefix("GCST", "db"))
{
  db_dir <- file.path(db_dir, id)
  if (!dir.exists(db_dir)) {
    dir.create(db_dir, recursive = TRUE)
  }
  base_url <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
  parse_as_xml <- function(url) {
    XML::htmlParse(RCurl::getURL(url))
  }
  xml <- parse_as_xml(base_url)
  find <- function(xml, xpath) {
    unlist(XML::xpathApply(xml, xpath, XML::xmlValue))
  }
  dirs <- find(xml, "//td/a[@href]")
  dirs <- dirs[ grpl(dirs, "^GCST") ]
  categories <- tibble::tibble(dirs = dirs)
  categories <- tidyr::separate(
    categories, dirs, c("from", "to"), "-", remove = FALSE
  )
  get_seq <- function(x) as.integer(strx(x, "[0-9]+"))
  categories <- dplyr::mutate(
    categories, from = get_seq(from), to = get_seq(to)
  )
  seq_id <- get_seq(id)
  cate_id <- dplyr::filter(categories, seq_id > from, seq_id < to)
  if (!nrow(cate_id)) {
    stop('!nrow(cate_id), can not match the data categary.')
  } else if (nrow(cate_id) > 1) {
    stop('nrow(cate_id) > 1, too many matched.')
  }
  dir_url <- paste0(base_url, "/", cate_id$dirs, "/", id, "/")
  xml <- parse_as_xml(dir_url)
  file <- find(xml, "//td/a[@href]")
  file <- grpf(file, "\\.tsv\\.gz$")
  if (!length(file)) {
    stop('!length(file), no such file.')
  }
  if (length(file) > 1) {
    file <- file[menu(file, title = "Download which?")]
  }
  file_local <- file.path(db_dir, file)
  file_url <- paste0(dir_url, "/", file)
  if (!file.exists(file_local)) {
    download.file(file_url, file_local)
  }
  res <- nl(id, file_local)
  attr(res, "url") <- file_url
  return(res)
}
