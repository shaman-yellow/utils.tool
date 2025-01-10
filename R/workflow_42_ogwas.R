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
    vcf_dir = .prefix("ogwas_vcf", "db"), save_dir = .prefix("ogwas_data", "db"))
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
    return(x)
  })
