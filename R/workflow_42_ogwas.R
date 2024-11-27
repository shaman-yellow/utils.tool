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

job_ogwas <- function(traits, ...)
{
  x <- .job_ogwas()
  object(x) <- as_tibble(e(MungeSumstats::find_sumstats(traits = traits)))
  lab(object(x)) <- "Traits in Open GWAS"
  x
}

setMethod("step0", signature = c(x = "job_ogwas"),
  function(x){
    step_message("Prepare your data with function `job_ogwas`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_ogwas"),
  function(x, ids, ref_genome = "GRCH38",
    vcf_dir = .prefix("ogwas_vcf", "db"), save_dir = .prefix("ogwas_data", "db"))
  {
    step_message("Import gwas summary data")
    if (!dir.exists(vcf_dir)) {
      dir.create(vcf_dir)
    }
    if (!dir.create(save_dir)) {
      dir.create(save_dir)
    }
    .suggest_bio_package("SNPlocs.Hsapiens.dbSNP155.GRCh38")
    .suggest_bio_package("BSgenome.Hsapiens.NCBI.GRCh38")
    x$db <- e(MungeSumstats::import_sumstats(ids = ids, ref_genome = ref_genome,
        vcf_dir = vcf_dir, save_dir = save_dir))
    return(x)
  })
