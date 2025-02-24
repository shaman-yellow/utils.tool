# ==========================================================================
# workflow of annova
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_annova <- setClass("job_annova", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: "),
    cite = "[@AnnovarFunctiWang2010]",
    method = "`Annova` used for mutation annotation",
    tag = "annova",
    analysis = "Annota 变异注释"
    ))

setGeneric("asjob_annova", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_annova"))

setMethod("step0", signature = c(x = "job_annova"),
  function(x){
    step_message("Prepare your data with function `job_annova`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_annova"),
  function(x){
    step_message("Quality control (QC).
      This do:
      "
    )
    return(x)
  })

annovar_variation <- function(file_input,
  prefix = tools::file_path_sans_ext(basename(file_input)),
  ref = "hg19", dbSNP = "avsnp151", db_species = "humandb")
{
  annovar_path <- pg("annovar")
  pg <- file.path(annovar_path, "annotate_variation.pl")
  if (!file.exists(pg)) {
    stop(glue::glue('!file.exists("{pg}"), no annovar downloaded?'))
  }
  db_path <- file.path(annovar_path, db_species)
  db_files <- file.path(
    db_path, paste0(ref, "_", dbSNP, ".txt.idx")
  )
  if (!file.exists(db_files)) {
    message("No buildver file found, download that.")
    cdRun(glue::glue("{pg} -downdb -webfrom annovar -buildver {ref} {dbSNP} {db_path} --verbose"))
  }
  types <- c("dropped", "filtered")
  outputs <- paste0(prefix, glue::glue(".{ref}_{dbSNP}_{types}"))
  if (!any(file.exists(outputs))) {
    cdRun(
      glue::glue("{pg} {file_input} {db_path} -filter -build {ref} -dbtype {dbSNP} --outfile {prefix}")
    )
  }
  outputs
}
