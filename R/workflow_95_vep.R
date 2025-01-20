# ==========================================================================
# workflow of vep
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_vep <- setClass("job_vep", 
  contains = c("job"),
  prototype = prototype(
    pg = "vep",
    info = c("https://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html"),
    cite = "[@The_Ensembl_Var_McLare_2016]",
    method = "",
    tag = "vep",
    analysis = "VEP 变异注释"
    ))

job_vep <- function(input, version = 113, species = "homo_sapiens",
  assembly = c("GRCh37", "GRCH38"), 
  save = "vep", hasID = FALSE, check_allele_columns = TRUE, sort = TRUE)
{
  if (dir.exists(save)) {
    if (sureThat("dir.exists(save) == TRUE, remove That?")) {
      unlink(save, TRUE)
    }
  }
  dir.create(save, FALSE)
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    if (file.size(input) > 1e8) {
      message("Detected large files, loading...")
    }
    input <- ftibble(input)
  }
  # chrom, pos, ID, ref, alt
  if (is.data.frame(input)) {
    if (!hasID) {
      message("The starts of columns should be: chrom, pos, ref, alt.")
      input <- dplyr::select(input, 1:4)
      input <- dplyr::mutate(input, ID = seq_len(nrow(input)), .after = 2)
    } else {
      message("The starts of columns should be: chrom, pos, ID, ref, alt.")
      input <- dplyr::select(input, 1:5)
    }
    if (!is.integer(input[[2]])) {
      stop('!is.integer(input[[2]]), "pos" columns should be "integer".')
    }
    if (check_allele_columns) {
      input_head <- head(input)
      pattern <- "[ATCG-]"
      if (!(all(grpl(input_head[[4]], pattern)) && all(grpl(input_head[[5]], pattern)))) {
        stop(
          '!(all(grpl(input_head[[4]], pattern)) && all(grpl(input_head[[5]], pattern))), not ATCG?'
        )
      }
      if (grpl(colnames(input)[4], "effect.*allele", TRUE) &&
        grpl(colnames(input)[5], "other.*allele", TRUE))
      {
        message("Detected columns names of 'allele'; Relocate columns of 4 and 5.")
        input <- dplyr::relocate(input, 1, 2, 3, 5, 4)
      }
    }
    colnames(input) <- c("#CHROM", "POS", "ID", "REF", "ALT")
  }
  if (sort) {
    message("Sort by chromosome and by location.")
    input <- dplyr::arrange(input, `#CHROM`, `POS`)
  }
  x <- .job_vep(object = input)
  x$version <- version
  x$species <- species
  assembly <- match.arg(assembly)
  x$assembly <- assembly
  x$save <- save
  return(x)
}

setMethod("step0", signature = c(x = "job_vep"),
  function(x){
    step_message("Prepare your data with function `job_vep`.")
  })

setMethod("step1", signature = c(x = "job_vep"),
  function(x, force = FALSE, fork = 4,
    fields = c("Uploaded_variation", "Existing_variation"), ...)
  {
    step_message("Run annotation.")
    vcf <- file.path(x$save, "input.vcf")
    writeVcf <- TRUE
    if (file.exists(vcf) && (writeVcf <- sureThat("Removing exists file?"))) {
      file.remove(vcf)
    }
    if (writeVcf) {
      vcf_heading <- c(
        "##fileformat=VCFv4.0", paste0(colnames(object(x)), collapse = "\t")
      )
      writeLines(vcf_heading, vcf)
      data.table::fwrite(
        object(x), vcf, append = TRUE, sep = "\t", quote = FALSE
      )
    }
    x$vcf <- vcf
    x$output <- file.path(x$save, "vep_output.tsv")
    x$fields <- fields
    vep(
      x$vcf, x$output, x$version, x$species,
      x$assembly, force = force, fork = fork, fields = x$fields, ...
    )
    # x$anno <- ftibble(x$output, skip = "#")
    return(x)
  })

setMethod("step2", signature = c(x = "job_vep"),
  function(x, mode = "Existing_variation", id_which = 3L, mode_which = 2L){
    step_message("Get annotation results.")
    mode <- match.arg(mode)
    anno_id <- data.table::fread(
      x$output, skip = "#", sep = "\t", select = id_which
    )
    if (mode == "Existing_variation") {
      if (!any(x$fields == "Existing_variation")) {
        stop('!any(x$fields == "Existing_variation")')
      }
      anno_var <- data.table::fread(
        x$output, skip = "#", sep = "|", select = mode_which
      )
      if (nrow(anno_id) != nrow(anno_var)) {
        stop('nrow(anno_id) != nrow(anno_var).')
      }
      anno_var <- dplyr::rename(anno_var, rsID = 1)
      anno_var <- dplyr::mutate(
        anno_var, rsID = strx(rsID, "rs[0-9]+")
      )
      anno <- dplyr::bind_cols(anno_id, anno_var)
    }
    x$anno <- anno
    return(x)
  })

vep <- function(input, output = "vep_output.tsv",
  version = 113, species = "homo_sapiens",
  assembly = "GRCh37", fields = c(
    "Uploaded_variation", "Existing_variation"
  ), force = FALSE, fork = 4)
{
  pg <- pg('vep')
  dir <- pg('vep_cache')
  input <- normalizePath(input)
  if (!file.exists(output) || force) {
    cdRun(
      glue::glue(
        "{pg} -i {input} -o {output} --offline \\
        --dir {dir} --species {species} --db_version {version} \\
        --assembly {assembly} \\
        --check_existing --force_overwrite --fork {fork} \\
        --vcf --fields '{paste0(fields, collapse = ',')}' \\
        --no_stats --verbose"
      )
    )
  } else {
    message("VEP output exists, skip running.")
  }
  return(output)
}
