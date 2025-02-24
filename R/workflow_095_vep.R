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

setGeneric("asjob_vep", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_vep"))

setMethod("asjob_vep", signature = c(x = "job_ogwas"),
  function(x, which = 1L, ...){
    metadata <- dplyr::slice(object(x), which)
    file <- x$db[[which]]
    if (!length(file) || length(file) > 1) {
      stop('!length(file) || length(file) > 1.')
    }
    x <- job_vep(
      file, save = file.path(dirname(file), "vep"), ...
    )
    x$metadata <- metadata
    x$from_ogwas <- TRUE
    return(x)
  })

job_vep <- function(input, version = 113, species = "homo_sapiens",
  assembly = c("GRCh37", "GRCH38"), 
  save = "vep", hasID = FALSE, check_allele_columns = TRUE, 
  sort = TRUE, check = TRUE)
{
  if (dir.exists(save)) {
    if (sureThat("dir.exists(save) == TRUE, remove That?")) {
      unlink(save, TRUE)
    }
  }
  dir.create(save, FALSE)
  add_to_file <- NULL
  if (is.character(input) && length(input) == 1 && file.exists(input)) {
    temp <- ftibble(input, nrow = 100)
    print(temp)
    if (check) {
      if (hasID) {
        message("The starts of columns should be: chrom, pos, ID, ref, alt.")
      } else {
        message("The starts of columns should be: chrom, pos, ref, alt.")
      }
      if (!sureThat("Is this data suitable?")) {
        stop("...")
      }
    }
    add_to_file <- input
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
  x$add_to_file <- add_to_file
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
    x <- methodAdd(x, "以 `Ensembl-vep` 对 SNP 注释 {cite_show('The_Ensembl_Var_McLare_2016')}，获取 rsID。")
    x <- snapAdd(
      x, "以 `VEP` (根据 chromosome, position, other allele (REF), effect allele (ALT)) 获取 rsID，"
    )
    # x$anno <- ftibble(x$output, skip = "#")
    return(x)
  })

setMethod("step2", signature = c(x = "job_vep"),
  function(x, mode = "Existing_variation", id_which = 3L){
    step_message("Get annotation results.")
    mode <- match.arg(mode)
    anno_id <- data.table::fread(
      x$output, skip = "#", sep = "\t", select = id_which
    )
    if (mode == "Existing_variation") {
      if (!any(x$fields == "Existing_variation")) {
        stop('!any(x$fields == "Existing_variation")')
      }
      x$file_variation <- file.path(x$save, "Existing_variation.tsv")
      if (!file.exists(x$file_variation)) {
        if (Sys.which("sed") == "") {
          stop('Sys.which("sed") == "", can not found `sed`.')
        }
        cdRun(glue::glue("sed 's/.*|rs\\([0-9]*\\).*/rs\\1/' {x$output} > {x$file_variation}"))
      } else {
        message(
          glue::glue("file.exists: {x$file_variation}, skip formatting.")
        )
      }
      anno_var <- data.table::fread(x$file_variation)
      for (i in seq_len(100)) {
        if (!grpl(anno_var[i, 1], "^#")) {
          break
        }
      }
      anno_var <- anno_var[-(seq_len(i - 1)), ]
      if (nrow(anno_id) != nrow(anno_var)) {
        stop('nrow(anno_id) != nrow(anno_var).')
      }
      colnames(anno_var) <- "rsID"
      anno <- dplyr::bind_cols(anno_id, anno_var)
    }
    x$anno <- anno
    return(x)
  })

setMethod("step3", signature = c(x = "job_vep"),
  function(x, add_to_file = x$add_to_file,
    newfile = add_filename_suffix(add_to_file, "rsID"),
    col = "rsID", clear = TRUE)
  {
    if (file.exists(newfile)) {
      message(glue::glue("file.exists: {newfile}."))
      x$newfile <- newfile
      if (clear) {
        object(x) <- NULL
        x$anno <- NULL
      }
      return(x)
    }
    if (is.null(x$anno)) {
      stop('is.null(x$anno), no annotation data.')
    }
    message("Sort `x$anno` by `ID` column.")
    x$anno <- dplyr::arrange(x$anno, ID)
    if (is.null(add_to_file)) {
      stop('is.null(add_to_file), no file?')
    }
    if (!file.exists(add_to_file)) {
      stop('!file.exists(add_to_file).')
    }
    input <- data.table::fread(add_to_file)
    message("Read successfully.")
    if (nrow(input) != nrow(x$anno)) {
      stop('nrow(input) != nrow(x$anno), can not match row numbers.')
    }
    message(glue::glue("Add '{col}' into '{add_to_file}'."))
    input[[ col ]] <- x$anno[[ col ]]
    message("Write output as '{newfile}'...")
    data.table::fwrite(input, newfile, sep = "\t")
    x$newfile <- newfile
    if (clear) {
      object(x) <- NULL
      x$anno <- NULL
    }
    return(x)
  })

setGeneric("asjob_fusion", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_fusion"))

setMethod("asjob_fusion", signature = c(x = "job_vep"),
  function(x, file_gwas = x$newfile, SNP = "rsID", 
    A1 = "effect_allele", A2 = "other_allele", Z = "Z", 
    BETA = "beta", P = "p_value", N = NULL)
  {
    args <- x@params[ names(x@params) %in% c("assembly", "species", "metadata") ]
    if (!is.null(args$metadata) && !is.null(x$from_ogwas) && is.null(N)) {
      N <- args$metadata$sample_size
      message(glue::glue("Use N as sample_size: {N}"))
    }
    x <- job_fusion(
      file_gwas = file_gwas, SNP = SNP, A1 = A1, A2 = A2, 
      Z = Z, N = N, P = P, BETA = BETA
    )
    x@params <- c(x@params, args)
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
