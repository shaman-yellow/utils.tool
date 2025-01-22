# ==========================================================================
# workflow of fusion
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_fusion <- setClass("job_fusion", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "fusion",
    info = c("Tutorial: http://gusevlab.org/projects/fusion/"),
    cite = "[@IntegrativeAppGusev2016]",
    method = "",
    tag = "fusion",
    analysis = "FUSION TWAS全转录组关联研究"
    ))

job_fusion <- function(file_gwas, SNP = "rsID",
  A1 = "effect_allele", A2 = "other_allele", Z = "Z", 
  BETA = "beta", P = "p_value", N = NULL, nThread = 4L, use = NULL)
{
  args <- as.list(environment())
  x <- .job_fusion(params = args)
  if (!file.exists(file_gwas)) {
    stop('!file.exists(file_gwas).')
  }
  message(glue::glue("Read heading rows of the file: {file_gwas}."))
  temp <- ftibble(file_gwas, nrow = 10)
  print(temp)
  .check_columns(temp, c(SNP, A1, A2, P), file_gwas)
  if (is.null(use)) {
    if (any(colnames(temp) == BETA)) {
      x$use <- use <- BETA
    } else if (any(colnames(temp) == "Z")) {
      x$use <- use <- "Z"
    } else {
      stop('any(colnames(temp) == BETA); any(colnames(temp) == Z)')
    }
    message(glue::glue("Use 'signed-sumstats' as: '{use}'"))
  }
  newfile <- add_filename_suffix(file_gwas, "format")
  if (use == "Z") {
    if (is.function(Z)) {
      if (!file.exists(newfile)) {
        expr <- parse(text = glue::glue("data[, Z := {Z()}]"))
        message(glue::glue("Test function in heading of file..."))
        data <- data.table::as.data.table(temp)
        res <- try(eval(expr), FALSE)
        if (inherits(res, "try-error")) {
          stop('inherits(res, "try-error"), can not found the columns?')
        }
        message(glue::glue("Successfully for heading of file."))
        print(as_tibble(data))
        message(glue::glue("Computing Z for: {file_gwas}"))
        data <- data.table::fread(file_gwas)
        eval(expr)
        message(glue::glue("Saving to {newfile}"))
        data <- dplyr::select(
          data, !!rlang::sym(SNP), 
          Z, !!rlang::sym(A1), !!rlang::sym(A2), !!rlang::sym(P)
        )
        data.table::fwrite(data, newfile, nThread = nThread, sep = "\t")
      } else {
        message(glue::glue("zScore file exists: {newfile}."))
      }
      message(glue::glue("Use file of: {newfile}"))
      x$Z <- Z <- "Z"
    } else if (!any(colnames(temp) == Z)) {
      message(glue::glue("No columns of '{Z}', test columns of `beta` and `se`."))
      if (any(colnames(temp) %in% c("beta", "standard_error"))) {
        message("Detected `beta` or `se`, you can format to get `Z`.")
      }
      stop("...")
    }
  }
  if (!file.exists(newfile)) {
    data <- data.table::fread(file_gwas)
    message(glue::glue("Saving to {newfile}"))
    data <- dplyr::select(
      data, !!rlang::sym(SNP), 
      !!rlang::sym(use), !!rlang::sym(A1), !!rlang::sym(A2),
      !!rlang::sym(P)
    )
    data.table::fwrite(data, newfile, nThread = nThread, sep = "\t")
  }
  x$file_gwas <- newfile
  return(x)
}

setMethod("step0", signature = c(x = "job_fusion"),
  function(x){
    step_message("Prepare your data with function `job_fusion`.")
  })

setMethod("step1", signature = c(x = "job_fusion"),
  function(x){
    step_message("Format as .sumstats file.")
    prefix <- tools::file_path_sans_ext(x$file_gwas, TRUE) 
    x$file_sumstats <- paste0(prefix, ".sumstats.gz")
    if (!file.exists(x$file_sumstats)) {
      format_sumstats(x$file_gwas, x$SNP, x$A1, x$A2, x$use, x$N, prefix)
    }
    x <- methodAdd(x, "以 Python 工具 LDSC (`munge_sumstats.py`) (<https://github.com/bulik/ldsc>) {cite_show('LD_Score_regres_Bulik_2015')} 将 GWAS summary 文件检查并格式化为 .sumstats 格式。")
    x <- snapAdd(x, "以 `ldsc` 将 GWAS summary 转化为 .sumstats 格式。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_fusion"),
  function(x, tissue = "Whole_Blood", from = "GTEx",
    weights_dir = .prefix("fusion_weights", "db"))
  {
    step_message("Compute TWAS score.")
    if (!dir.exists(weights_dir)) {
      dir.create(weights_dir)
    }
    file_candidates <- file.path(
      weights_dir, paste0("items_", from, ".rds")
    )
    if (!file.exists(file_candidates)) {
      url <- "http://gusevlab.org/projects/fusion/"
      str_html <- RCurl::getURL(url)
      lst <- get_table.html(str_html, elFun = tryGetLink.plantaedb)
      if (from == "GTEx") {
        cols <- c("Tissue", "All Samples", "link", "EUR Samples")
      } else {
        stop('from == "GTEx", no any other ...')
      }
      isThat <- vapply(lst, 
        function(x) {
          all(cols %in% colnames(x))
        }, logical(1))
      data_candidates <- lst[ isThat ][[1]]
      saveRDS(data_candidates, file_candidates)
    } else {
      data_candidates <- readRDS(file_candidates)
    }
    x$data_candidates <- data_candidates
    if (missing(tissue)) {
      message(showStrings(data_candidates[[ "Tissue" ]], FALSE, FALSE))
      stop("Please provide the `Tissue`.")
    }
    return(x)
  })

format_sumstats <- function(input, SNP, A1, A2, Sign, N, output = "output")
{
  # https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format
  pg_py <- pg("ldscPython")
  pg_path <- pg("ldsc")
  pg_munge <- file.path(pg_path, "munge_sumstats.py")
  input <- normalizePath(input)
  # columns: SNP, N, Z (sign), A1, A2
  cdRun(glue::glue("{pg_py} {pg_munge} --sumstats {input} \\
      --snp {SNP} \\
      --a1 {A1} \\
      --a2 {A2} \\
      --signed-sumstats {Sign},0 \\
      --N {N} \\
      --out {output}"))
}

