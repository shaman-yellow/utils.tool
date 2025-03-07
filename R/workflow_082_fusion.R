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
  function(x, tissue = "Whole_Blood", from = "GTEx", type = "all_samples_link",
    dir_weights = .prefix("fusion_weights", "db"), force = FALSE)
  {
    step_message("Prepare TWAS weight file.")
    if (!dir.exists(dir_weights)) {
      dir.create(dir_weights)
    }
    file_candidates <- file.path(
      dir_weights, paste0("items_", from, ".rds")
    )
    if (!file.exists(file_candidates) || force) {
      url <- "http://gusevlab.org/projects/fusion/"
      str_html <- RCurl::getURL(url)
      lst <- get_table.html(str_html, elFun = tryGetLink.plantaedb)
      if (from == "GTEx") {
        cols <- c("Tissue", "All Samples", "link", "EUR Samples")
        isThat <- vapply(lst, 
          function(x) {
            all(cols %in% colnames(x))
          }, logical(1))
        data_candidates <- lst[ isThat ][[1]]
        colnames(data_candidates) <- c(
          "tissue", "all_samples", "all_samples_link", "eur_samples", "eur_samples_link"
        )
        data_candidates <- as_tibble(data_candidates)
        data_candidates <- dplyr::mutate(
          data_candidates, dplyr::across(
            dplyr::ends_with("link"),
            function(x) sub(".* ### ", "", x)
          )
        )
        saveRDS(data_candidates, file_candidates)
      } else {
        stop('from == "GTEx", no any other ...')
      }
    } else {
      data_candidates <- readRDS(file_candidates)
    }
    if (missing(tissue)) {
      message(showStrings(data_candidates[[ "tissue" ]], FALSE, FALSE))
      stop("Please provide the `Tissue`.")
    }
    candidate <- dplyr::filter(
      data_candidates, grpl(tissue, !!tissue)
    )
    if (!nrow(candidate)) {
      stop('!nrow(candidate), can not match.')
    }
    if (nrow(candidate) > 1) {
      stop('nrow(candidate) > 1, too many candidate.')
    }
    file_weight <- file.path(dir_weights, basename(candidate[[ type ]]))
    dir_untar <- tools::file_path_sans_ext(file_weight, TRUE)
    x <- methodAdd(x,
      "获取 {candidate$tissue} 组织的表达权重文件 (Expression Weights) (<{candidate[[type]]}>)。"
    )
    if (!dir.exists(dir_untar)) {
      if (!file.exists(file_weight)) {
        download.file(candidate[[ type ]], file_weight)
      }
      message(glue::glue("Decompressed '{file_weight}' to '{dir_untar}'."))
      utils::untar(
        normalizePath(file_weight, mustWork = FALSE), 
        exdir = normalizePath(dir_untar)
      )
    }
    file_pos <- file.path(
      dir_untar, paste0(basename(dir_untar), ".pos")
    )
    if (!file.exists(file_pos)) {
      stop('!file.exists(file_pos), file missing?')
    }
    x$weight_pos <- file_pos
    x$weight_dir <- dir_untar
    return(x)
  })

setMethod("step3", signature = c(x = "job_fusion"),
  function(x, chrs = "all", dir_output = "fusion", 
    cl = 3L, use = c("adjust.P", "P"), use.cut = .05)
  {
    step_message("TWAS: Performing the expression imputation.")
    if (!requireNamespace("plink2R", quietly = TRUE)) {
      stop('!requireNamespace("plink2R", quietly = TRUE), not installed?')
    }
    dir_output <- touch_dir(dir_output)
    pg_dir <- pg("fusion")
    file_script <- file.path(pg_dir, "FUSION.assoc_test.R")
    dir_ref_ld <- file.path(pg_dir, "LDREF")
    allChrs <- gs(list.files(dir_ref_ld, "bim$"), ".*\\.([0-9]+)\\.[^0-9]*", "\\1")
    if (identical(chrs, "all")) {
      message(glue::glue("Use all chromosome: {bind(allChrs)}"))
      chrs <- allChrs
    } else if (!any(as.character(chrs) %in% allChrs)) {
      stop('!any(as.character(chrs) %in% allChrs), some not match?')
    }
    ref_ld_chr <- file.path(dir_ref_ld, "1000G.EUR.")
    file_sumstats <- x$file_sumstats
    outputs <- pbapply::pblapply(chrs, cl = cl,
      function(chr) {
        output <- file.path(
          dir_output, paste0("chr_", chr, ".dat")
        )
        if (!file.exists(output)) {
          cdRun(glue::glue("Rscript {file_script} \\
              --sumstats {file_sumstats} \\
              --weights {x$weight_pos} \\
              --weights_dir {x$weight_dir} \\
              --ref_ld_chr {ref_ld_chr} \\
              --chr {chr} \\
              --out {output}"
              ))
          if (!file.exists(output)) {
            stop('!file.exists(output), not successfully?')
          }
        } else {
          message(glue::glue("File {output} exists, skip imputation."))
        }
        return(output)
      })
    outputs <- frbind(lapply(outputs, 
      function(file) {
        x <- ftibble(file)
        dplyr::mutate(x, TWAS.P.adjust = p.adjust(TWAS.P, "fdr"))
      }))
    outputs <- map_gene(outputs, "ID", "ENSEMBL")
    use <- match.arg(use)
    use <- switch(use, adjust.P = "TWAS.P.adjust", P = "TWAS.P")
    outputs <- dplyr::arrange(outputs, !!rlang::sym(use))
    sigOutputs <- dplyr::filter(outputs, !!rlang::sym(use) < !!use.cut)
    sigOutputs <- setLegend(sigOutputs, "为 TWAS 显著统计表 ({use} &lt; {use.cut})。")
    feature(x) <- sigOutputs$SYMBOL[ !is.na(sigOutputs$SYMBOL) ]
    outputs <- setLegend(outputs, "为 TWAS 基因与疾病关联性统计结果，显著性 TWAS.P.adjust 由 TWAS.P 以染色体对应的基因数 (Expression Weights) FDR 校正计算。该表格的解释请参考 <http://gusevlab.org/projects/fusion/>。")
    x <- tablesAdd(x, TWAS_statistic = outputs, TWAS_significant = sigOutputs)
    x <- methodAdd(x, "以 `FUSION` {cite_show('IntegrativeAppGusev2016')} (<http://gusevlab.org/projects/fusion/>) 进行 TWAS 预测，得到基因与疾病之间的关联统计。")
    chrs <- sort(as.integer(chrs))
    if (identical(chrs, min(chrs):max(chrs))) {
      chrs <- paste0(min(chrs), "-", max(chrs))
    }
    x <- snapAdd(x, "以 `FUSION` 预测基因与疾病之间的关联 (chromosome: {bind(chrs)})。(TWAS 能够提供 SNP 如何通过调控基因表达来影响表型的机制)")
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_fusion"),
  function(x, ref, pattern = NULL, cut.pct = .1, use = "contrasts", 
    plot_heatmap = FALSE, group.by = x$group.by)
  {
    message("Filter by Cell (`pattern`).")
    if (ref@step < 3L) {
      stop('ref@step < 3L...')
    }
    if (is.null(twas <- ref@tables$step3$TWAS_significant)) {
      stop('is.null(ref@table$step3$TWAS_significant).')
    }
    twas <- dplyr::select(twas, SYMBOL, TWAS.Z, TWAS.P)
    if (is.null(degs <- x[[ use ]])) {
      stop('is.null(x[[ use ]]), no this DEGs list computed by `diff`?')
    }
    x <- snapAdd(
      x, "分析{snap(feature(ref))}是否在各个细胞群中差异表达。", step = "fusion"
    )
    x$.map_heading <- "Seurat 细胞群中的 TWAS 风险相关基因"
    fusion_degsAll <- degs <- tbmerge(degs, twas, by.x = "gene", by.y = "SYMBOL")
    fusion_degsAll <- .set_lab(fusion_degsAll, sig(x), "TWAS associated genes of Cell Cluster DEGs")
    x$fusion_degsAll <- setLegend(fusion_degsAll, "TWAS 风险相关的细胞群 DEGs。")
    cell <- "细胞群"
    if (!is.null(pattern)) {
      degs <- dplyr::filter(degs, grpl(contrast, pattern))
      cell <- "{pattern}"
    }
    filter <- ""
    if (!is.null(cut.pct)) {
      degs <- dplyr::filter(degs, pct.1 > cut.pct | pct.2 > cut.pct)
      filter <- glue::glue("(筛选至少有检出率 {cut.pct * 100}% 的细胞表达该基因)")
    }
    degs <- .set_lab(degs, sig(x), "Filtered TWAS associated genes of Cell Cluster DEGs")
    degs <- setLegend(degs, "TWAS 风险相关的{cell} DEGs {filter}。")
    x <- snapAdd(
      x, "筛选 TWAS 风险相关的 DEGs {filter}：{try_snap(degs, 'contrast', 'gene')}。",
      step = "fusion"
    )
    x$fusion_degs <- degs
    if (plot_heatmap) {
      x$fusion_p.hp <- e(Seurat::DoHeatmap(object(x), features = unique(degs$gene),
          group.by = x$group.by, raster = TRUE, group.colors = color_set(),
          slot = "data", label = FALSE))
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

