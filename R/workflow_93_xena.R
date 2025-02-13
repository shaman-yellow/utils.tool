# ==========================================================================
# workflow of xena
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_xena <- setClass("job_xena", 
  contains = c("job"),
  prototype = prototype(
    pg = "xena",
    info = c("https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html#workflow"),
    cite = "[@The_UCSCXenaToo_Wang_2019]",
    method = "",
    tag = "xena",
    analysis = "UCSCXenaTools 癌症相关数据获取"
    ))

job_xena <- function(mode = "TcgaTargetGtex", dir = .prefix("xena", "db"), use_hiplot = TRUE)
{
  mode <- match.arg(mode)
  options(use_hiplot = use_hiplot)
  dir.create(dir, FALSE)
  x <- .job_xena(params = list(mode = mode, dir = dir))
  if (mode == "TcgaTargetGtex") {
    files <- file.path(dir, paste0("TcgaTargetGtex_gene_expected_count", c("", ".gz")))
    if (!any(file.exists(files))) {
      data <- dplyr::filter(
        UCSCXenaTools::XenaData, grpl(
          XenaDatasets, "TcgaTargetGtex", TRUE
        ), DataSubtype %in% c("phenotype", "gene expression RNAseq")
      )
      data <- dplyr::filter(data, grpl(XenaDatasets, "gene_expected_count|phenotype"))
      xena <- e(UCSCXenaTools::XenaGenerate(data))
      query <- e(UCSCXenaTools::XenaQuery(xena))
      e(
        UCSCXenaTools::XenaDownload(
          query, dir, method = "wget", extra = "-c", force = FALSE, download_probeMap = TRUE
        )
      )
    }
    file_clinical <- file.path(dir, "TcgaTargetGTEX_phenotype.txt.gz")
    x$clinical <- ftibble(file_clinical)
    file_genes <- file.path(
      dir, "probeMap", "gencode.v23.annotation.gene.probemap"
    )
    x$genes <- ftibble(file_genes)
    x$file_count <- file.path(dir, "TcgaTargetGtex_gene_expected_count.gz")
    if (.Platform$OS.type == "unix") {
      tmp <- tempfile("head_log2_count", fileext = ".tsv")
      cdRun(glue::glue("zcat {x$file_count} | head -n 10 > {tmp}"))
      x$head_counts <- ftibble(tmp)
    }
  }
  x$project <- mode
  x <- methodAdd(x, "以 R 包 `UCSCXenaTools` ({packageVersion('UCSCXenaTools')}) {cite_show('The_UCSCXenaToo_Wang_2019')} 获取 `{mode}` 类型数据。")
  x <- snapAdd(x, "获取 UCSC Xena 的 `{mode}` 数据。")
  return(x)
}

setMethod("step0", signature = c(x = "job_xena"),
  function(x){
    step_message("Prepare your data with function `job_xena`.")
  })

setMethod("step1", signature = c(x = "job_xena"),
  function(x, cancer, site, Normal = TRUE, cancer_types = c("Primary Tumor", "Metastatic"),
    group = "guess", mode = c(
      "SKCM", "OV", "COAD"
    ), add_batch = TRUE)
  {
    step_message("Filter metadata (clinical data).")
    if (!missing(mode)) {
      mode <- match.arg(mode)
      if (mode == "SKCM") {
        cancer <- "Skin Cutaneous Melanoma"
        site <- "Skin"
      } else if (mode == "OV") {
        cancer <- "Ovarian Serous Cystadenocarcinoma"
        site <- "Ovary"
      } else if (mode == "COAD") {
        cancer <- "Colon Adenocarcinoma"
        site <- "Colon"
      }
      x$project <- paste0(x$project, "-", mode)
    } else {
      x$project <- paste0(x$project, "-", cancer)
    }
    metadata <- x$clinical
    data_cancer <- dplyr::filter(
      metadata, `_primary_site` == !!site, grpl(
        detailed_category, cancer, TRUE
      )
    )
    if (!is.null(cancer_types)) {
      data_cancer <- dplyr::filter(data_cancer, `_sample_type` %in% !!cancer_types)
    }
    if (Normal) {
      data_normal <- dplyr::filter(
        metadata, `_primary_site` == !!site, `_sample_type` == "Normal Tissue"
      )
    }
    metadata <- frbind(lapply(ls(pattern = "^data_"), get, envir = environment()))
    metadata <- dplyr::filter(metadata, sample %in% !!colnames(x$head_counts))
    if (identical(group, "guess")) {
      metadata <- dplyr::mutate(
        metadata, group = ifelse(
          `_sample_type` %in% !!cancer_types, !!mode, "Normal"
        ), .after = 1
      )
    }
    if (add_batch) {
      metadata <- dplyr::mutate(metadata, batch = `_study`)
    }
    x <- snapAdd(x, "共 {nrow(metadata)} 个数据。样本组织为：{try_snap(metadata$`_primary_site`)}。样本类型为：{try_snap(metadata$`_sample_type`)}。数据来源为：{try_snap(metadata$`_study`)}。")
    x$metadata <- metadata
    return(x)
  })

setMethod("step2", signature = c(x = "job_xena"),
  function(x){
    step_message("Load the data via data.table::fread (Use: select).")
    if (!is.null(x$counts)) {
      message('!is.null(x$counts), skip read.')
    } else {
      counts <- e(
        data.table::fread(
          x$file_count, select = c("sample", x$metadata$sample), 
          header = TRUE, sep = "\t"
        )
      )
      x$counts <- as_tibble(counts)
    }
    return(x)
  })

setMethod("asjob_limma", signature = c(x = "job_xena"),
  function(x, raw_counts = TRUE){
    message("Convert job_xena as job_limma.")
    counts <- x$counts
    project <- x$project
    if (raw_counts) {
      snapAdd_onExit(
        "x", "将 {project} 数据转化为 counts 型数据 (原数据为 log2(counts + 1)) 。"
      )
      counts <- dplyr::mutate(
        counts, dplyr::across(dplyr::where(is.double), function(x) 2 ^ x - 1)
      )
    }
    genes <- dplyr::rename(x$genes, hgnc_symbol = gene)
    x <- job_limma(new_dge(x$metadata, counts, genes))
    x$project <- project
    return(x)
  })
