# ==========================================================================
# workflow of estimate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_estimate <- setClass("job_estimate", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "estimate",
    info = c(""),
    cite = "[@Inferring_tumou_Yoshih_2013]",
    method = "",
    tag = "",
    analysis = "estimate 免疫评分"
    ))

setGeneric("asjob_estimate", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_estimate"))

setMethod("asjob_estimate", 
  signature = c(x = "job_limma"),
  function(x, use = if (x$isTcga) "gene_name" else "hgnc_symbol", dir = "estimate")
  {
    if (x@step < 1) {
      stop("The data should be normalized.")
    }
    object <- x$normed_data
    rownames(object) <- object$genes[[ use ]]
    data <- object$E
    x <- job_estimate(data)
    x <- snapAdd(x, "以数据集 ({x$project}, dataset: {x@sig}) 进行 ESTIMATE 免疫评分计算。")
    x
  })

job_estimate <- function(data, dir = "estimate")
{
  if (!dir.exists(dir)) {
    dir.create(dir)
  }
  inputfile <- file.path(dir, "input.tsv")
  if (any(duplicated(rownames(data)))) {
    message("`data` has duplicated rownames, deduplicated herein.")
    data <- data[!duplicated(rownames(data)), ]
  }
  meta <- tibble::tibble(sample = colnames(data), mutate = make.names(colnames(data)))
  data.table::fwrite(data, inputfile, row.names = T, sep = "\t")
  .job_estimate(object = inputfile, params = list(dir = dir, meta = meta))
}

setMethod("step0", signature = c(x = "job_estimate"),
  function(x){
    step_message("Prepare your data with function `job_estimate`.")
  })

setMethod("step1", signature = c(x = "job_estimate"),
  function(x, platform = c("illumina", "affymetrix", "agilent")){
    step_message("Estimate score.")
    require(estimate)
    platform <- match.arg(platform)
    file_filter <- file.path(x$dir, "filter.tsv")
    e(estimate::filterCommonGenes(object(x), file_filter, id = "GeneSymbol"))
    x$file_score <- file_score <- file.path(x$dir, "score.tsv")
    e(estimate::estimateScore(file_filter, file_score, platform))
    return(x)
  })

setMethod("step2", signature = c(x = "job_estimate"),
  function(x, metadata = NULL, sig.test = NULL)
  {
    step_message("Collate results")
    t.immuneScores <- ftibble(x$file_score, skip = 2)
    t.immuneScores <- dplyr::select(t.immuneScores, -Description)
    t.immuneScores <- tidyr::pivot_longer(t.immuneScores, -NAME,
      names_to = "mutate", values_to = "score")
    t.immuneScores <- map(t.immuneScores, "mutate", x$meta, "mutate", "sample")
    t.immuneScores <- dplyr::group_by(t.immuneScores, NAME)
    t.immuneScores <- dplyr::mutate(t.immuneScores,
      Group = ifelse(score > median(score), "High", "Low")
    )
    t.immuneScores <- dplyr::ungroup(t.immuneScores)
    if (!is.null(metadata) && !is.null(sig.test)) {
      t.immuneScores <- map(t.immuneScores, "sample", metadata, "sample", sig.test, col = sig.test)
      p.immuneScoresPlot <- wrap(.map_boxplot2(t.immuneScores, T,
        x = "Group", y = sig.test, xlab = "Group", ylab = sig.test, ids = "NAME"
      ), 8, 3.5)
    } else {
      p.immuneScores <- NULL
    }
    x <- tablesAdd(x, t.immuneScores)
    x <- plotsAdd(x, p.immuneScoresPlot)
    x <- methodAdd(x, "以 R 包 `estimate` ({packageVersion('estimate')}) {cite_show('Inferring_tumou_Yoshih_2013')} 预测数据集的 stromal, immune, estimate 得分。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_estimate"),
  function(x, metadata, group = "group",
    file_tisidb = .prefix(file.path("TISIDB", "immunomodulator.txt"), "db"))
  {
    step_message("Use genes in TISIDB...")
    data <- dplyr::rename(ftibble(object(x)), gene = 1)
    res <- .plot_imMod_boxplot(data, file_tisidb, metadata, group)
    p.Top10ImmuneRelatedGenes <- res$p.Top10ImmuneRelatedGenes
    t.SignificantImmuneRelatedGenes <- res$dataSig
    x <- plotsAdd(x, p.Top10ImmuneRelatedGenes)
    x <- tablesAdd(x, t.SignificantImmuneRelatedGenes)
    x <- methodAdd(x, "从 TISIDB {cite_show('TISIDB_an_inte_Ru_Be_2019')} 数据库下载的 178 个基因 (genes encoding immunomodulators and chemokines) 比较表达量差异。")
    return(x)
  })

.plot_imMod_boxplot <- function(data, file_tisidb, metadata, group = "group") {
  if (!file.exists(file_tisidb)) {
    dir.create(dirname(file_tisidb), F)
    immunoModulator <- data.table::fread(
      text = RCurl::getURL("http://cis.hku.hk/TISIDB/data/immunomodulator.txt"), header = F)
    data.table::fwrite(immunoModulator, file_tisidb, sep = "\t")
  } else {
    immunoModulator <- ftibble(file_tisidb)
  }
  immunoModulator <- dplyr::select(immunoModulator, Type = V2, Gene = V3)
  data <- dplyr::filter(data, gene %in% !!immunoModulator$Gene)
  genesFound <- nrow(data)
  data <- tidyr::pivot_longer(data, -gene, names_to = "sample", values_to = "expr")
  data <- map(data, "gene", immunoModulator, "Gene", "Type", col = "Type")
  data <- map(data, "sample", metadata, "sample", group, col = group)
  dataSig <- tidyr::pivot_wider(data, names_from = !!rlang::sym(group), values_from = expr)
  dataSig <- dplyr::group_by(dataSig, gene)
  fun <- function(x, y) wilcox.test(x, y)$p.value
  dataSig <- dplyr::summarise(dataSig, p.value = fun(High, Low))
  dataSig <- dplyr::mutate(dataSig, p.adjust = p.adjust(p.value, "fdr"))
  dataSig <- dplyr::arrange(dataSig, p.value)
  dataSig <- dplyr::filter(dataSig, p.value < .05)
  p.Top10ImmuneRelatedGenes <- .map_boxplot2(
    dplyr::filter(data, gene %in% head(dataSig$gene, n = 10)),
    T, x = group, xlab = group, y = "expr", ylab = "Expression", ids = "gene", nrow = 1
  )
  p.Top10ImmuneRelatedGenes <- wrap(p.Top10ImmuneRelatedGenes, 10, 6)
  namel(dataSig, p.Top10ImmuneRelatedGenes)
}
