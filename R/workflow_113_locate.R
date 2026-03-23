# ==========================================================================
# workflow of locate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_locate <- setClass("job_locate", 
  contains = c("job"),
  prototype = prototype(
    pg = "locate",
    info = c("https://mrcuizhe.github.io/interacCircos_documentation/html/users_from_rcircos.html",
      "http://www.rnalocate.org/"),
    cite = "",
    method = "",
    tag = "locate",
    analysis = "染色体定位和亚细胞定位"
    ))

setGeneric("asjob_locate",
  function(x, ...) standardGeneric("asjob_locate"))

setMethod("asjob_locate", signature = c(x = "feature"),
  function(x, ...){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_locate(object = fea)
    dir.create("tmp", FALSE)
    x$gene_data <- expect_local_data(
      "tmp", "biomart", get_gene_positions, list(genes = fea), ...
    )
    return(x)
  })

setMethod("step0", signature = c(x = "job_locate"),
  function(x){
    step_message("Prepare your data with function `job_locate`.")
  })

setMethod("step1", signature = c(x = "job_locate"),
  function(x){
    step_message("Quality control (QC).")
    p.locateChr <- funPlot(plot_genes_in_RCircos, list(gene_data = x$gene_data))
    p.locateChr <- set_lab_legend(
      p.locateChr,
      glue::glue("{x@sig} Chromosome localization"),
      glue::glue("基因于染色体定位|||外圈数字表示染色体（1-22表示1-22条人类染色体，XY对应性染色体）")
    )
    snap <- glue::glue(
      "基因 {x$gene_data$Gene} 位于 {s(x$gene_data$Chromosome, 'chr', '')} 染色体上"
    )
    x <- snapAdd(x, "如图所示{aref(p.locateChr)}，{bind(snap)}。")
    x <- methodAdd(x, "染色体定位分析可有效揭示基因在染色体上的分布特征。以 R 包 `RCircos` ({packageVersion('RCircos')}) 生成基因染色体定位图谱。")
    x <- plotsAdd(x, p.locateChr)
    return(x)
  })

setMethod("step2", signature = c(x = "job_locate"),
  function(x, org = "Homo sapiens"){
    step_message("Get location data.")
    data <- get_mRNA_subcellular_data(org = org)
    data <- dplyr::filter(data, RNA_Symbol %in% object(x))
    if (any(whichNot <- !object(x) %in% data$RNA_Symbol)) {
      message(glue::glue("Not got: {object(x)[ whichNot ]}"))
    }
    data <- dplyr::arrange(data, dplyr::desc(RNALocate_Score))
    data <- dplyr::distinct(
      data, RNA_Symbol, Subcellular_Localization, .keep_all = TRUE
    )
    x$locateData <- data
    p.locateScore <- wrap_scale(
      .plot_subcellular_scores(data), 
      length(unique(data$Subcellular_Localization)), 
      length(object(x)[!whichNot]), 
      pre_height = 3.5, min_width = 2
    )
    x <- methodAdd(x, "从 RNALocate v3.0 (<http://www.rnalocate.org/>) 获取 mRNA 亚细胞定位数据，并用 R 将定位和得分数据可视化。")
    p.locateScore <- set_lab_legend(
      p.locateScore,
      glue::glue("{x@sig} RNA Subcellular Localization Distribution"),
      glue::glue("RNA 亚细胞定位分布|||纵坐标为不同基因，横坐标为的预测的蛋白质亚细胞定位得分：RNA 亚细胞定位关联信息来自不同类型的资源，包括实验证据和预测证据；实验证据对置信度评分的贡献应该比预测证据更大；强有力的实验证据应该比薄弱的实验证据提供更可靠的证据；有更多证据支持的 RNA 亚细胞定位关联应比证据较少支持的关联具有更高的置信度评分 (<http://www.rnalocate.org/help>)。")
    )
    top <- dplyr::distinct(data, RNA_Symbol, .keep_all = TRUE)
    snap <- glue::glue("蛋白 {top$RNA_Symbol} 分布于 {top$Subcellular_Localization}")
    x <- snapAdd(x, "如图{aref(p.locateScore)}，最有证据证明 {bind(snap)}。")
    x <- plotsAdd(x, p.locateScore)
    return(x)
  })

get_mRNA_subcellular_data <- function(org = "Homo sapiens", dir = .prefix("mRNA_subcellular", "db"))
{
  filename <- "mRNA subcellular localization information.txt"
  file <- file.path(dir, filename)
  if (!file.exists(file)) {
    dir.create(dir, FALSE)
    url <- "http://www.rnalocate.org/static/download/mRNA%20subcellular%20localization%20information.zip"
    zipfile <- file.path(dir, "mRNA.zip")
    utils::download.file(url, zipfile)
    unzip(zipfile, exdir = dir)
  }
  dplyr::filter(ftibble(file), Species == !!org)
}

.plot_subcellular_scores <- function(data) {
  p <- ggplot(data, aes(x = Subcellular_Localization, 
      y = reorder(RNA_Symbol, RNALocate_Score), 
      size = RNALocate_Score,
      color = RNALocate_Score)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Subcellular Localization", 
    y = "RNA",
    size = "Score",
    color = "Score")
  p
}

plot_genes_in_RCircos <- function(gene_data) {
  require(RCircos)
  data(UCSC.HG38.Human.CytoBandIdeogram, package = "RCircos")
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos::RCircos.Set.Core.Components(
    cyto.info, 
    chr.exclude = NULL,
    tracks.inside = 2,
    tracks.outside = 0
  )
  params <- RCircos::RCircos.Get.Plot.Parameters()
  params$text.size <- 1
  params$point.size <- 1.2
  RCircos::RCircos.Reset.Plot.Parameters(params)
  RCircos::RCircos.Set.Plot.Area()
  RCircos::RCircos.Chromosome.Ideogram.Plot()
  RCircos::RCircos.Gene.Connector.Plot(gene_data, track.num = 1, side = "in")
  RCircos::RCircos.Gene.Name.Plot(gene_data, name.col = 4, track.num = 2, side = "in")
}


get_gene_positions <- function(genes) {
  ensembl <- new_biomart()
  gene_positions <- biomaRt::getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
    filters = "hgnc_symbol",
    values = genes,
    mart = ensembl
  )
  valid_chrs <- c(as.character(1:22), "X", "Y")
  gene_positions <- gene_positions[ gene_positions$chromosome_name %in% valid_chrs, ]
  gene_positions$chromosome <- paste0("chr", gene_positions$chromosome_name)
  result <- data.frame(
    Chromosome = gene_positions$chromosome,
    Start = gene_positions$start_position,
    End = gene_positions$end_position,
    Gene = gene_positions$hgnc_symbol,
    stringsAsFactors = FALSE
  )
  chr_order <- order(gene_positions$chromosome_name)
  result <- result[chr_order, ]
  if (nrow(result) < length(genes)) {
    not_found <- setdiff(genes, result$Gene)
    message(glue::glue("Not found: {bind(not_found)}"))
  }
  return(result)
}
