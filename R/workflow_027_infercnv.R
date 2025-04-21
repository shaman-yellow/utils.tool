# ==========================================================================
# workflow of infercnv
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_infercnv <- setClass("job_infercnv",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://bioconductor.org/packages/release/bioc/html/infercnv.html"),
    method = "Package inferCNV used for CNV anlysis and cancer cell prediction",
    tag = "scrna:cancer",
    analysis = "InferCNV 变异拷贝数分析"
    ))

setGeneric("asjob_infercnv", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_infercnv"))

setMethod("asjob_infercnv", signature = c(x = "job_seurat"),
  function(x, ref, groups = NULL, subset = "seurat_clusters",
    group.by = x$group.by, outdir = glue::glue("infercnv_{x@sig}"))
  {
    if (missing(ref)) {
      stop('missing(ref).')
    }
    metadata <- dplyr::select(
      as_tibble(object(x)@meta.data), rownames, !!rlang::sym(group.by), !!rlang::sym(subset)
    )
    if (any(!ref %in% metadata[[ group.by ]])) {
      stop('any(!ref %in% metadata[[ group.by ]]).')
    }
    counts <- SeuratObject::LayerData(object(x), "count")
    rownames(counts) <- gname(rownames(counts))
    counts <- counts[ !duplicated(rownames(counts)), ]
    if (!is.null(groups)) {
      if (any(!ref %in% groups)) {
        allGroups <- unique(c(groups, ref))
      }
      message(glue::glue("Before Cells filter: {bind(dim(counts))}"))
      metadata <- dplyr::filter(metadata, !!rlang::sym(group.by) %in% allGroups)
      counts <- counts[, colnames(counts) %in% metadata$rownames]
      message(glue::glue("After Cells filter: {bind(dim(counts))}"))
    }
    if (!is.null(subset)) {
      metadata <- dplyr::mutate(
        metadata, group_subset = paste0(
          !!rlang::sym(group.by), "_", !!rlang::sym(subset)
          ),
        group_subset = ifelse(
          !!rlang::sym(group.by) %in% !!ref, 
          as.character(!!rlang::sym(group.by)), group_subset
        )
      )
      message(glue::glue("\n{showStrings(metadata$group_subset, trunc = FALSE)}"))
      metadata <- dplyr::select(metadata, rownames, group_subset)
      group.by <- "group_subset"
    }
    ranges <- get_gene_ranges(rownames(counts))
    counts <- counts[rownames(counts) %in% ranges$symbols, ]
    ranges <- ranges[match(rownames(counts), ranges$symbols), ]
    genes <- dplyr::select(
      ranges, symbols, seqnames, start, end
    )
    dir.create(outdir, FALSE)
    tmp.metadata <- file.path(outdir, "metadata.tsv")
    write_tsv(metadata, tmp.metadata, col.names = FALSE)
    tmp.genes <- file.path(outdir, "genes.tsv")
    write_tsv(genes, tmp.genes, col.names = FALSE)
    cell_groups <- unique(as.character(metadata[[ group.by ]]))
    obj.cnv <- e(infercnv::CreateInfercnvObject(counts,
        tmp.genes, tmp.metadata, ref))
    x <- .job_infercnv(object = obj.cnv)
    x$outdir <- outdir
    x$tmp.metadata <- tmp.metadata
    x$tmp.genes <- tmp.genes
    x <- methodAdd(x, "以 R 包 `infercnv` ({packageVersion('infercnv')}) 探索肿瘤单细胞 RNA 测序数据，以识别体细胞大规模染色体拷贝数的变异。")
    x <- snapAdd(x, "使用 `inferCNV` 识别肿瘤细胞的染色体拷贝数变异，选择 {bind(ref)} 为参考细胞 (正常细胞)，识别 {bind(groups)} 中的拷贝数变异。分析中，{bind(groups)} 以 {subset} 标识次级聚类。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Prepare your data with function `asjob_infercnv`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_infercnv"),
  function(x, workers = 4, cutoff = .1, hmm = FALSE, ...){
    step_message("Run inferCNV.")
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L, ...,
        {
          x <- step1(
            x, workers = "{workers}",
            cutoff = "{cutoff}", hmm = "{hmm}"
          )
        }
      )
      return(x)
    }
    if (dir.exists(x$outdir)) {
      unlink(x$outdir, TRUE, TRUE)
    }
    options(scipen = 100)
    object(x) <- e(infercnv::run(
        object(x),
        # cutoff = 1 works well for Smart-seq2
        # and cutoff = 0.1 works well for 10x Genomics
        num_threads = workers, cluster_by_groups = TRUE,
        cutoff = cutoff, denoise = TRUE, HMM = hmm,
        out_dir = x$outdir, save_rds = FALSE, save_final_rds = FALSE
        ))
    return(x)
  })

setMethod("step2", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Got results.")
    if (is.remote(x)) {
      dir <- file.path(x$map_local, x$outdir)
      if (!dir.exists(dir)) {
        get_file_from_remote(
          x$outdir, x$wd, x$map_local, recursive = TRUE
        )
      }
    } else {
      dir <- x$outdir
    }
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_infercnv"),
  function(x, wd = glue::glue("~/infercnv_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })


get_gene_ranges <- function(symbols, version = c("hg38", 
  "hg19"), gname = TRUE)
{
  raw <- symbols
  if (gname) {
    symbols <- gname(symbols)
  }
  version <- match.arg(version)
  name_fun <- glue::glue("TxDb.Hsapiens.UCSC.{version}.knownGene")
  if (!requireNamespace(name_fun, quietly = TRUE)) {
    BiocManager::install(name_fun)
  }
  db <- get_fun(name_fun, asNamespace(name_fun), "S4")
  ranges <- e(GenomicFeatures::genes(db))
  # if (!requireNamespace("EnsDb.Hsapiens.v86")) {
  #   BiocManager::install("EnsDb.Hsapiens.v86")
  #   # EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  # }
  entrez <- e(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "ALIAS",
    keys = symbols,
    column = c("ENTREZID")
  ))
  hasThats <- !is.na(entrez) & (entrez %in% ranges$gene_id)
  entrez <- entrez[ hasThats ]
  ranges <- ranges[ match(entrez, ranges$gene_id) ]
  ranges <- data.frame(ranges)
  ranges$symbols <- raw[ hasThats ]
  return(ranges)
}

