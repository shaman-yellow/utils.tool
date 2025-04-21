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
  function(x, ref.pattern = "Macrophage", group.by = x$group.by, outdir = "infercnv")
  {
    if (is.null(x$mart)) {
      mart <- new_biomart()
    } else {
      mart <- x$mart
    }
    counts <- SeuratObject::LayerData(object(x))
    rownames(counts) <- gs(rownames(counts), "\\.[0-9]*", "")
    genes <- filter_biomart(mart, general_attrs(), "hgnc_symbol", rownames(counts))
    genes <- dplyr::select(genes,
      hgnc_symbol, chromosome_name, start_position, end_position
    )
    counts <- counts[rownames(counts) %in% genes$hgnc_symbol, ]
    metadata <- dplyr::select(as_tibble(object(x)@meta.data), rownames, !!rlang::sym(group.by))
    # x$infercnv_used <- namel(metadata, counts, genes)
    tmp.metadata <- tempfile("cnv_metadata", fileext = ".tsv")
    write_tsv(metadata, tmp.metadata, col.names = FALSE)
    tmp.genes <- tempfile("cnv_metadata", fileext = ".tsv")
    write_tsv(genes, tmp.genes, col.names = FALSE)
    cells <- metadata[[ group.by ]]
    ref <- cells[ grepl(ref.pattern, cells, ignore.case = TRUE) ]
    obj.cnv <- e(infercnv::CreateInfercnvObject(as.matrix(counts),
        tmp.genes, tmp.metadata, ref, chr_exclude = c("X", "Y", "M")))
    x <- .job_infercnv(object = obj.cnv)
    x$outdir <- outdir
    x$tmp.metadata <- tmp.metadata
    x$tmp.genes <- tmp.genes
    return(x)
  })

setMethod("step0", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Prepare your data with function `asjob_infercnv`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_infercnv"),
  function(x, hmm = FALSE, cutoff = .1){
    step_message("Run inferCNV.")
    unlink(x$outdir, TRUE, TRUE)
    object(x) <- e(infercnv::run(
        object(x),
        # cutoff = 1 works well for Smart-seq2
        # and cutoff = 0.1 works well for 10x Genomics
        cutoff = cutoff, denoise = TRUE, HMM = hmm, plot_steps = FALSE,
        no_prelim_plot = TRUE, out_dir = x$outdir, output_format = "pdf",
        save_rds = FALSE, save_final_rds = FALSE
        ))
    return(x)
  })
