# ==========================================================================
# workflow of WGCNA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_wgcna <- setClass("job_wgcna", 
  contains = c("job"),
  representation = representation(
    object = "elist",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html"),
    cite = "[@WgcnaAnRPacLangfe2008]"
    ))

setGeneric("asjob_wgcna", 
  function(x, ...) standardGeneric("asjob_wgcna"))

setMethod("asjob_wgcna", signature = c(x = "job_seurat"),
  function(x, features = NULL, cells = NULL){
    step_message("Use SeuratObject:::subset.Seurat to subset the data.")
    hasIt <- features %in% rownames(object(x))
    message("Features found:")
    print(prop.table(table(hasIt)))
    sub <- e(suppressWarnings(SeuratObject:::subset.Seurat(object(x),
        features = features[ hasIt ], cells = cells
        )))
    log_counts <- as_tibble(sub[[ SeuratObject::DefaultAssay(sub) ]]@scale.data)
    metadata <- as_tibble(sub@meta.data)
    gene_annotation <- tibble::tibble(gene = rownames(sub))
    job_wgcna(metadata, log_counts, gene_annotation, "gene")
  })

job_wgcna <- function(metadata, log_counts,
  gene_annotation)
{
  elist <- new_elist(metadata, log_counts, gene_annotation)
  datExpr0 <- as_wgcData(elist)
  .job_wgcna(object = elist, params = list(datExpr0 = datExpr0))
}

setMethod("step0", signature = c(x = "job_wgcna"),
  function(x){
    step_message("Prepare your data with function `job_wgcna`. ",
      "Note that the ", crayon::red("first column"),
      " in each data (metadata, counts, genes annotation)",
      " were used as ID column of ",
      "corresponding content. \n",
      crayon::red("metadata:"), " with 'sample' and 'group' (optional).\n",
      crayon::red("counts:"), " with id column and then expression columns.\n",
      "genes: ", crayon::red("Traits data"), " could be placed herein; ",
      "a step would serve all numeric data in `genes` as trait data then ",
      "perform statistic."
    )
  })

setMethod("step1", signature = c(x = "job_wgcna"),
  function(x){
    step_message("Cluster sample tree.",
    "This do:",
    "generate `x@params$raw_sample_tree`; `x@plots[[ 1 ]]`"
    )
    raw_sample_tree <- draw_sampletree(params(x)$datExpr0)
    raw_sample_tree.p <- wrap(recordPlot(), nrow(params(x)$datExpr0) * .8, nrow(params(x)$datExpr0))
    x@params$raw_sample_tree <- raw_sample_tree
    x@plots[[ 1 ]] <- list(raw_sample_tree = raw_sample_tree.p)
    return(x)
  })

setMethod("step2", signature = c(x = "job_wgcna"),
  function(x, size = NULL, height = NULL){
    step_message("Cut sample tree with `height` and `size`. ",
      "This do: ",
      "clip `x@object`; generate `x@params$datExpr`; ",
      "generate `x@params$allTraits`. "
    )
    if (!is.null(size)) {
      datExpr <- exclude(params(x)$datExpr0,
        cut_tree(params(x)$raw_sample_tree, size, height))
      x@params$datExpr <- datExpr
      object(x) <- clip_data(object(x), datExpr)
    }
    x@params$allTraits <- as_wgcTrait(object(x))
    return(x)
  })

setMethod("step3", signature = c(x = "job_wgcna"),
  function(x, powers = 1:50, cores = 7){
    step_message("Analysis of network topology for soft-thresholding powers. ",
      "This do: ",
      "Generate x@params$sft; plots in `x@plots[[ 3 ]]`. "
    )
    e(WGCNA::enableWGCNAThreads(cores))
    sft <- cal_sft(params(x)$datExpr, powers = powers)
    x@params$sft <- sft
    x@plots[[ 3 ]] <- list(sft = wrap(plot_sft(sft), 10, 5))
    return(x)
  })

setMethod("step4", signature = c(x = "job_wgcna"),
  function(x, power = x@params$sft$powerEstimate, cores = 7, ...){
    step_message("One-step network construction and module detection.
      Extra parameters would passed to `cal_module`.
      This do: Generate `x@params$MEs`; plots (net) in `x@plots[[ 4 ]]`.
      By default, red{{x@params$sft$powerEstimate}} is used
      as `power` for WGCNA calculation.
      "
    )
    e(WGCNA::enableWGCNAThreads(cores))
    net <- cal_module(params(x)$datExpr, power, ...)
    if (!is(net, "wgcNet"))
      net <- .wgcNet(net)
    x@plots[[ 4 ]] <- list(net = net)
    x@params$MEs <- get_eigens(net)
    return(x)
  })

setMethod("step5", signature = c(x = "job_wgcna"),
  function(x, traits = NULL){
    step_message("Correlation test for modules with trait data. ",
      "This do:",
      "Generate plots in `x@plots[[ 5 ]]`; ",
      "tables in `x@tables[[ 5 ]]`"
    )
    if (!is.null(traits)) {
      .check_columns(traits, c("sample"))
      traits <- traits[match(rownames(x@params$datExpr), traits$sample), ]
      rownames <- traits$sample
      traits <- dplyr::select_if(traits, is.numeric)
      traits <- data.frame(traits)
      rownames(traits) <- traits$sample
      x@params$allTraits <- .wgcTrait(traits)
    }
    if (is.null(params(x)$allTraits))
      stop("is.null(params(x)$allTraits) == T")
    if (ncol(params(x)$allTraits) == 0)
      stop("ncol(params(x)$allTraits) == 0, no data in `allTraits`.")
    hps_corp <- new_heatdata(params(x)$MEs, params(x)$allTraits)
    hps_corp <- callheatmap(hps_corp)
    x@plots[[ 5 ]] <- list(hps_corp = hps_corp)
    x@tables[[ 5 ]] <- list(corp = hps_corp@data_long)
    return(x)
  })

setMethod("step6", signature = c(x = "job_wgcna"),
  function(x, use.trait = NULL, use = c("adj.pvalue", "pvalue")){
    step_message("Calculate gene significance (GS) and module membership (MM).",
      "This do:",
      "Generate `x@params$mm`, `x@params$gs`; ",
      "tables (filter by pvalue < 0.05) `x@tables[[ 6 ]]`"
    )
    use <- match.arg(use)
    mm <- cal_corp(params(x)$datExpr, params(x)$MEs, "gene", "module")
    mm.s <- mutate(as_tibble(mm), adj.pvalue = p.adjust(pvalue, "BH"))
    mm.s <- dplyr::filter(mm.s, !!rlang::sym(use) < .05)
    gs <- cal_corp(params(x)$datExpr, params(x)$allTraits, "gene", "trait")
    gs.s <- dplyr::mutate(tibble::as_tibble(gs), adj.pvalue = p.adjust(pvalue, "BH"))
    gs.s <- dplyr::filter(gs.s, !!rlang::sym(use) < .05)
    if (!is.null(use.trait)) {
      gs.s <- dplyr::filter(gs.s, trait %in% dplyr::all_of(use.trait))
    }
    p.mm_gs <- new_upset(gs = gs.s$gene, mm = mm.s$gene)
    show(p.mm_gs)
    p.mm_gs <- wrap(recordPlot(), 3, 3)
    dev.off()
    x@params$mm <- mm
    x@params$gs <- gs
    x@params$ins.mm_gs <- intersect(gs.s$gene, mm.s$gene)
    x@tables[[ 6 ]] <- list(mm = mm.s, gs = gs.s)
    x@plots[[ 6 ]] <- namel(p.mm_gs)
    return(x)
  })


