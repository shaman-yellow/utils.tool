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
    info = c("Tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html")
    ))

job_wgcna <- function(metadata, log_counts,
  gene_annotation, rename_id = "hgnc_symbol")
{
  elist <- new_elist(metadata, log_counts, gene_annotation, rename_id)
  datExpr0 <- as_wgcData(elist)
  .job_wgcna(object = elist, params = list(datExpr0 = datExpr0))
}

setMethod("step0", signature = c(x = "job_wgcna"),
  function(x){
    callNextMethod()
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
    x <- callNextMethod()
    step_message("Cluster sample tree.",
    "This do:",
    "generate `x@params$raw_sample_tree`; `x@plots[[ 1 ]]`"
    )
    raw_sample_tree <- draw_sampletree(params(x)$datExpr0)
    raw_sample_tree.p <- wrap(recordPlot(), 14, 7)
    x@params$raw_sample_tree <- raw_sample_tree
    x@plots[[ 1 ]] <- list(raw_sample_tree = raw_sample_tree.p)
    return(x)
  })

setMethod("step2", signature = c(x = "job_wgcna"),
  function(x, size, height){
    x <- callNextMethod(x)
    step_message("Cut sample tree with `height` and `size`. ",
      "This do: ",
      "clip `x@object`; generate `x@params$datExpr`; ",
      "generate `x@params$allTraits`. "
    )
    if (missing(size) | missing(height))
      stop("missing(size) | missing(height)")
    datExpr <- exclude(params(x)$datExpr0,
      cut_tree(params(x)$raw_sample_tree, size, height))
    x@params$datExpr <- datExpr
    object(x) <- clip_data(object(x), datExpr)
    x@params$allTraits <- as_wgcTrait(object(x))
    return(x)
  })

setMethod("step3", signature = c(x = "job_wgcna"),
  function(x){
    x <- callNextMethod(x)
    step_message("Analysis of network topology for soft-thresholding powers. ",
      "This do: ",
      "Generate x@params$sft; plots in `x@plots[[ 3 ]]`. "
    )
    sft <- cal_sft(params(x)$datExpr)
    x@params$sft <- sft
    x@plots[[ 3 ]] <- list(sft = wrap(plot_sft(sft), 10, 5))
    return(x)
  })

setMethod("step4", signature = c(x = "job_wgcna"),
  function(x, power = x@params$sft$powerEstimate, ...){
    x <- callNextMethod(x)
    step_message("One-step network construction and module detection.
      Extra parameters would passed to `cal_module`.
      This do: Generate `x@params$MEs`; plots (net) in `x@plots[[ 4 ]]`.
      By default, red{{x@params$sft$powerEstimate}} is used
      as `power` for WGCNA calculation.
      "
    )
    net <- cal_module(params(x)$datExpr, power, ...)
    x@plots[[ 4 ]] <- list(net = net)
    x@params$MEs <- get_eigens(net)
    return(x)
  })

setMethod("step5", signature = c(x = "job_wgcna"),
  function(x){
    x <- callNextMethod(x)
    step_message("Correlation test for modules with trait data. ",
      "This do:",
      "Generate plots in `x@plots[[ 5 ]]`; ",
      "tables in `x@tables[[ 5 ]]`"
    )
    hps_corp <- new_heatdata(params(x)$MEs, params(x)$allTraits)
    hps_corp <- callheatmap(hps_corp)
    x@plots[[ 5 ]] <- list(hps_corp = hps_corp)
    x@tables[[ 5 ]] <- list(corp = hps_corp@data_long)
    return(x)
  })

setMethod("step6", signature = c(x = "job_wgcna"),
  function(x){
    x <- callNextMethod(x)
    step_message("Calculate gene significance (GS) and module membership (MM).",
      "This do:",
      "Generate `x@params$mm`, `x@params$gs`; ",
      "tables (filter by pvalue < 0.05) `x@tables[[ 6 ]]`"
    )
    mm <- cal_corp(params(x)$datExpr, params(x)$MEs, "gene", "module")
    mm.s <- dplyr::filter(tibble::as_tibble(mm), pvalue < .05)
    gs <- cal_corp(params(x)$datExpr, params(x)$allTraits, "gene", "trait")
    gs.s <- dplyr::mutate(tibble::as_tibble(gs), p.adjust = p.adjust(pvalue, "BH"))
    gs.s <- dplyr::filter(gs.s, pvalue < .05)
    x@params$mm <- mm
    x@params$gs <- gs
    x@tables[[ 6 ]] <- list(mm = mm.s, gs = gs.s)
    return(x)
  })

