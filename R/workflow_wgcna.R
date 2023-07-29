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
  prototype = NULL)

job_wgcna <- function(metadata, log_counts,
  gene_annotation, rename_id = "hgnc_symbol")
{
  elist <- new_elist(metadata, log_counts, gene_annotation, rename_id)
  datExpr0 <- as_wgcData(elist)
  .job_wgcna(object = elist, params = list(datExpr0 = datExpr0))
}

setMethod("step0", signature = c(x = "job_wgcna"),
  function(x){
    step0()
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
    step_message("Cluster sample tree.")
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
    datExpr <- exclude(params(x)$datExpr0, cut_tree(tree, size, height))
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
    step_message("One-step network construction and module detection. ",
      "Extra parameters would passed to `cal_module`. ",
      "Generate `x@params$MEs`; plots (net) in `x@plots[[ 4 ]]`. "
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
      "Generate plots in `x@plots[[ 4 ]]`. "
    )
    hps_corp <- new_heatdata(params(x)$MEs, params(x)$allTraits)
    hps_corp <- callheatmap(hps_corp)
    x@plots[[ 5 ]] <- list(hps_corp = hps_corp)
    return(x)
  })
