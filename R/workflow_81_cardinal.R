# ==========================================================================
# workflow of cardinal
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# http://maldi-msi.org/index.php?option=com_content&view=article&id=189&Itemid=69#ImzMLConverter

.job_cardinal <- setClass("job_cardinal", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "cardinal",
    info = c("https://cardinalmsi.org/"),
    cite = "[@CardinalV3ABemis2023]",
    method = "The R package `Cardinal` used for analyzing mass spectrometry imaging datasets",
    tag = "cardinal",
    analysis = "Cardinal 空间代谢组数据分析"
    ))

job_cardinal <- function(files, names = NULL, groups = NULL, ...)
{
  if ( !is.null(names) ) {
    if ( length(files) != length(names) ) {
      stop("files and the names must be the same length.")
    }
  }
  if ( is.null(groups) ) {
    stop("`groups` can not be NULL.")
  }
  n <- 0L
  objects <- pbapply::pblapply(files,
    function(file) {
      n <<- n + 1L
      object <- Cardinal::readMSIData(file, ...)
      if ( !is.null(names) ) {
        Cardinal::pixelData(object)$run <- factor(names[ n ])
      }
      Cardinal::pixelData(object)$group <- factor(groups[ n ])
      object
    })
  if (length(objects) > 1) {
    object <- do.call(c, objects)
  } else {
    object <- objects[[ 1 ]]
  }
  .job_cardinal(object = object)
}

setMethod("step0", signature = c(x = "job_cardinal"),
  function(x){
    step_message("Prepare your data with function `job_cardinal`.")
  })

setMethod("step1", signature = c(x = "job_cardinal"),
  function(x, test = 1:6)
  {
    step_message("Preprocessing.")
    object(x) <- e(Cardinal::normalize(object(x)))
    p.normalize <- focus(x, i = test)
    object(x) <- e(Cardinal::smooth(object(x)))
    p.smooth <- focus(x, i = test)
    object(x) <- e(Cardinal::reduceBaseline(object(x)))
    p.reduceBaseline <- focus(x, i = test)
    object(x) <- e(Cardinal::peakPick(object(x), SNR = 2))
    p.peakPick <- focus(x, i = test)
    object(x) <- e(Cardinal::peakAlign(object(x)))
    plots <- namel(p.normalize, p.smooth, p.reduceBaseline, p.peakPick)
    plots <- lapply(plots,
      function(x) {
        wrap(x, 10, 5)
      })
    t.features <- as_tibble(data.frame(Cardinal::featureData(x@object)))
    t.features <- dplyr::mutate(t.features, i = seq_len(nrow(t.features)), .before = mz)
    t.features <- .set_lab(t.features, sig(x), "all features")
    x@plots[[ 1 ]] <- plots
    x@tables[[ 1 ]] <- namel(t.features)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cardinal"),
  function(x){
    step_message("Statistic (PCA).")
    x$pca <- e(Cardinal::PCA(object(x)))
    p.pca <- Cardinal::plot(x$pca, type = "x", groups = x$pca@pixelData$group)
    p.pca <- wrap(p.pca, 8, 6)
    p.pca <- .set_lab(p.pca, sig(x), "PCA plot")
    x@plots[[ 2 ]] <- namel(p.pca)
    return(x)
  })

setMethod("step3", signature = c(x = "job_cardinal"),
  function(x, formula = "~ group", show = 10,
    use.segment = FALSE, workers = 5, ...)
  {
    step_message("Segment-based means testing.")
    if ( use.segment ) {
      object <- Cardinal::spatialDGMM(object(x), groups = Cardinal::run(object(x)),
        BPPARAM = BiocParallel::MulticoreParam(workers), ...)
    } else {
      object <- object(x)
    }
    x$mtest <- Cardinal::meansTest(object, as.formula(formula), samples = Cardinal::run(object(x)))
    t.tops <- as_tibble(data.frame(Cardinal::topFeatures(x$mtest)))
    t.tops <- dplyr::filter(t.tops, fdr < .05)
    t.tops <- .set_lab(t.tops, sig(x), "Significant differences features")
    x@tables[[ 3 ]] <- namel(t.tops)
    get_data.mtest <- function(x, i) {
      x <- lapply(x[ i ],
        function(obj) {
          obj$data
        })
      names(x) <- paste0("I = ", i)
      x <- frbind(x, idcol = "features")
      x <- dplyr::mutate(x, features = sortCh(features, FALSE, TRUE))
      x
    }
    data <- get_data.mtest(x$mtest, tops <- head(t.tops$i, show))
    p.boxTops <- .map_boxplot2(data, FALSE, x = "group", y = "intensity",
      ids = "features", scales = "free_y",
      xlab = "Group", ylab = "Normalized Itensity"
    )
    p.boxTops <- wrap(p.boxTops)
    p.boxTops$tops <- tops
    p.boxTops <- .set_lab(p.boxTops, sig(x), "boxplot of top features")
    x@plots[[ 3 ]] <- namel(p.boxTops)
    return(x)
  })

setMethod("map", signature = c(x = "job_cardinal", ref = "df"),
  function(x, ref, order_by = "i", distinct = TRUE, tol = .01, ...)
  {
    message("Merge identification data into feature data and top tables.")
    if (x@step < 1L) {
      stop("The value got TRUE: `x@step < 1L`")
    }
    fun_merge <- function(x, ref, distinct, ...) {
      x <- as_tibble(tol_merge(x, ref,
          main_col = "mz", sub_col = "mz", tol = tol, ...
          ))
      x <- dplyr::arrange(x, !!!rlang::syms(as.list(order_by)))
      if (distinct) {
        x <- dplyr::distinct(x, i, .keep_all = TRUE)
      }
      x
    }
    if (!is.null(x$t.featuresRaw)) {
      x@tables$step1$t.features <- x$t.featuresRaw
    } else {
      x$t.featuresRaw <- x@tables$step1$t.features
    }
    x@tables$step1$t.features <- fun_merge(x@tables$step1$t.features, ref, distinct, ...)
    if (x@step >= 3L) {
      if (!is.null(x$t.topsRaw)) {
        x@tables$step3$t.tops <- x$t.topsRaw
      } else {
        x$t.topsRaw <- x@tables$step3$t.tops
      }
      x@tables$step3$t.tops <- fun_merge(x@tables$step3$t.tops, ref, distinct, ...)
    }
    return(x)
  })

setMethod("focus", signature = c(x = "job_cardinal"),
  function(x, ...)
  {
    Cardinal::plot(object(x), log2(intensity + 1) ~ mz,
      ylab = "Log2(intensity + 1)", zlab = "m/z", ...)
  })

setMethod("vis", signature = c(x = "job_cardinal"),
  function(x, i = NULL, mz = NULL, layout = c(1, length(ids(x))), ids = NULL, use.ids = NULL, ...)
  {
    if (!is.null(ids) && !is.null(use.ids)) {
      meta <- dplyr::filter(x@tables$step1$t.features, !!rlang::sym(use.ids) %in% ids)
      meta <- .set_lab(meta, sig(x), "metadata-of-visualized-metabolites")
      p.lst <- lapply(meta$i,
        function(i) {
          wrap(Cardinal::image(object(x), i = i, layout = layout, ...),
            1.5 + 1.5 * layout[2], 2.5 * layout[1])
        })
      names(p.lst) <- paste0("Feature_", meta$i, "_", meta[[ use.ids ]])
      p.lst <- .set_lab(p.lst, sig(x), names(p.lst), "image-visualization")
      lab(p.lst) <- "Feature image visualizations"
      return(namel(p.lst, meta))
    }
    p <- wrap(Cardinal::image(object(x), i = i, mz = mz, layout = layout, ...),
      1.5 + 1.5 * layout[2], 2.5 * layout[1])
    p <- .set_lab(p, sig(x), paste0("image plot of Feature ", i))
    p
  })

setMethod("ids", signature = c(x = "job_cardinal"),
  function(x, which = "run", unique = TRUE){
    ids <- object(x)[[ which ]]
    if ( unique ) {
      if ( is.factor(ids) ) {
        ids <- levels(ids)
      } else {
        ids <- unique(ids)
      }
    }
    ids
  })
