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

job_cardinal <- function(files, names = NULL, ...)
{
  if ( !is.null(names) ) {
    if ( length(files) != length(names) ) {
      stop("files and the names must be the same length.")
    }
  }
  n <- 0L
  objects <- pbapply::pblapply(files,
    function(file) {
      n <<- n + 1L
      object <- Cardinal::readMSIData(file, ...)
      if ( !is.null(names) ) {
        Cardinal::pixelData(object)$run <- factor(names[ n ])
      }
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
    x@plots[[ 1 ]] <- plots
    x@tables[[ 1 ]] <- namel(t.features)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cardinal"),
  function(x){
    step_message("Statistic.")
    
    return(x)
  })

setMethod("focus", signature = c(x = "job_cardinal"),
  function(x, ...)
  {
    Cardinal::plot(object(x), log2(intensity + 1) ~ mz,
      ylab = "Log2(intensity + 1)", zlab = "m/z", ...)
  })

setMethod("vis", signature = c(x = "job_cardinal"),
  function(x, mz = NULL, layout = c(1, length(ids(x))), ...)
  {
    wrap(Cardinal::image(object(x), mz = mz, layout = layout, ...),
      1.5 + 1.5 * layout[2], 2.5 * layout[1])
  })

setMethod("ids", signature = c(x = "job_cardinal"),
  function(x, which = "run", unique = T){
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
