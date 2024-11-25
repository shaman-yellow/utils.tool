# ==========================================================================
# workflow of mfuzz
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mfuzz <- setClass("job_mfuzz", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "mfuzz",
    info = c("https://bioconductor.org/packages/release/bioc/vignettes/Mfuzz/inst/doc/Mfuzz.pdf"),
    cite = "[@Mfuzz_a_softwa_Kumar_2007]",
    method = "",
    tag = "mfuzz",
    analysis = "Mfuzz 聚类分析"
    ))

job_mfuzz <- function()
{
  .job_mfuzz()
}

setMethod("step0", signature = c(x = "job_mfuzz"),
  function(x){
    step_message("Prepare your data with function `job_mfuzz`.")
  })

setMethod("step1", signature = c(x = "job_mfuzz"),
  function(x, order = NULL, row.names = 1) {
    step_message("Convert to `ExpressionSet`.")
    data <- object(x)
    if (is.null(order)) {
      order <- seq_along(colnames(data))[-row.names]
    }
    data <- as.matrix(data.frame(data[, order], row.names = data[[ row.names ]]))
    data <- new('ExpressionSet', exprs = data)
    object(x) <- data
    return(x)
  })

setMethod("step2", signature = c(x = "job_mfuzz"),
  function(x, centers = 8, mflow = c(2, 4), alpha = .1){
    step_message("Clustering of genes.")
    require(Mfuzz)
    m <- e(Mfuzz::mestimate(object(x)))
    clusters <- e(Mfuzz::mfuzz(object(x), c = centers, m = m))
    Mfuzz::mfuzz.plot(object(x), clusters, mflow, min.mem = alpha)
    p.clusters <- wrap(recordPlot(), 10, 6)
    p.clusters <- .set_lab(p.clusters, sig(x), "Mfuzz clusters")
    x$m <- m
    x$clusters <- clusters
    x$centers <- centers
    x@plots[[ 2 ]] <- namel(p.clusters)
    meth(x)$step2 <- glue::glue("以 R 包 `Mfuzz` ({packageVersion('Mfuzz')}) {cite_show('Mfuzz_a_softwa_Kumar_2007')} 对基因聚类分析，设定 fuzzification 参数为 {m} (以 `Mfuzz::mestimate` 预估) ，得到 {centers} 个聚类。")
    return(x)
  })

setMethod("asjob_enrich", signature = c(x = "job_mfuzz"),
  function(x, clusters, ...){
    step_message("Extract genes for selected clusters.")
    clusters <- lapply(clusters,
      function(n) {
        alls <- x$clusters$cluster
        names(alls)[alls %in% n]
      })
    job_enrich(clusters, ...)
  })

setGeneric("asjob_mfuzz", 
  function(x, ...) standardGeneric("asjob_mfuzz"))

setMethod("asjob_mfuzz", 
  signature = c(x = "job_limma"),
  function(x){
    if (x@step < 3L) {
      stop("`x@step` should not less than 3L.")
    }
    data <- x@plots$step3$p.hp@data@data
    data <- dplyr::group_by(data, group, genes)
    data <- dplyr::summarize(data, expression = mean(expression))
    data <- tidyr::pivot_wider(data, names_from = group, values_from = expression)
    message("Please modify the columns of `object(x)`.")
    message(glue::glue("Now the order is \n{paste(colnames(data)[-1], collapse = '\n')}"))
    .job_mfuzz(object = data)
  })
