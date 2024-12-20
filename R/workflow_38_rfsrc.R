# ==========================================================================
# workflow of rfsrc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_rfsrc <- setClass("job_rfsrc", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "",
    method = "R package `randomForestSRC` for feature selection",
    tag = "fea",
    analysis = "RandomForestSRC 随机森林"
    ))

job_rfsrc <- function()
{
  .job_rfsrc()
}

setGeneric("asjob_rfsrc", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_rfsrc"))

setMethod("asjob_rfsrc", signature = c(x = "job_maf"),
  function(x){
    data <- e(maftools::genesToBarcodes(object(x), object(x)@gene.summary[[1]], justNames = T))
    data <- as_tibble(as_df.lst(data))
    data <- dplyr::mutate(data, value = 1L)
    data <- e(tidyr::spread(data, type, value))
    rownames <- data$name
    data <- map(data, "name", x$clinical, "Tumor_Sample_Barcode", "vital_status")
    data <- e(dplyr::mutate_all(data, function(x) ifelse(is.na(x), 0L, x)))
    data <- as.data.frame(data)
    rownames(data) <- rownames
    data <- e(dplyr::mutate(data, vital_status = as.factor(vital_status)))
    data <- dplyr::rename(data, class = vital_status)
    .job_rfsrc(object = data)
  })

setMethod("step0", signature = c(x = "job_rfsrc"),
  function(x){
    step_message("Prepare your data with function `job_rfsrc`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_rfsrc"),
  function(x, top = 30){
    step_message("Classification importance.")
    if (is.null(x$imp)) {
      model <- e(randomForestSRC::rfsrc(class ~ ., data = object(x), splitrule = "gini"))
      imp <- e(randomForestSRC::vimp(model))
      x$imp <- imp
      x$model <- model
    } else {
      imp <- x$imp
    }
    t.imp <- as_tibble(imp$importance)
    t.imp <- dplyr::arrange(t.imp, dplyr::desc(all))
    t.imp <- .set_lab(t.imp, sig(x), "Variable importance")
    data <- dplyr::slice(as_tibble(t.imp), 1:top)
    p.imp <- ggplot(data) +
      geom_point(aes(x = reorder(rownames, all), y = all), shape = 21) +
      coord_flip() +
      labs(y = "Variable importance (Gini)", x = "") +
      ggthemes::theme_calc()
    p.imp <- .set_lab(wrap(p.imp), sig(x), "Variable importance")
    x@tables[[ 1 ]] <- namel(t.imp)
    x@plots[[ 1 ]] <- namel(p.imp)
    return(x)
  })
