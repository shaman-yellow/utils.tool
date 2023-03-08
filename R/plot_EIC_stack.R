# ==========================================================================
# plot extracted ions chromatograph (EIC) for features using `MSnbase`
# to extract mass data.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#' @aliases plot_EIC_stack
#'
#' @title Draw extracted ions chromatography for 'features'
#'
#' @description Use quantification table (with peak start time and end time)
#' exported by MCmine to draw EIC plot.
#'
#' @name plot_EIC_stack
NULL
#> NULL

#' @export plot_EIC_stack
#' @aliases plot_EIC_stack
#' @description \code{plot_EIC_stack}: ...
#' @rdname plot_EIC_stack
plot_EIC_stack <- 
  function(
    idset,
    metadata,
    quant.path,
    mzml.path,
    palette = ggsci::pal_npg()(10),
    mz.tol = 0.01,
    rt.tol = 0.1,
    cl = NULL,
    data = NULL)
  {
    if (is.null(data)) {
      .suggest_bio_package("MSnbase")
      ## metadata
      .check_columns(metadata, c("file", "sample", "group"), "metadata")
      metadata <- dplyr::arrange(metadata, sample)
      feature <- data.table::fread(quant.path)
      .check_columns(
        feature, c("row ID",	"row m/z", "row retention time"),
        "data.table::fread(quant.path)"
      )
      feature <- dplyr::select(
        feature, .features_id = 1, mz = 2, rt = 3,
        dplyr::contains(metadata$sample) & dplyr::contains("Peak RT")
      )
      feature <- dplyr::mutate(feature, .features_id = as.character(.features_id))
      feature <- dplyr::filter(feature, .features_id %in% idset)
      feature <- tidyr::gather(feature, type, time, -.features_id, -mz, -rt)
      feature <- feature.rt.during <- dplyr::mutate(
        feature, time = ifelse(time == 0, NA, time),
        sub.type = stringr::str_extract(type, "(?<=RT ).*?$"),
        sample = gsub("\\.mz.{-}ML Peak RT.*$", "", type)
      )
      feature <- dplyr::group_by(feature, .features_id, mz, rt, sub.type)
      feature <- dplyr::summarize(
        feature, sub.type.min = min(time, na.rm = T),
        sub.type.max = max(time, na.rm = T),
        .groups = "drop_last"
      )
      feature <- dplyr::mutate(
        feature, time = ifelse(sub.type == "start", sub.type.min, sub.type.max)
      )
      feature <- dplyr::select(feature, -contains("sub.type."))
      feature <- tidyr::spread(feature, sub.type, time)
      ## read data
      if (!is.null(cl))
        bioc.par(cl)
      data <- MSnbase::readMSData(
        paste0(mzml.path, "/", metadata$file),
        pdata = new("NAnnotatedDataFrame", metadata),
        mode = "onDisk"
      )
      ## extract EIC
      rt.tol.sec <- rt.tol * 60
      if (!is.null(cl))
        bioc.par(cl)
      eic.list <- pbapply::pbapply(
        feature, 1,
        function(vec){
          ## mz range for EIC
          mz <- as.numeric(vec[["mz"]])
          mz.range <- c(mz - mz.tol, mz + mz.tol)
          ## rt range for EIC
          rt.range <- c(vec[["start"]], vec[["end"]])
          rt.range <- as.numeric(rt.range) * 60
          rt.range <- c(rt.range[1] - rt.tol.sec, rt.range[2] + rt.tol.sec)
          ms1.vec <- MSnbase::chromatogram(data, msLevel = 1L, mz = mz.range,
            rt = rt.range, aggregationFun = "max")
          data.list <- lapply(unlist(ms1.vec),
            function(chr){
              int <- MSnbase::intensity(chr)
              rt <- MSnbase::rtime(chr)
              data.frame(real.time = rt, int = int)
            })
          names(data.list) <- metadata$sample
          df <- data.table::rbindlist(data.list, idcol = T)
          df <- dplyr::rename(df, sample = .id)
          dplyr::mutate(df, .features_id = vec[[".features_id"]])
        })
      ## define whether the peak belong to the feature
      eic.df <- data.table::rbindlist(eic.list)
      eic.df <- merge(
        eic.df, feature.rt.during, by = c(".features_id", "sample"), allow.cartesian = T
      )
      eic.df <- dplyr::select(eic.df, -type)
      eic.df <- tidyr::spread(eic.df, key = sub.type, value = time)
      eic.df <- merge(eic.df, metadata, by = "sample", all.x = T)
      eic.df <- dplyr::mutate(
        eic.df, real.time.min = real.time / 60,
        feature = ifelse(real.time.min >= start & real.time.min <= end,
          sample, "Non feature"),
        fill = ifelse(feature == "Non feature", feature, group),
        mz = round(mz, 4),
        anno.mz = paste("Precursor m/z:\n  ", mz - mz.tol, "~", mz + mz.tol),
        anno.rt = paste("RT (min):", round(rt, 1)),
        anno = paste0(anno.mz, "\n", anno.rt)
      )
      ## annotation (mz and rt)
      anno <- dplyr::select(eic.df, .features_id, int, real.time.min, contains("anno"))
      anno <- dplyr::group_by(anno, .features_id)
      anno <- dplyr::summarize(
        anno, anno.x = min(real.time.min, na.rm = T),
        anno.y = max(int, na.rm = T) * 3 / 4,
        anno = unique(anno)
      )
      data <- list(eic.df = eic.df, anno = anno)
    }
    if (!any(names(palette) == "Non feature")) {
      palette[[ "Non feature" ]] <- "grey95"
    }
    data$p <- ggplot(data[[ "eic.df" ]]) +
      geom_line(
        aes(x = real.time.min,
          y = int,
          group = sample,
          color = fill),
        lineend = "round") +
      labs(color = "Peak attribution", x = "RT (min)", y = "Intensity") +
      geom_text(data = data[[ "anno" ]],
        aes(x = anno.x, y = anno.y, label = anno),
        hjust = 0, fontface = "bold", family = .font) +
      scale_y_continuous(labels = scales::scientific) +
      facet_wrap( ~ paste("ID:", .features_id), scales = "free") +
      theme_minimal() +
      scale_color_manual(values = palette) +
      theme(text = element_text(family = .font),
        plot.background = element_rect(fill = "white", size = 0, color = "transparent"),
        strip.text = element_text(size = 12)) +
      geom_blank()
    return(data)
  }

#' @export bioc.par
#' @aliases bioc.par
#' @description \code{bioc.par}: ...
#' @rdname plot_EIC_stack
bioc.par <-
  function(cl = 4){
    BiocParallel::register(
      BiocParallel::bpstart(
        BiocParallel::MulticoreParam(cl)
      )
    )
  }
