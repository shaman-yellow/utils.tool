# ==========================================================================
# workflow of qiime
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_qiime <- setClass("job_qiime", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("Tutorial: ",
      "https://docs.qiime.org/2023.7/tutorials/moving-pictures-usage/",
      "https://docs.qiime.org/2023.7/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq",
      "https://docs.qiime.org/2023.7/data-resources/",
      collapse = "\n"
    )
    ))

.qzv <- setClass("qzv", 
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

new_qzv <- function(..., lst = NULL, path) {
  if (missing(path)) {
    x <- get("x", envir = parent.frame(1))
    path <- x@params$wd
  }
  if (is.null(lst))
    files <- unlist(list(...))
  else 
    files <- unlist(lst)
  if (!is.null(path))
    files <- paste0(path, "/", files)
  names(files) <- get_realname(files)
  sapply(files, simplify = F,
    function(file) {
      .qzv(file)
    })
}

setMethod("show", signature = c(object = "qzv"),
  function(object){
    system(paste0("qiime tools view ", as.character(object)))
  })

job_qiime <- function(metadata, wd)
{
  x <- .job_qiime(object = metadata)
  x@params$wd <- wd
  write_tsv(metadata, meta_file <- paste0(x@params$wd, "/metadata.tsv"))
  x@params$meta_file <- meta_file
  return(x)
}

setMethod("step0", signature = c(x = "job_qiime"),
  function(x){
    step_message("Prepare your data with function `job_qiime`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_qiime"),
  function(x, env_pattern = "qiime", env_path = "~/miniconda3/envs/", conda = "~/miniconda3/bin/conda")
  {
    step_message("Standby Qiime2 environment, and import data")
    conda_env <- e(filter(reticulate::conda_list(), grepl(env_pattern, name))$name)
    e(base::Sys.setenv(RETICULATE_PYTHON = paste0(env_path, "/", conda_env, "/bin/python")))
    e(reticulate::py_config())
    e(reticulate::use_condaenv(conda_env, conda, required = TRUE))
    platform <- e(reticulate::import("platform"))
    x@params$platform <- platform
    if (!is_qiime_file_exists("demux.qza")) {
      E(cdRun(
          "qiime tools import",
          " --type 'SampleData[PairedEndSequencesWithQuality]'",
          " --input-format PairedEndFastqManifestPhred33V2",
          " --input-path metadata.tsv",
          " --output-path demux.qza",
          path = x@params$wd)
      )
    }
    if (!is_qiime_file_exists("demux.qzv")) {
      E(cdRun("qiime demux summarize ",
          " --i-data demux.qza ",
          " --o-visualization demux.qzv",
          path = x@params$wd)
      )
    }
    x@plots[[ 1 ]] <- new_qzv("demux.qzv")
    return(x)
  })

setMethod("step2", signature = c(x = "job_qiime"),
  function(x, len_f, len_r, left_f = 0, left_r = 0, workers = 7){
    step_message("Time consumed running.")
    if (!is_qiime_file_exists("table.qza")) {
      E(cdRun(
          "qiime dada2 denoise-paired ",
          " --i-demultiplexed-seqs demux.qza ",
          " --p-n-threads ", workers,
          " --p-trim-left-f ", left_f, " --p-trim-left-r ", left_r,
          " --p-trunc-len-f ", len_f, " --p-trunc-len-r ", len_r,
          " --o-table table.qza ",
          " --o-representative-sequences rep-seqs.qza ",
          " --o-denoising-stats denoising-stats.qza",
          path = x@params$wd)
      )
    }
    x@params$table_qza <- paste0(path, "/table.qza")
    return(x)
  })

setMethod("step3", signature = c(x = "job_qiime"),
  function(x){
    step_message("Some visualization and 'align-to-tree-mafft-fasttree'")
    E(cdRun(
        "qiime metadata tabulate ",
        " --m-input-file denoising-stats.qza ",
        " --o-visualization denoising-stats.qzv",
        path = x@params$wd)
    )
    E(cdRun(
        "qiime feature-table summarize ",
        " --i-table table.qza ",
        " --m-sample-metadata-file metadata.tsv ",
        " --o-visualization table.qzv",
        path = x@params$wd)
    )
    E(cdRun(
        "qiime feature-table tabulate-seqs ",
        " --i-data rep-seqs.qza ",
        " --o-visualization rep-seqs.qzv",
        path = x@params$wd)
    )
    if (!is_qiime_file_exists("tree")) {
      E(cdRun(
          "qiime phylogeny align-to-tree-mafft-fasttree ",
          " --i-sequences rep-seqs.qza ",
          " --output-dir tree",
          path = x@params$wd)
      )
    }
    x@plots[[ 3 ]] <- new_qzv("denoising-stats.qzv", "table.qzv", "rep-seqs.qzv")
    return(x)
  })

setMethod("step4", signature = c(x = "job_qiime"),
  function(x, min){
    step_message("Before diversity of alpha and beta.")
    x@params$min <- min
    if (!is_qiime_file_exists("diversity")) {
      E(cdRun(
          "qiime diversity core-metrics-phylogenetic ",
          " --i-phylogeny tree/rooted_tree.qza ",
          " --i-table table.qza ",
          " --p-sampling-depth ", x@params$min,
          " --m-metadata-file metadata.tsv ",
          " --output-dir diversity",
          path = x@params$wd)
      )
    }
    files <- list.files(paste0(x@params$wd, "/diversity"), "\\.qzv$", full.names = T)
    x@plots[[ 4 ]] <- new_qzv(lst = files, path = NULL)
    return(x)
  })

setMethod("step5", signature = c(x = "job_qiime"),
  function(x, max, group = "group"){
    step_message("Alpha and Beta diversity.")
    x@params$max <- max
    dir.create(ga <- paste0(x@params$wd, "/diversity_alpha_", group, "_significant"))
    dir.create(gb <- paste0(x@params$wd, "/diversity_beta_", group, "_significant"))
    ## alpha
    E(pbapply::pblapply(c("faith_pd", "shannon", "observed_features", "evenness"),
        function(method) {
          cdRun(
            "qiime diversity alpha-group-significance ",
            " --i-alpha-diversity diversity/", method, "_vector.qza ",
            " --m-metadata-file metadata.tsv ",
            " --o-visualization diversity_alpha_", group, "_significant/", method, ".qzv",
            path = x@params$wd
          )
        })
    )
    E(cdRun(
        "qiime diversity alpha-rarefaction ",
        " --i-table table.qza ",
        " --i-phylogeny tree/rooted_tree.qza ",
        " --p-max-depth ", max,
        " --m-metadata-file metadata.tsv ",
        " --o-visualization alpha-rarefaction.qzv",
        path = x@params$wd)
    )
    ## beta
    E(pbapply::pblapply(c("unweighted_unifrac", "bray_curtis", "weighted_unifrac", "jaccard"),
        function(index) {
          cdRun(
            "qiime diversity beta-group-significance ",
            " --i-distance-matrix diversity/", index, "_distance_matrix.qza ",
            " --m-metadata-file metadata.tsv ",
            " --m-metadata-column ", group,
            " --p-pairwise ",
            " --o-visualization ", "diversity_beta_", group, "_significant/",
            index, "_", group, "_significance.qzv",
            path = x@params$wd
          )
        })
    )
    qzvs <- sapply(get_filename(c(ga, gb)), simplify = F,
      function(dir) {
        files <- list.files(paste0(x@params$wd, "/", dir), "\\.qzv$", full.names = T)
        new_qzv(lst = files, path = NULL)
      })
    x@plots[[ 5 ]] <- qzvs
    return(x)
  })

setMethod("step6", signature = c(x = "job_qiime"),
  function(x, classifier = "../weighted_silva_2023_07.qza")
  {
    step_message("Time consumed steps.")
    if (!file.exists(classifier)) {
      E(cdRun("wget -O ", classifier,
          " https://data.qiime2.org/2023.7/common/silva-138-99-nb-weighted-classifier.qza"))
    }
    if (!is_qiime_file_exists("taxonomy.qza")) {
      E(cdRun(
          "qiime feature-classifier classify-sklearn ",
          " --i-classifier ", normalizePath(classifier),
          " --i-reads rep-seqs.qza ",
          " --o-classification taxonomy.qza",
          path = x@params$wd)
      )
    }
    E(cdRun(
        "qiime metadata tabulate ",
        " --m-input-file taxonomy.qza ",
        " --o-visualization taxonomy.qzv",
        path = x@params$wd
      )
    )
    E(cdRun(
        "qiime taxa barplot ",
        " --i-table table.qza ",
        " --i-taxonomy taxonomy.qza ",
        " --m-metadata-file metadata.tsv ",
        " --o-visualization taxa-bar-plots.qzv",
        path = x@params$wd
      )
    )
    x@plots[[ 6 ]] <- new_qzv("taxa-bar-plots.qzv", "taxonomy.qzv")
    x@params$taxonomy <- paste0(x@params$wd, "/taxonomy.qza")
    return(x)
  })

setMethod("step7", signature = c(x = "job_qiime"),
  function(x, levels = 4:6, table = "table.qza"){
    step_message("Differential analysis.")
    files <- pbapply::pblapply(levels,
      function(level) {
        message()
        ancom_test(level, table, path = x@params$wd)
      })
    x@plots[[ 7 ]] <- new_qzv(lst = files, path = NULL)
    return(x)
  })

ancom_test <- function(level = 6, table = "table.qza", path = "./sra_data") {
  if (!is.null(level)) {
    E(cdRun(
        "qiime taxa collapse ",
        " --i-table ", table,
        " --i-taxonomy taxonomy.qza ",
        " --p-level ", level,
        " --o-collapsed-table ", ntable <- paste0("table_level_", level, ".qza"),
        path = "./sra_data"
        ))
    table <- ntable
  }
  E(cdRun(
      "qiime composition add-pseudocount ",
      " --i-table ", table,
      " --o-composition-table ", com_table <- paste0("comp_table_level_", level, ".qza"),
      path = path
      ))
  E(cdRun(
      "qiime composition ancom ",
      " --i-table ", com_table,
      " --m-metadata-file metadata.tsv ",
      " --m-metadata-column ", "group",
      " --o-visualization ", res <- paste0("ancom_test_group_level_", level, ".qzv"),
      path = path
      ))
  paste0(path, "/", res)
}

qiime_vis <- function(file) {
  E(cdRun("qiime tools view ", file))
}

is_qiime_file_exists <- function(file, path) {
  if (missing(path)) {
    x <- get("x", envir = parent.frame(1))
    path <- x@params$wd
  }
  file.exists(paste0(path, "/", file))
}
