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
      collapse = "\n"),
    cite = "[@ReproducibleIBolyen2019; @TheBiologicalMcdona2012; @Dada2HighResCallah2016; @ErrorCorrectinHamday2008; @MicrobialCommuHamday2009]",
    method = "`Qiime2` used for gut microbiome 16s rRNA analysis"
    ))

.qzv <- setClass("qzv", 
  contains = c("character"),
  representation = representation(),
  prototype = NULL)

new_qzv <- function(..., lst = NULL, path) {
  if (is.null(lst))
    files <- unlist(list(...))
  else 
    files <- unlist(lst)
  if (missing(path)) {
    x <- get("x", envir = parent.frame(1))
    if (!is.null(x@params$set_remote)) {
      path <- x@params$map_local
      pbapply::pblapply(files,
        function(file) {
          system(paste0("scp -r ", x@params$remote, ":", x@params$wd, "/",
              file, " ", path))
        })
    } else {
      path <- x@params$wd
    }
  }
  files <- paste0(path, "/", files)
  files <- lapply(files,
    function(file) {
      if (dir.exists(file)) {
        list.files(file, "*.qzv", full.names = T, recursive = T)
      } else {
        file
      }
    })
  files <- unlist(files)
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

job_qiime <- function(metadata, wd = "qiime_data")
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
    x@params$platform <- activate_qiime(env_pattern, env_path, conda)
    if (is.null(x@params$cdRun)) {
      x@params$cdRun <- cdRun
    }
    if (!is_qiime_file_exists("demux.qza")) {
      E(x@params$cdRun(
          "qiime tools import",
          " --type 'SampleData[PairedEndSequencesWithQuality]'",
          " --input-format PairedEndFastqManifestPhred33V2",
          " --input-path metadata.tsv",
          " --output-path demux.qza",
          path = x@params$wd)
      )
    }
    if (!is_qiime_file_exists("demux.qzv")) {
      E(x@params$cdRun(
          "qiime demux summarize ",
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
      E(x@params$cdRun(
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
    x@params$table_qza <- paste0(x@params$wd, "/table.qza")
    return(x)
  })

setMethod("step3", signature = c(x = "job_qiime"),
  function(x){
    step_message("Some visualization and 'align-to-tree-mafft-fasttree'")
    E(x@params$cdRun(
        "qiime metadata tabulate ",
        " --m-input-file denoising-stats.qza ",
        " --o-visualization denoising-stats.qzv",
        path = x@params$wd)
    )
    E(x@params$cdRun(
        "qiime feature-table summarize ",
        " --i-table table.qza ",
        " --m-sample-metadata-file metadata.tsv ",
        " --o-visualization table.qzv",
        path = x@params$wd)
    )
    E(x@params$cdRun(
        "qiime feature-table tabulate-seqs ",
        " --i-data rep-seqs.qza ",
        " --o-visualization rep-seqs.qzv",
        path = x@params$wd)
    )
    if (!is_qiime_file_exists("tree")) {
      E(x@params$cdRun(
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
      E(x@params$cdRun(
          "qiime diversity core-metrics-phylogenetic ",
          " --i-phylogeny tree/rooted_tree.qza ",
          " --i-table table.qza ",
          " --p-sampling-depth ", x@params$min,
          " --m-metadata-file metadata.tsv ",
          " --output-dir diversity",
          path = x@params$wd)
      )
    }
    x@plots[[ 4 ]] <- new_qzv("/diversity")
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
          x@params$cdRun(
            "qiime diversity alpha-group-significance ",
            " --i-alpha-diversity diversity/", method, "_vector.qza ",
            " --m-metadata-file metadata.tsv ",
            " --o-visualization diversity_alpha_", group, "_significant/", method, ".qzv",
            path = x@params$wd
          )
        })
    )
    E(x@params$cdRun(
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
          x@params$cdRun(
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
        new_qzv(dir)
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
    classifier <- normalizePath(classifier)
    if (!is_qiime_file_exists("taxonomy.qza")) {
      E(x@params$cdRun(
          "qiime feature-classifier classify-sklearn ",
          " --i-classifier ", classifier,
          " --i-reads rep-seqs.qza ",
          " --o-classification taxonomy.qza",
          path = x@params$wd)
      )
    }
    E(x@params$cdRun(
        "qiime metadata tabulate ",
        " --m-input-file taxonomy.qza ",
        " --o-visualization taxonomy.qzv",
        path = x@params$wd
      )
    )
    E(x@params$cdRun(
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
        ancom_test(level, table, path = x@params$wd, x = x)
      })
    x@plots[[ 7 ]] <- new_qzv(lst = files)
    return(x)
  })

ancom_test <- function(level = 6, table = "table.qza", path = "./sra_data", x) {
  if (!is.null(level)) {
    E(x@params$cdRun(
        "qiime taxa collapse ",
        " --i-table ", table,
        " --i-taxonomy taxonomy.qza ",
        " --p-level ", level,
        " --o-collapsed-table ", ntable <- paste0("table_level_", level, ".qza"),
        path = path
        ))
    table <- ntable
  }
  E(x@params$cdRun(
      "qiime composition add-pseudocount ",
      " --i-table ", table,
      " --o-composition-table ", com_table <- paste0("comp_table_level_", level, ".qza"),
      path = path
      ))
  E(x@params$cdRun(
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
    if (!is.null(x@params$set_remote)) {
      message(crayon::yellow("Check remote ..."))
      res <- remote_file.exists(paste0(path, "/", file), remote = x@params$remote)
      return(res)
    }
  }
  file.exists(paste0(path, "/", file))
}

remote_file.exists <- function(file, remote = "remote") {
  if (missing(remote)) {
    x <- get("x", envir = parent.frame(1))
    remote <- x@params$remote
  }
  res <- system(paste0("ssh ", remote, " '", expr_sys.file.exists(file), "'"), intern = T)
  if (res == "T") T
  else F
}

expr_sys.file.exists <- function(file) {
  paste0("if [ -e ", file, " ]; then echo T; else echo F; fi")
}

setMethod("set_remote", signature = c(x = "job_qiime"),
  function(x, path, wd = path,
    pattern = if (is.null(x@params$pattern_fq)) "fastq\\.gz$" else x@params$pattern_fq,
    tmpdir = "/data/hlc/tmp", map_local = "qiime_local", remote = "remote")
  {
    ## must be here
    x@params$remote <- remote
    files <- list.remote_rf(paste0(path, "/", object(x)$Run), pattern)
    metadata <- try_fqs_meta(object(x), files, filter = T)
    print(metadata)
    object(x) <- metadata
    meta_file <- tempfile("metadata", fileext = "tsv")
    write_tsv(metadata, meta_file)
    system(paste0("scp ", meta_file, " ", remote, ":", wd, "/metadata.tsv"))
    x@params$meta_file <- paste0(wd, "/metadata.tsv")
    x@params$cdRun <- remoteRun
    x@params$set_remote <- T
    x@params$map_local <- map_local
    x@params$postfix <- function(x) {
      x[1] <- gs(x[1], "^qiime", "~/miniconda3/bin/conda run -n qiime2 qiime")
      x
    }
    x@params$wd <- wd
    x@params$tmpdir <- tmpdir
    return(x)
  })

try_fqs_meta <- function(metadata, filepath, filter = F) {
  if (any(duplicated(metadata[[ "dirs" ]]))) {
    filepath <- mapply(get_realname(metadata$reports), filepath, SIMPLIFY = F,
      FUN = function(name, files) {
        files[ grpl(files, name, fixed = T) ]
      })
  }
  fun <- function(lst, n) {
    vapply(lst,
      function(vec) {
        vec <- sort(vec)
        if (length(vec) >= n)
          vec[n]
        else
          character(1)
      }, character(1))
  }
  metadata[[ "forward-absolute-filepath" ]] <- fun(filepath, 1)
  metadata[[ "reverse-absolute-filepath" ]] <- fun(filepath, 2)
  if (filter) {
    metadata <- filter(metadata, `forward-absolute-filepath` != "")
  }
  metadata
}

activate_qiime <- function(env_pattern = "qiime", env_path = "~/miniconda3/envs/", conda = "~/miniconda3/bin/conda") {
  activate_base(env_pattern, env_path, conda)
}
