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
    pg = "qiime",
    info = paste0("Tutorial: ",
      "https://docs.qiime.org/2023.7/tutorials/moving-pictures-usage/",
      "https://docs.qiime.org/2023.7/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq",
      "https://docs.qiime.org/2023.7/data-resources/",
      collapse = "\n"),
    cite = "[@ReproducibleIBolyen2019; @TheBiologicalMcdona2012; @Dada2HighResCallah2016; @ErrorCorrectinHamday2008; @MicrobialCommuHamday2009]",
    method = "`Qiime2` used for gut microbiome 16s rRNA analysis",
    tag = "16s:qiime2+mp",
    analysis = "Qiime2 16s-RNAseq 分析"
    ))

.qzv <- setClass("qzv", 
  contains = c("character"),
  representation = representation(pg = "character"),
  prototype = prototype(pg = "qiime"))

new_qzv <- function(..., lst = NULL, path, x, pg = NULL) {
  if (missing(x)) {
    x <- get("x", envir = parent.frame(1))
  }
  if (!is(x, "job")) {
    if (is.null(pg)) {
      stop("Do not know how to run Qiime2. At least `pg` or `x` (job) should passsed one.")
    }
  } else {
    pg <- pg(x)
  }
  if (is.null(lst))
    files <- unlist(list(...))
  else 
    files <- unlist(lst)
  if (missing(path)) {
    if (is.remote(x)) {
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
        list.files(file, "*.qzv", full.names = TRUE, recursive = TRUE)
      } else {
        file
      }
    })
  files <- unlist(files)
  names(files) <- get_realname(files)
  sapply(files, simplify = FALSE,
    function(file) {
      .qzv(file, pg = pg)
    })
}

setMethod("show", signature = c(object = "qzv"),
  function(object){
    cdRun(object@pg, " tools view ", as.character(object))
  })

job_qiime <- function(metadata, wd = "qiime_data", export = "qiime_export")
{
  x <- .job_qiime(object = metadata)
  x@params$wd <- wd
  x@params$export <- export
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
  function(x)
    # env_pattern = "qiime", env_path = "~/miniconda3/envs/", conda = "~/miniconda3/bin/conda"
  {
    step_message("Standby Qiime2 environment, and import data")
    # x@params$platform <- activate_qiime(env_pattern, env_path, conda)
    if (is.null(x@params$cdRun)) {
      x@params$cdRun <- cdRun
    }
    if (!is_qiime_file_exists("demux.qza")) {
      E(x@params$cdRun(
          pg(x), " tools import",
          " --type 'SampleData[PairedEndSequencesWithQuality]'",
          " --input-format PairedEndFastqManifestPhred33V2",
          " --input-path metadata.tsv",
          " --output-path demux.qza",
          path = x@params$wd)
      )
    }
    if (!is_qiime_file_exists("demux.qzv")) {
      E(x@params$cdRun(
          pg(x), " demux summarize ",
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
          pg(x), " dada2 denoise-paired ",
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
        pg(x), " metadata tabulate ",
        " --m-input-file denoising-stats.qza ",
        " --o-visualization denoising-stats.qzv",
        path = x@params$wd)
    )
    E(x@params$cdRun(
        pg(x), " feature-table summarize ",
        " --i-table table.qza ",
        " --m-sample-metadata-file metadata.tsv ",
        " --o-visualization table.qzv",
        path = x@params$wd)
    )
    E(x@params$cdRun(
        pg(x), " feature-table tabulate-seqs ",
        " --i-data rep-seqs.qza ",
        " --o-visualization rep-seqs.qzv",
        path = x@params$wd)
    )
    if (!is_qiime_file_exists("tree")) {
      E(x@params$cdRun(
          pg(x), " phylogeny align-to-tree-mafft-fasttree ",
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
          pg(x), " diversity core-metrics-phylogenetic ",
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
            pg(x), " diversity alpha-group-significance ",
            " --i-alpha-diversity diversity/", method, "_vector.qza ",
            " --m-metadata-file metadata.tsv ",
            " --o-visualization diversity_alpha_", group, "_significant/", method, ".qzv",
            path = x@params$wd
          )
        })
    )
    E(x@params$cdRun(
        pg(x), " diversity alpha-rarefaction ",
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
            pg(x), " diversity beta-group-significance ",
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
    qzvs <- sapply(basename(c(ga, gb)), simplify = FALSE,
      function(dir) {
        new_qzv(dir)
      })
    x@plots[[ 5 ]] <- qzvs
    return(x)
  })

setMethod("step6", signature = c(x = "job_qiime"),
  function(x, classifier = .prefix("weighted_silva_2023_07.qza", "db"))
  {
    step_message("Time consumed steps.")
    if (!file.exists(classifier)) {
      E(cdRun("wget -O ", classifier,
          " https://data.qiime2.org/2023.7/common/silva-138-99-nb-weighted-classifier.qza"))
    }
    classifier <- normalizePath(classifier)
    if (!is_qiime_file_exists("taxonomy.qza")) {
      E(x@params$cdRun(
          pg(x), " feature-classifier classify-sklearn ",
          " --i-classifier ", classifier,
          " --i-reads rep-seqs.qza ",
          " --o-classification taxonomy.qza",
          path = x@params$wd)
      )
    }
    E(x@params$cdRun(
        pg(x), " metadata tabulate ",
        " --m-input-file taxonomy.qza ",
        " --o-visualization taxonomy.qzv",
        path = x@params$wd
      )
    )
    E(x@params$cdRun(
        pg(x), " taxa barplot ",
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
  function(x, levels = 2:6, table = "table.qza", force = FALSE)
  {
    step_message("Differential analysis.")
    if (force || is.null(x@plots$step7$qzv_raw)) {
      files <- pbapply::pblapply(levels,
        function(level) {
          message()
          # The explanation about the test results:
          # https://forum.qiime2.org/t/how-to-interpret-ancom-results/1958
          # https://forum.qiime2.org/t/specify-w-cutoff-for-anacom/1844
          ancom_test(level, table, path = x@params$wd, x = x)
        })
      qzv_raw <- new_qzv(lst = files)
    } else {
      qzv_raw <- x@plots$step7$qzv_raw
    }
    lst <- lapply(qzv_raw,
      function(qzv) {
        dir <- write(qzv, output_dir = x$export, pg = suppressMessages(pg(x)))
        data <- ftibble(paste0(dir, "/data.tsv"))
        res <- ftibble(paste0(dir, "/ancom.tsv"))
        p <- plot_volcano.ancom(data, res)
        quant <- ftibble(paste0(dir, "/percent-abundances.tsv"))
        namel(p, quant)
      })
    t.quant <- lapply(lst,
      function(x) {
        data <- x$quant
        metadata <- data[1:2, ]
        data <- data[-(1:2), ]
        colnames(data) <- c("taxon", paste0(metadata[2, ], "_", metadata[1, ])[-1])
        data <- tidyr::gather(data, "group_percent", "value", -taxon)
        data <- tidyr::separate(data, group_percent, c("group", "percentile"), sep = "_")
        data <- dplyr::mutate(data, percentile = as.double(percentile),
          value = as.double(value))
        data
      })
    p.export <- lapply(lst, function(x) x$p)
    t.ancom <- lapply(p.export, function(x) as_tibble(x$data))
    t.ancom <- .set_lab(t.ancom, sig(x), names(t.ancom))
    lab(t.ancom) <- "Ancom test results"
    p.export <- lapply(p.export, wrap)
    p.export <- .set_lab(p.export, sig(x), names(p.export), "volcano")
    lab(p.export) <- "Ancom test visualization"
    ## plot Percentile
    n <- 0L
    p.quant <- lapply(t.quant,
      function(data) {
        n <<- n + 1L
        sigs <- dplyr::filter(t.ancom[[ n ]], significant)$id
        data <- dplyr::mutate(data,
          significant = ifelse(taxon %in% !!sigs, "Sig.", "Non."),
          taxon = stringr::str_trunc(taxon, 50, "left"),
          group = paste0("Group: ", group),
          percentile = paste0("Percentile: ", percentile, "%"),
          percentile = factor(percentile, levels = unique(percentile))
        )
        if (length(unique(data$taxon)) > 30) {
          data <- dplyr::filter(data, significant == "Sig.")
        }
        if (nrow(data)) {
          p <- ggplot(data) +
            geom_col(aes(x = taxon, y = value, fill = significant)) +
            coord_flip() +
            labs(x = "Abundance", y = "Taxonomy", fill = "Significance") +
            facet_grid(group ~ percentile) +
            scale_fill_manual(values = rev(color_set2())) +
            theme_minimal()
          wrap(p, 14, 10)
        }
      })
    p.quant <- .set_lab(p.quant, sig(x), names(p.quant), "Percentile abundance")
    t.quant <- .set_lab(t.quant, sig(x), names(t.quant), "data Percentile abundance")
    x@plots[[ 7 ]] <- namel(qzv_raw, p.export, p.quant)
    x@tables[[ 7 ]] <- namel(t.ancom, t.quant)
    return(x)
  })

setMethod("res", signature = c(x = "job_qiime"),
  function(x, level = c("6", "5", "4", "3", "2"))
  {
    if (x@step < 7L) {
      stop("x@step < 7L")
    }
    level <- match.arg(level)
    data <- x@tables$step7$t.ancom[[ paste0("ancom_test_group_level_", level) ]]
    data <- dplyr::filter(data, significant)
    data$id
  })

setMethod("pattern", signature = c(x = "job_qiime"),
  function(x, level = c("6", "5", "4", "3", "2"), not = c("GB", "^un$", "^bacterium$", "^bacteria$"))
  {
    if (x@step < 7L) {
      stop("x@step < 7L")
    }
    level <- match.arg(level)
    data <- x@tables$step7$t.ancom[[ paste0("ancom_test_group_level_", level) ]]
    data <- dplyr::filter(data, significant)
    if (!nrow(data)) {
      stop("No significant microbiota found.")
    }
    lst <- unlist(strsplit(data$id, ";"))
    lst <- gs(lst, "^.__([^_]+).*", "\\1")
    lst <- lst[ !grpl(lst, paste0(not, collapse = "|"), TRUE) ]
    lst[ !duplicated(lst) ]
  })

plot_volcano.ancom <- function(data, res) {
  data <- map(data, "id", res, "V1", "Reject null hypothesis", col = "significant")
  data <- dplyr::mutate(data, label = stringr::str_trunc(id, 30, side = "left"))
  p <- ggplot(data, aes(x = clr, y = W)) +
    geom_point(aes(size = W, color = significant), alpha = .5, stroke = 0) +
    ggrepel::geom_text_repel(aes(label = label)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    scale_color_manual(values = color_set()) +
    labs(x = "FALSE-statistic (clr)", y = "W-value")
  p
}

setMethod("write", signature = c(x = "qzv"),
  function(x, output = NULL, output_dir = "qiime_export", pg = get_fun("pg")("qiime", FALSE), overwrite = TRUE)
  {
    file <- as.character(x)
    if (is.null(output)) {
      output <- get_realname(file)
    }
    dir.create(output_dir, FALSE)
    cdRun(pg, " tools export ",
    " --input-path ", file,
    " --output-path ", name <- paste0(output_dir, "/", output))
    return(name)
  })

ancom_test <- function(level = 6, table = "table.qza", path = "./qiime_data", x) {
  if (!is.null(level)) {
    E(x@params$cdRun(
        pg(x), " taxa collapse ",
        " --i-table ", table,
        " --i-taxonomy taxonomy.qza ",
        " --p-level ", level,
        " --o-collapsed-table ", ntable <- paste0("table_level_", level, ".qza"),
        path = path
        ))
    table <- ntable
  }
  E(x@params$cdRun(
      pg(x), " composition add-pseudocount ",
      " --i-table ", table,
      " --o-composition-table ", com_table <- paste0("comp_table_level_", level, ".qza"),
      path = path
      ))
  E(x@params$cdRun(
      pg(x), " composition ancom ",
      " --i-table ", com_table,
      " --m-metadata-file metadata.tsv ",
      " --m-metadata-column ", "group",
      " --o-visualization ", res <- paste0("ancom_test_group_level_", level, ".qzv"),
      path = path
      ))
  return(res)
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
  res <- system(paste0("ssh ", remote, " '", expr_sys.file.exists(file), "'"), intern = TRUE)
  if (res == "TRUE") TRUE
  else FALSE
}

expr_sys.file.exists <- function(file) {
  paste0("if [ -e ", file, " ]; then echo TRUE; else echo FALSE; fi")
}

setMethod("set_remote", signature = c(x = "job_qiime"),
  function(x, path, wd = path,
    pattern = if (is.null(x@params$pattern_fq)) "fastq\\.gz$" else x@params$pattern_fq)
  {
    ## must be here
    files <- list.remote_rf(paste0(path, "/", object(x)$Run), pattern)
    metadata <- try_fqs_meta(object(x), files, filter = TRUE)
    print(metadata)
    object(x) <- metadata
    meta_file <- tempfile("metadata", fileext = "tsv")
    write_tsv(metadata, meta_file)
    system(paste0("scp ", meta_file, " ", remote, ":", wd, "/metadata.tsv"))
    x@params$meta_file <- paste0(wd, "/metadata.tsv")
    x@params$cdRun <- remoteRun
    x@params$wd <- wd
    return(x)
  })

try_fqs_meta <- function(metadata, filepath, filter = FALSE) {
  if (any(duplicated(metadata[[ "dirs" ]]))) {
    filepath <- mapply(get_realname(metadata$reports), filepath, SIMPLIFY = FALSE,
      FUN = function(name, files) {
        files <- files[ grpl(files, name, fixed = TRUE) ]
        if (length(files) > 2) {
          files <- files[ grpl(files, paste0("/", name, "[^/]+$")) ]
        }
        files
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

activate_qiime <- function(env_pattern = "qiime", env_path = pg("conda_env"), conda = pg("conda")) {
  activate_base(env_pattern, env_path, conda)
}

get_taxon_data <- function(file = .prefix("qiime2/taxonomy.tsv", "db")) {
  if (!file.exists(file)) {
    # system("wget https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza")
    stop("Please download the taxonomy file from qiime2 website.")
  }
  ftibble(file)[, -1]
}

query_class <- function(x, level = c("g", "p", "c", "o", "f", "s"),
  use.first = TRUE, fun_data = get_taxon_data, strict = TRUE)
{
  if (is.null(ref <- getOption("taxon_data", NULL))) {
    message("Got and set global taxonomy data for querying")
    ref <- fun_data()
    options(taxon_data = ref)
  }
  if (use.first) {
    x <- strx(x, "[a-zA-Z0-9\\-]{2,}")
  }
  level <- match.arg(level)
  ref <- dplyr::mutate(ref, Taxon = gs(Taxon, paste0("(.*", level, "__", "[^_]*?;).*"), "\\1"))
  ref <- dplyr::distinct(ref)
  xUnique <- unique(x)
  res <- .find_and_sort_strings(ref$Taxon,
    if (strict) paste0("__", xUnique, ";") else xUnique, TRUE)
  res <- vapply(res,
    function(x) {
      if (identical(x, character(0))) {
        NA_character_ 
      } else {
        x <- unique(x)
        if (length(x) > 1) {
          warning("length(x) > 1, use first: ", x[1])
          x <- x[1]
        }
        x
      }
    },
    FUN.VALUE = character(1), USE.NAMES = FALSE
  )
  res[ match(x, xUnique) ]
}


