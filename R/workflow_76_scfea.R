# ==========================================================================
# workflow of scfea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scfea <- setClass("job_scfea", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "scfea",
    info = c("https://github.com/changwn/scFEA/tree/master"),
    cite = "[@AGraphNeuralAlgham2021]",
    method = "The `scFEA` (python) was used to estimate cell-wise metabolic via single cell RNA-seq data",
    tag = "scrna:flux"
    ))

setGeneric("asjob_scfea", 
  function(x, ...) standardGeneric("asjob_scfea"))

setMethod("asjob_scfea", signature = c(x = "job_seurat"),
  function(x, cells = NULL, org = c("mouse", "human"), assay = "RNA", dir = "scfea", ...)
  {
    org <- match.arg(org)
    message("Use org: ", org)
    if (org == "mouse") {
      moduleGene_file <- "module_gene_complete_mouse_m168.csv"
    } else if (org == "human") {
      moduleGene_file <- "module_gene_m168.csv"
    }
    if (T) {
      file <- paste0(pg("scfea_db"), "/", moduleGene_file)
      fun_test <- function(file) {
        x <- readLines(file)[-1]
        x <- lapply(strsplit(x, ","), function(x) x[-1])
        unique(unlist(x))
      }
      refs <- fun_test(file)
    }
    data <- object(x)@assays[[ assay ]]@counts
    data <- data[ rownames(data) %in% refs, ]
    if (!is.null(cells)) {
      data <- data[, cells]
    }
    message("Dim: ", paste0(dim(data), collapse = ", "))
    message("Max: ", max(apply(data, 2, max)))
    if (dir.exists(dir)) {
      if (usethis::ui_yeah("Directory of `dir` exists, remove that?")) {
        unlink(dir, recursive = T)
      }
    }
    dir.create(dir, F)
    expr_file <- paste0(dir, "/data.csv")
    data.table::fwrite(data.frame(data, check.names = F), expr_file, row.names = T)
    job_scfea(expr_file, org = org)
  })

job_scfea <- function(expr_file, org = c("mouse", "human"), test = F)
{
  org <- match.arg(org)
  x <- .job_scfea()
  x$dir <- get_path(expr_file)
  x$input_file <- get_filename(expr_file)
  if (org == "mouse") {
    x$moduleGene_file <- "module_gene_complete_mouse_m168.csv"
    x$stoichiometry_matrix <- "cmMat_complete_mouse_c70_m168.csv"
    x$cName_file <- "cName_complete_mouse_c70_m168.csv"
    x$module_annotation <- "Human_M168_information.symbols.csv"
  } else if (org == "human") {
    x$moduleGene_file <- "module_gene_m168.csv"
    x$stoichiometry_matrix <- "cmMat_c70_m168.csv"
    x$cName_file <- "cName_c70_m168.csv"
    x$module_annotation <- "Human_M168_information.symbols.csv"
  }
  if (test) {
    file <- paste0(pg("scfea_db"), "/", x$moduleGene_file)
    fun_test <- function(file) {
      x <- readLines(file)[-1]
      x <- lapply(strsplit(x, ","), function(x) x[-1])
      unique(unlist(x))
    }
    refs <- fun_test(file)
    data <- ftibble(expr_file)
    print(table(refs %in% data[[1]]))
  }
  return(x)
}

setMethod("step0", signature = c(x = "job_scfea"),
  function(x)
  {
    step_message("Prepare your data with function `asjob_scfea`.")
  })

setMethod("step1", signature = c(x = "job_scfea"),
  function(x)
  {
    step_message("Run prediction.")
    if (!is.remote(x)) {
      dir.create("output")
    }
    rem_run(
      pg("scfeaPython", is.remote(x)), " ", pg("scfea", is.remote(x)),
      " --data_dir ", pg("scfea_db", is.remote(x)),
      " --input_dir ", x$dir,
      " --test_file ", x$input_file,
      " --res_dir ", x$dir,
      " --moduleGene_file ", x$moduleGene_file,
      " --stoichiometry_matrix ", x$stoichiometry_matrix,
      " --cName_file ", x$cName_file,
      " --output_flux_file ", "flux.csv",
      " --output_balance_file ", "balance.csv"
    )
    if (is.remote(x)) {
      cdRun("scp -r ", x$remote, ":", x$wd, "/* ", x$map_local)
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_scfea"),
  function(x){
    step_message("Collate results.")
    if (is.remote(x)) {
      dir <- x$map_local
    } else {
      dir <- x$dir
    }
    t.balance <- ftibble(paste0(dir, "/", "balance.csv"))
    t.flux <- ftibble(paste0(dir, "/", "flux.csv"))
    fun <- function(file) {
      lst <- strsplit(readLines(file)[-1], ",")
      names(lst) <- lapply(lst, function(x) x[1])
      lst <- lapply(lst, function(x) x[-1])
      as_df.lst(lst, "module", "gene")
    }
    dir_db <- pg("scfea_db")
    t.moduleGenes <- fun(paste0(dir_db, "/", x$moduleGene_file))
    maybePlot <- list.files(dir, "loss_[0-9]+-[0-9]+\\.png", full.names = T)
    if (length(maybePlot)) {
      p.loss <- .file_fig(maybePlot)
      p.loss <- .set_lab(p.loss, sig(x), "Convergency of the loss terms during training")
    } else {
      p.loss <- NULL
    }
    t.anno <- ftibble(paste0(dir_db, "/", x$module_annotation))
    t.anno <- dplyr::mutate(t.anno, name = paste0(Compound_IN_name, " -> ", Compound_OUT_name))
    x@tables[[ 2 ]] <- namel(t.flux, t.anno, t.balance, t.moduleGenes)
    x@plots[[ 2 ]] <- namel(p.loss)
    return(x)
  })

setMethod("cal_corp", signature = c(x = "job_limma", y = "job_seurat"),
  function(x, y, from = tops$name, to = tops$gene, names = NULL,
    tops = y@tables$step2$tops[[ 1 ]])
  {
    if (is.null(x$from_scfea)) {
      stop("The 'job_limma' is not converted from job_scfea.")
    }
    if (is(to, "list")) {
      groupTo <- to
      to <- unique(unlist(to))
    }
    stop("...")
  })

setMethod("asjob_limma", signature = c(x = "job_scfea"),
  function(x, metadata, group, scale_sample = T, scale_var = F, ...)
  {
    if (missing(metadata) || missing(group)) {
      stop("missing(metadata) || missing(group)")
    }
    counts <- data.frame(x@tables$step2$t.flux)
    rownames(counts) <- counts[[1]]
    counts <- counts[, -1]
    if (scale_sample) {
      counts <- scale(counts, ...)
    }
    counts <- t(counts)
    if (scale_var) {
      counts <- scale(counts, ...)
    }
    counts <- as_tibble(counts, .name_repair = "minimal")
    genes <- x@tables$step2$t.anno
    annoExtra <- reframe_col(as_tibble(x@tables$step2$t.moduleGenes), "gene", list)
    genes <- map(genes, colnames(genes)[1], annoExtra, "module", "gene", col = "gene")
    if (!any(metadata[[1]] %in% colnames(counts))) {
      stop("Use first column of `metadata` as sample name (cell name), but not match any.")
    } else {
      message("Use first column of `metadata` as sample name (cell name).")
      metadata <- dplyr::filter(metadata, !!rlang::sym(colnames(metadata)[1]) %in% !!colnames(counts))
      metadata <- dplyr::rename(metadata, sample = !!rlang::sym(colnames(metadata)[1]),
        group = !!rlang::sym(group))
      metadata <- dplyr::relocate(metadata, sample, group)
    }
    cli::cli_alert_info("job_limma_normed")
    x <- job_limma_normed(counts, metadata, genes = genes)
    x$from_scfea <- T
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_scfea"),
  function(x, wd = "/data/hlc/scfea", remote = "remote")
  {
    if (!is.remote(x)) {
      x$set_remote <- T
      x$remote <- remote
      if (rem_file.exists(wd)) {
        isThat <- usethis::ui_yeah("Dir exists, remove all files ?")
        if (isThat) {
          cdRun("ssh ", remote, " 'rm -r ", wd, "'")
        }
      }
      rem_dir.create(wd)
      x$wd <- wd
      rem_dir.create("output")
      cdRun("scp ", paste0(x$dir, "/data.csv"), " ", remote, ":", wd)
      x$map_local <- x$dir
      x$dir <- "."
    }
    return(x)
  })
