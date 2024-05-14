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
    method = "The `scFEA` (python) was used to estimate cell-wise metabolic via single cell RNA-seq data"
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
    data.table::fwrite(data.frame(data), expr_file, row.names = T)
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
  } else if (org == "human") {
    x$moduleGene_file <- "module_gene_m168.csv"
    x$stoichiometry_matrix <- "cmMat_c70_m168.csv"
    x$cName_file <- "cName_c70_m168.csv"
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
