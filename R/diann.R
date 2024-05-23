# ==========================================================================
# Diann
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

diann <- function(targets, fasta,
  savepath = paste0(Sys.Date(), "_Diann"),
  out = file.path(savepath, "report.tsv"), lib = NULL,
  out.lib = file.path(savepath, "report-lib.tsv"),
  min.fr.mz = 200, max.fr.mz = 1800, cut = "K*,R*",
  missed.cleavages = 1, min.pep.len = 7, max.pep.len = 30,
  min.pr.mz = 300, max.pr.mz = 1800, min.pr.charge = 1, max.pr.charge = 4,
  qvalue = .01, threads = 10, verbose = 1,
  command = "/usr/diann/1.8.1/diann-1.8.1",
  pattern = "\\.mzML|\\.wiff", recursive = F,
  notRun = NULL, more.args = NULL)
{
  if (length(targets == 1) & dir.exists(targets[1])) {
    if (dir.exists(targets)) {
      files <- list.files(targets, pattern, full.names = T, recursive = recursive)
      if (!length(files)) {
        stop("No files found of pattern: ", pattern)
      }
    }
  } else {
    if (!all(file.exists(targets))) {
      stop("Check the input files, some not exists.")
    }
    files <- targets
  }
  if (!dir.exists(savepath)) {
    dir.create(savepath)
  }
  normp <- function(x) normalizePath(x, mustWork = T)
  savepath <- normp(savepath)
  if (any(file.exists(c(out, out.lib)))) {
    res <- menu(c("Overwrite the files.", "Stop running."))
    if (res == 2) {
      stop("Stop running.")
    }
  }
  ifrun <- function(x) {
    if (is.null(notRun)) {
      return(x)
    }
    check <- vapply(notRun, function(i) grepl(x, i), logical(1))
    if (any(check)) NULL else x
  }
  expr <- paste0(command, " ",
    paste0(paste0(" --f ", normp(files)), collapse = " "),
    if (is.null(lib)) NULL else paste0("  --lib ", normp(lib)),
    "  --threads ", threads,
    "  --verbose ", verbose,
    "  --out ", out,
    "  --qvalue ", qvalue,
    ifrun("  --matrices "),
    "  --out-lib ", out.lib,
    ifrun("  --gen-spec-lib "),
    ifrun("  --predictor "),
    ifrun("  --reannotate "),
    "  --fasta ", normp(fasta),
    ifrun("  --fasta-search "),
    "  --min-fr-mz ", min.fr.mz,
    "  --max-fr-mz ", max.fr.mz,
    ifrun("  --met-excision "),
    "  --cut ", cut,
    "  --missed-cleavages ", missed.cleavages,
    "  --min-pep-len ", min.pep.len,
    "  --max-pep-len ", max.pep.len,
    "  --min-pr-mz ", min.pr.mz,
    "  --max-pr-mz ", max.pr.mz,
    "  --min-pr-charge ", min.pr.charge,
    "  --max-pr-charge ", max.pr.charge,
    ifrun("  --unimod4 "),
    ifrun("  --reanalyse "),
    ifrun("  --relaxed-prot-inf "),
    ifrun("  --smart-profiling "),
    ifrun("  --peak-center "),
    ifrun("  --no-ifs-removal ")
  )
  if (!is.null(more.args)) {
    expr <- paste0(expr, " ", paste0(more.args, collapse = " "))
  }
  res <- try(system(expr), T)
  runlog <- tempfile(pattern = "diann_script", fileext = ".sh")
  writeLines(expr, runlog)
  message("\nThe script for Diann running: ", runlog)
}
