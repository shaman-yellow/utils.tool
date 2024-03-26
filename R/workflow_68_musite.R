# ==========================================================================
# workflow of musite
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_musite <- setClass("job_musite", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "musitePython",
    info = c("https://github.com/duolinwang/MusiteDeep_web"),
    cite = "[@MusitedeepADWang2020]",
    method = "Python tool `MusiteDeep` was used for protein post-translational modification site prediction and visualization"
    ))

job_musite <- function(symbols)
{
  .job_musite(object = symbols)
}

setMethod("step0", signature = c(x = "job_musite"),
  function(x){
    step_message("Prepare your data with function `job_musite`.")
  })

setMethod("step1", signature = c(x = "job_musite"),
  function(x, org = c("hsa", "mmu"), save = "Seq")
  {
    step_message("Prepare protein sequences.")
    org <- match.arg(org)
    x$mart <- new_biomart(org)
    x$seqs <- get_seq.pro(object(x), x$mart)
    x$seqs_file <- write(x$seqs$fasta, save, max = NULL)
    return(x)
  })

setMethod("step2", signature = c(x = "job_musite"),
  function(x,
    type = c("Hydroxylysine", "Hydroxyproline", "Methylarginine",
      "Methyllysine", "N-linked_glycosylation", "N6-acetyllysine",
      "O-linked_glycosylation", "Phosphoserine_Phosphothreonine",
      "Phosphotyrosine", "Pyrrolidone_carboxylic_acid", "S-palmitoyl_cysteine",
      "SUMOylation", "Ubiquitination"
      ),
    cutoff = .5,
    output = timeName("musite"))
  {
    step_message("Compute PTM position.")
    if (is.null(x$output)) {
      x$output <- output
      dir.create(x$output)
    }
    x$type <- match.arg(type)
    message("Use PTM type: ", x$type)
    models <- paste0(paste0(pg("musiteModel"), "/", x$type), collapse = ";")
    message("Use models: ", models)
    output_file <- paste0(x$output, "/res_results.txt")
    if (!file.exists(output_file)) {
      cdRun(pg(x), " ", pg("musitePTM"),
        " -input ", x$seqs_file,
        " -output ", x$output, "/res",
        " -model-prefix \"", models, "\""
      )
    }
    lst <- sep_list(readLines(output_file), "^>", T)
    t.data <- lapply(lst[-1],
      function(x) {
        data.table::fread(text = x[-1])
      })
    names(t.data) <- gs(vapply(lst[-1], function(x) x[1], character(1)), "^>", "")
    t.data <- frbind(t.data, idcol = "Sequence_name")
    colnames(t.data)[-1] <- c("Position",	"Residue", "PTMscores", "Cutoff")
    t.data <- dplyr::select(t.data, -Cutoff)
    t.data <- tidyr::separate(t.data, PTMscores, c("PTM_type", "PTM_score"), ":")
    t.data <- dplyr::mutate(t.data, PTM_score = as.double(PTM_score))
    t.data <- dplyr::relocate(t.data, Sequence_name, PTM_type)
    t.data <- dplyr::arrange(t.data, dplyr::desc(PTM_score))
    t.data <- .set_lab(t.data, sig(x), "prediction PTM of", x$type)
    t.cutoff <- dplyr::filter(t.data, PTM_score > .5)
    t.cutoff <- .set_lab(t.cutoff, sig(x), "high score prediction PTM of", x$type)
    x@tables[[ 2 ]] <- namel(t.data, t.cutoff)
    return(x)
  })



