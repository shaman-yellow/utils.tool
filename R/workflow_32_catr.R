# ==========================================================================
# workflow of catr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_catr <- setClass("job_catr", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://service.tartaglialab.com/static_files/shared/documentation_omics2.html"),
    cite = "[@ICatIRapidArmaos2021]",
    method = "The `catRAPID omics` v2.1 used for protein binding with RNA prediction.",
    tag = "rbp",
    analysis = "CatRAPID omics RBP 预测"
    ))

job_catr <- function(protein, rna)
{
  .job_catr(object = list(protein = gname(protein), rna = gname(rna)))
}

setMethod("step0", signature = c(x = "job_catr"),
  function(x){
    step_message("Prepare your data with function `job_catr`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_catr"),
  function(x, wd = timeName("catRapid"), rna_type = "gene_exon_intron"){
    step_message("Prepare upload files.")
    if (is.null(x$mart)) {
      mart <- new_biomart()
      x$mart <- mart
    } else {
      mart <- x$mart
    }
    dir.create(wd, FALSE)
    x$wd <- wd
    protein_seq <- get_seq.pro(object(x)$protein, mart)
    rna_seq <- get_seq.rna(object(x)$rna, mart, to = rna_type)
    write(protein_seq$fasta, paste0(wd, "/protein"))
    write(rna_seq$fasta, paste0(wd, "/rna"))
    x$protein_seq <- protein_seq
    x$rna_seq <- rna_seq
    return(x)
  })

setMethod("step2", signature = c(x = "job_catr"),
  function(x, num = 1, dir = "~/Downloads/", pattern = "output_full.*\\.zip"){
    step_message("Deparse the results files. red{{Make sure the results files exists in the `dir`.}}")
    collateFiles(paste0("candidate_", num), pattern, from = dir, to = x$wd, suffix = ".zip")
    res <- lapply(list.files(x$wd, "\\.zip$", full.names = TRUE),
      function(file) {
        ftibble(unzip(file, exdir = x$wd))
      })
    res <- do.call(dplyr::bind_rows, res)
    ## filter
    top <- dplyr::arrange(res, Protein_ID, RNA_ID, dplyr::desc(Ranking))
    top <- dplyr::distinct(top, Protein_ID, RNA_ID, .keep_all = TRUE)
    data <- dplyr::mutate(top, .id = seq_len(nrow(top)))
    args <- list(
      rlang::quo(RBP_Propensity == 1),
      rlang::quo(Interaction_Propensity > 0),
      rlang::quo(numof.RNA.Binding_Domains_Instances > 0),
      rlang::quo(numof.RNA_Binding_Motifs_Instances > 0)
    )
    sets <- lapply(args,
      function(arg) {
        dplyr::filter(data, !!arg)$.id
      })
    names(sets) <- c(
      "RBP_Propensity",
      "Interaction_Propensity",
      "numof.RNA.Binding_Domains_Instances",
      "numof.RNA_Binding_Motifs_Instances"
    )
    p.upset <- new_upset(lst = sets, trunc = NULL)
    x@tables[[ 2 ]] <- list(res = res, top = top)
    x@plots[[ 2 ]] <- list(p.upset = p.upset)
    x$sets_filtered <- sets
    x$sets_ins <- ins(lst = sets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_catr"),
  function(x, ref = NULL, top = NULL, group = NULL, group.use = 1, label.auto = TRUE, label.freq = NULL)
  {
    step_message("Select data for Sankey plot.")
    if (!is.null(ref)) {
      dat.ref <- dplyr::slice(x@tables$step2$top, ref)
    } else {
      dat.ref <- x@tables$step2$top
    }
    if (!is.null(top)) {
      dat.ref <- dplyr::arrange(dat.ref, dplyr::desc(Ranking))
      dat.ref <- head(dat.ref, n = top)
    }
    ## p.alluvial
    data <- dplyr::select(dat.ref, 1:2)
    if (!is.null(group)) {
      names <- rep(names(group), lengths(group))
      values <- unlist(group, use.names = FALSE)
      data$group <- names[match(data[[ group.use ]], values)]
    } else {
      data$group <- "group1"
    }
    p.allu <- new_allu(data, 3, label.auto = label.auto, label.freq = label.freq)
    x@tables[[ 3 ]] <- namel(dat.ref)
    x@plots[[ 3 ]] <- namel(p.allu)
    return(x)
  })

timeName <- function(prefix) {
  paste0(prefix, gs(Sys.time(), " |:", "_"))
}
