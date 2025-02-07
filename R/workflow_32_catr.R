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
    step_message("Prepare your data with function `job_catr`. ")
  })

setMethod("step1", signature = c(x = "job_catr"),
  function(x, wd = "catRapid", rna_type = c("transcript_exon_intron", "gene_exon_intron")){
    step_message("Prepare upload files.")
    if (is.null(x$mart)) {
      mart <- new_biomart()
      x$mart <- mart
    } else {
      mart <- x$mart
    }
    if (dir.exists(wd)) {
      if (sureThat("dir.exists, remove that?")) {
        unlink(wd, TRUE)
      }
    }
    dir.create(wd, FALSE)
    x$wd <- wd
    protein_seq <- get_seq.pro(object(x)$protein, mart)
    rna_type <- match.arg(rna_type)
    rna_seq <- get_seq.rna(object(x)$rna, mart, to = rna_type)
    x <- methodAdd(x, "以 R 包 `biomaRt` ({packageVersion('biomaRt')}) {cite_show('MappingIdentifDurinc2009')} 获取 RNA 或蛋白质的序列 (`biomaRt::getSequence`)。")
    x <- snapAdd(x, "获取 {bind(object(x)$rna)} 的 {rna_type} 序列，以及 {bind(object(x)$protein)} 的 peptide 序列。")
    write(protein_seq$fasta, "protein")
    write(rna_seq$fasta, "rna")
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

get_seq.pro <- function(ids, mart, unique = TRUE, fasta = TRUE, from = "hgnc_symbol", to = "peptide") {
  ids <- unique(ids)
  data <- e(biomaRt::getSequence(id = ids, type = from, seqType = to, mart = mart))
  data <- filter(data, !grepl("unavailable", !!rlang::sym(to)))
  if (unique) {
    data <- dplyr::mutate(data, long = nchar(!!rlang::sym(to)))
    data <- dplyr::arrange(data, dplyr::desc(long))
    data <- dplyr::distinct(data, hgnc_symbol, .keep_all = TRUE)
  }
  if (fasta) {
    fasta <- apply(data, 1,
      function(vec) {
        c(paste0(">", vec[[2]]), vec[[1]])
      })
    fasta <- .fasta(unlist(fasta))
    namel(data, fasta)
  } else {
    data
  }
}

get_seq.rna <- function(ids, mart, unique = TRUE, fasta = TRUE, from = "hgnc_symbol",
  to = c("transcript_exon_intron", "gene_exon_intron", "coding"))
{
  # mrefseq <- biomaRt::getSequence(id = "NM_001621", type = "refseq_mrna", seqType = "coding", mart = mart)
  # coding <- biomaRt::getSequence(id = "AHR", type = "hgnc_symbol", seqType = "coding", mart = mart)
  # identical(mrefseq[[1]], coding[[1]])
  # ref: <https://www.sangon.com/customerCenter/classOnline/help_gene>
  to <- match.arg(to)
  do.call(get_seq.pro, as.list(environment()))
}


