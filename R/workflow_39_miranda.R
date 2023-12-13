# ==========================================================================
# workflow of miranda
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_miranda <- setClass("job_miranda", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@Enrigh2003; @MirbaseFromMKozoma2019]",
    method = paste0("The tool of `miRanda` used for predicting mRNA targets for microRNAs (miRNA) ",
      "and `miRBase` used for getting sequence of miRNA")
    ))

job_miranda <- function(mirna, rna)
{
  .job_miranda(object = list(mirna = gname(mirna), rna = gname(rna)))
}

setMethod("step0", signature = c(x = "job_miranda"),
  function(x){
    step_message("Prepare your data with function `job_miranda`."
    )
  })

setMethod("step1", signature = c(x = "job_miranda"),
  function(x, wd = timeName("miranda"), rna_type = "gene_exon_intron"){
    step_message("Prepare computional files.")
    if (is.null(x$mart)) {
      mart <- new_biomart()
      x$mart <- mart
    } else {
      mart <- x$mart
    }
    dir.create(wd, F)
    x$wd <- wd
    # mirna_seq <- get_seq.rna(object(x)$mirna, mart)
    # mirna_file <- write(mirna_seq$fasta, paste0(wd, "/mirna"))[[1]]
    rna_seq <- get_seq.rna(object(x)$rna, mart, to = rna_type)
    rna_file <- write(rna_seq$fasta, paste0(wd, "/rna"))[[1]]
    # x$mirna_file <- mirna_file
    # x$mirna_seq <- mirna_seq
    x$rna_file <- rna_file
    x$rna_seq <- rna_seq
    return(x)
  })
