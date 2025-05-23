# ==========================================================================
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

set.sig.wd <- function(gse)
  {
    path.work <- .prefix("geo_db/ahr_sig", "db")
    path.de.1 <- "/ftp.ncbi.nlm.nih.gov/geo/series/"
    gse.dir <- gsub("[0-9]{3}$", "nnn", gse)
    wd <- paste0(path.work, path.de.1, gse.dir, "/", gse, "/suppl")
    setwd(wd)
    cat("## The work directory at:", getwd(), "\n")
    print(list.files(all.files = TRUE))
    cat("## Done\n")
  }
 
set.initial.wd <- function(wd = .prefix("geo_db/ahr_sig", "db"))
  {
    setwd(wd)
  }

decomp_tar2txt <- function(){
    list.files(pattern = "_RAW\\.tar$") %>%
      utils::untar(exdir = ".")
    list.files(pattern = "\\.gz$") %>%
      lapply(R.utils::gunzip)
    file <- list.files(pattern = "^GSM.*txt$")
    df <- data.table::data.table(file = file)
    return(df)
  }
 
anno.gene.biomart <- 
  function(
    db = "hsapiens_gene_ensembl",
    host = NULL,
    attr = NA,
    ex.attr = NA
    )
  {
    if(is.null(host)){
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = db)
    }else{
      ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = db, host = host)
    }
    if(!is.character(attr)){
      attr <- c(
        "ensembl_gene_id",
        "entrezgene_id",
        "hgnc_symbol",
        "chromosome_name",
        "start_position",
        "end_position",
        "description"
      )
    }
    if(is.character(ex.attr))
      attr <- c(attr, ex.attr)
    anno <- biomaRt::getBM(attr, mart = ensembl)
    anno <- dplyr::as_tibble(anno)
    return(anno)
  }
 
anno.gene.edb <- 
  function(keys,
    package = "EnsDb.Hsapiens.v86",
    columns = c("GENEID", "SYMBOL", "EXONID",
      "EXONSEQSTART", "EXONSEQEND"),
    keytype = "GENEID",
    ex.attr = NA
    )
  {
    if(is.character(ex.attr))
      columns <- c(columns, ex.attr)
    obj <- paste0(package, "::", package)
    obj <- eval(parse(text = obj))
    genes <- AnnotationDbi::select(obj,
      keytype = keytype,
      keys = keys,
      columns = columns)
    genes <- dplyr::distinct(
      genes, dplyr::across(dplyr::all_of(keytype)), .keep_all = TRUE
    )
  }

anno.into.list <- 
  function(dge = dge.list, anno = gene.anno, col = "ensembl_gene_id")
  {
    genes <- data.frame(col = rownames(dge))
    colnames(genes) <- col
    genes <- merge(genes, anno, all.x = TRUE, by = col, sort = FALSE)
    genes <- dplyr::distinct(genes, dplyr::across(dplyr::all_of(col)), .keep_all = TRUE)
    cat("## is.na:\n")
    check.col <- colnames(genes)[2]
    check.na <- dplyr::filter(
      genes, is.na(!!!rlang::syms(check.col))) %>% 
      dplyr::as_tibble()
    print(check.na)
    dge$genes <- genes
    return(dge)
  }

re.sample.group <- 
  function(dge = dge.list, meta = meta.df)
  {
    colnames(dge) <- meta$sample
    dge$samples$group <- meta$group
    return(dge)
  }
# fpkm2count <- 
# function(
#          dge = dge.list,
#          eff.len = "eff.len",
#          lim.min = 0.0001
#          ){
#   num <- nrow(dge$counts)
#   eff.len <- dge$genes[[eff.len]]
#   ## counts...
#   cset <- apply(dge$counts, 2,
#                 function(vec){
#                   vec <- log(vec + lim.min, base = exp(1))
#                   vec <- vec - log(1e9) + log(eff.len)
#                   sum.counts.delog <- sum(vec) / (1 - num)
#                   vec <- vec + sum.counts.delog
#                   round(10 ^ vec)
#                 })
#   dge$counts <- cset
#   return(dge)
# }

fpkm_log2tpm <- 
  function(dge = dge.list)
  {
    cset <- apply(dge$counts, 2,
      function(vec){
        log2(
          # vec / sum(vec) * 1e6 + 1
          exp(log(vec) - log(sum(vec, na.rm = TRUE)) + log(1e6)) + 1
        )
      })
    dge$counts <- cset
    return(dge)
  }

get_gsm.data <- 
  function(info)
  {
    lapply(info, function(obj){
      gsm <- Biobase::phenoData(obj) %>% 
        Biobase::sampleNames() %>% 
        unlist() %>% 
        pbapply::pblapply(function(gsm){
          gsm.dir <- gsub("[0-9]{3}$", "nnn", gsm)
          ftp <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/",
            gsm.dir, "/", gsm, "/suppl/")
          system(paste("wget -np -m", ftp))
        })
    })
  }

# limma_downstream.eset <- 
#   function(
#     eset,
#     design,
#     contr.matrix,
#     cut.q = 0.05,
#     cut.fc = 0.3
#     ){
#     fit <- limma::lmFit(eset, design)
#     fit.cont <- limma::contrasts.fit(fit, contrasts = contr.matrix)
#     ebayes <- limma::eBayes(fit.cont)
#     res <- lapply(seq_len(ncol(contr.matrix)), function(coef){
#       results <- limma::topTable(ebayes, coef = coef, number = Inf) %>% 
#         dplyr::filter(adj.P.Val < cut.q, abs(logFC) > cut.fc) %>% 
#         dplyr::as_tibble() 
#       return(results)
#     })
#     names(res) <- colnames(contr.matrix)
#     return(res)
#   }

list.attr.biomart <- 
  function(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  {
    attr <- biomaRt::useEnsembl(biomart = biomart,
      dataset = dataset) %>% 
    biomaRt::listAttributes(mart = .) %>% 
    dplyr::as_tibble()
  return(attr)
  }

show_boxplot <- 
  function(df)
  {
    df <- dplyr::as_tibble(df) %>% 
      dplyr::mutate(id = rownames(df)) %>% 
      reshape2::melt(id.vars = "id", variable.name = "class", value.name = "value")
    p <- ggplot(df) +
      geom_boxplot(aes(x = class, y = value))
    p 
  }
