# ==========================================================================
# function for getting published data (tables)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

get_data.mrlx2022 <- function(file = .prefix("published_data/MendelianRandoLiuX2022_s1.xlsx"))
{
  fun <- function(data) {
    data <- dplyr::mutate(data,
      variant_id = paste0(Chr., "_", BP, "_", A1, "_", A2, "_b38")
    )
    data <- dplyr::relocate(data, variant_id)
    data <- dplyr::select(data, variant_id:Pdiscovery)
    data
  }
  snp_microbiota <- fxlsx(file, sheet = 3, startRow = 4)
  snp_microbiota <- fun(snp_microbiota)
  get_pattern <- function(x) {
    x <- ifelse(grpl(x, "^[a-z]_"), gs(x, "^[a-z]_", ""), NA)
    ifelse(is.na(x), NA, gs(x, "_", " "))
  }
  snp_microbiota <- dplyr::mutate(snp_microbiota,
    .gutmd_pattern = get_pattern(Microbiome.features)
  )
  snp_metabolite <- fxlsx(file, sheet = 7, startRow = 4)
  snp_metabolite <- fun(snp_metabolite)
  observational_correlation <- fxlsx(file, sheet = 10, startRow = 3)
  fun <- function(data) {
    which <- which(is.na(data$X2))
    x1 <- dplyr::mutate(data[ (which[1] + 1):(which[2] - 1), ], type = data$X1[which[1]])
    x2 <- dplyr::mutate(data[ (which[2] + 1):nrow(data), ], type = data$X1[which[2]])
    x2[, 1:2] <- x2[, 2:1]
    dplyr::bind_rows(x1, x2)
  }
  one_sample_mr <- fun(fxlsx(file, sheet = 11, startRow = 4))
  two_sample_mr <- fun(fxlsx(file, sheet = 12, startRow = 4))
  lst <- namel(snp_microbiota, snp_metabolite, observational_correlation, one_sample_mr, two_sample_mr)
  lst[1:2] <- lapply(lst[1:2], function(x) dplyr::select(x, -Gene))
  x <- .job_publish(object = lst, cite = "[@MendelianRandoLiuX2022]")
  x$get_pattern <- get_pattern
  x
}

get_data.aaog2014 <- function(file = .prefix("published_data/AnAtlasOfGenShin2014_s1.xlsx")) {
  .job_publish(object, cite = "[@AnAtlasOfGenShin2014]")
}

get_data.cacc2021 <- function(file = .prefix("published_data/ChangesAndCorChen2021_s2.xlsx"))
{
  data <- fxlsx(file)
  data <- dplyr::select(data, metabolite = 1, microbiota = 2, cor = 3, pvalue = 4, AdjPvalue)
  data <- add_anno(.corp(data))
  x <- .job_publish(cite = "[@ChangesAndCorChen2021]")
  x$heatmap <- wrap(callheatmap(new_heatdata(.corp(data))), 28, 12)
  x$heatmap <- .set_lab(x$heatmap, "PUBLISHED-ChangesAndCorChen2021-correlation-heatmap")
  object(x) <- dplyr::filter(data, pvalue < .05)
  object(x) <- .set_lab(object(x), "PUBLISHED-ChangesAndCorChen2021-significant-correlation")
  ids <- PubChemR::get_cids(unique(object(x)$metabolite))
  ids <- e(PubChemR::CIDs(ids))
  ids <- dplyr::mutate(ids, CID = as.integer(CID))
  object(x) <- map(object(x), "metabolite", ids, "Identifier", "CID", col = "cid")
  x
}

get_data.pmb2023 <- function(file = .prefix("published_data/ProteinMetabolBenson2023_mmc5_small.rds"))
{
  # https://github.com/aeisman/protein-metabolite
  # https://mbenson.shinyapps.io/protein-metabolite/
  data <- readRDS(file)
  data <- dplyr::filter(data, META_Q < .05)
  x <- .job_publish(object = data, cite = "[@ProteinMetabolBenson2023]")
  object(x) <- .set_lab(object(x), "PUBLISHED-ProteinMetabolBenson2023-significant-correlation")
  x
}

get_data.pps2022 <- function(dir = .prefix("published_data/ProteomicsProfShao2022"))
{
  metadata <- fxlsx(file.path(dir, "ProteomicsProfShao2022_s4.xlsx"))
  metadata <- .set_lab(metadata, "PUBLISHED-ProteomicsProfShao2022-metadata")
  data <- fxlsx(file.path(dir, "ProteomicsProfShao2022_s1.xlsx"), startRow = 2)
  colnames(data)[1:3] <- c("UniProt_Knowledgebase_ID", "Protein_name", "Gene_name")
  data <- .set_lab(data, "PUBLISHED-ProteomicsProfShao2022-data")
  obj <- namel(metadata, data)
  x <- .job_publish(object = obj, cite = "[@ProteomicsProfShao2022]")
  # The MS proteomics data are available on the iProX database with the project
  # ID: IPX0001414000 and the subproject ID: IPX000141400. The raw sequence
  # data have been uploaded to SRA with an ID: PRJNA598559.
  x$metadata <- dplyr::rename(metadata, sample = Patient.ID)
  x$counts <- dplyr::select(data, 1, tidyselect::matches("^P[0-9]+$"))
  message("The `x$counts` is log2 transformed from Raw data.")
  x$counts <- dplyr::mutate_if(x$counts, is.double, log2)
  return(x)
}

get_data.rot2016 <- function(file = .prefix("published_data/RefinementOfTLehman2016_s3.xlsx")) {
  # https://pubmed.ncbi.nlm.nih.gov/27310713/
  obj <- fxlsx(file, sheet = 2, startRow = 2)
  x <- .job(object = obj, cite = "[@RefinementOfTLehman2016]",
    method = "Selection of TNBC cases in TCGA-BRCA was adopted from results of study of Lehmann et al")
  return(x)
}

get_data.nmt2015 <- function(file = .prefix("published_data/NovelMarkersTJablon2015.rds"), update = FALSE)
{
  # scRNA-seq, mouse Macrophage, M0, M1, M2 signature
  if (file.exists(file) && !update) {
    obj <- readRDS(file)
  } else {
    # Top 10 up- and down-Regulated genes common to M1 and M2 macrophages.
    m1m2 <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4689374/table/pone.0145342.t003/?report=objectonly"
    # Genes increased more than 2-FC in M1 and decreased more than 2-FC in M2 macrophages.
    m1_up_m2_down <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4689374/table/pone.0145342.t004/?report=objectonly"
    # Genes increased more than 2-FC in M2 and decreased more than 2-FC in M1 macrophages.
    m2_up_m1_down <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4689374/table/pone.0145342.t005/?report=objectonly"
    fun_m1m2 <- function(url) {
      x <- as_tibble(get_table.html(RCurl::getURL(url))[[1]])
      x <- dplyr::filter(x, !grpl(V1, "Gene|Regulated"), grpl(V1, "^[A-Za-z0-9]+$"))
      colnames(x) <- c("gene", "FC.M1_vs_M0", "p.M1_vs_M0", "FC.M2_vs_M0", "p.M2_vs_M0")
      x <- dplyr::mutate(x, dplyr::across(2:5, as.double))
      x
    }
    fun_m <- function(url) {
      x <- as_tibble(get_table.html(RCurl::getURL(url))[[1]], .name_repair = "minimal")
      colnames(x) <- c("gene", "FC.M1_vs_M0", "p.M1_vs_M0", "FC.M2_vs_M0", "p.M2_vs_M0")
      x <- dplyr::mutate(x, dplyr::across(2:5, as.double))
      x
    }
    sig.m0 <- fun_m1m2(m1m2)
    sig.m1 <- fun_m(m1_up_m2_down)
    sig.m2 <- fun_m(m2_up_m1_down)
    obj <- namel(sig.m0, sig.m1, sig.m2)
    saveRDS(obj, file)
  }
  x <- .job(object = obj, cite = "[@NovelMarkersTJablon2015]",
    method = "The data in published article of Jablonski et al used for distinguishing macrophage phenotypes (M0/M1/M2)"
  )
  return(x)
}
