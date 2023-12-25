# ==========================================================================
# function for getting published data (tables)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

get_data.mrlx2022 <- function(file = "~/outline/lixiao/published_data/MendelianRandoLiuX2022_s1.xlsx")
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
  lst <- namel(snp_microbiota, snp_metabolite)
  lst <- lapply(lst, function(x) dplyr::select(x, -Gene))
  x <- .job_publish(object = lst, cite = "[@MendelianRandoLiuX2022]")
  x$get_pattern <- get_pattern
  x
}

get_data.aaog2014 <- function(file = "~/outline/lixiao/published_data/AnAtlasOfGenShin2014_s1.xlsx") {
  .job_publish(object, cite = "[@AnAtlasOfGenShin2014]")
}

get_data.cacc2021 <- function(file = "~/outline/lixiao/published_data/ChangesAndCorChen2021_s2.xlsx") {
  data <- fxlsx(file)
  data <- dplyr::select(data, metabolite = 1, microbiota = 2, cor = 3, pvalue = 4, AdjPvalue)
  data <- add_anno(.corp(data))
  x <- .job_publish(cite = "[@ChangesAndCorChen2021]")
  x$heatmap <- wrap(callheatmap(new_heatdata(.corp(data))), 28, 12)
  x$heatmap <- .set_lab(x$heatmap, "PUBLISHED-ChangesAndCorChen2021-correlation-heatmap")
  object(x) <- dplyr::filter(data, pvalue < .05)
  object(x) <- .set_lab(object(x), "PUBLISHED-ChangesAndCorChen2021-significant-correlation")
  x
}

