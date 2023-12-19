# ==========================================================================
# function for getting published data (tables)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

get_data.mrlx2022 <- function(file = "~/outline/lixiao/published_data/MendelianRandoLiuX2022_s1.xlsx") {
  fun <- function(data) {
    data <- dplyr::mutate(data,
      variant_id = paste0(Chr., "_", BP, "_", A1, "_", A2, "_b38")
    )
  }
  snp_microbiota <- fxlsx(file, sheet = 3, startRow = 4)
  snp_microbiota <- fun(snp_microbiota)
  snp_metabolite <- fxlsx(file, sheet = 7, startRow = 4)
  snp_metabolite <- fun(snp_metabolite)
  lst <- namel(snp_microbiota, snp_metabolite)
  .job_publish(object = lst, cite = "[@MendelianRandoLiuX2022]")
}

get_data.aaog2014 <- function(data = "~/outline/lixiao/published_data/AnAtlasOfGenShin2014_s1.xlsx") {
  .job_publish(object, cite = "[@AnAtlasOfGenShin2014]")
}



