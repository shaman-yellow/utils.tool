# ==========================================================================
# workflow of fe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_fe <- setClass("job_fe", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "fe",
    info = c("http://www.zhounan.org/ferrdb/current/"),
    cite = "[@FerrdbV2UpdaZhou2023]",
    method = "",
    tag = "fe",
    analysis = "FerrDb 铁死亡调控因子"
    ))

job_fe <- function()
{
  .job_fe()
}

setMethod("step0", signature = c(x = "job_fe"),
  function(x){
    step_message("Prepare your data with function `job_fe`.")
  })

setMethod("step1", signature = c(x = "job_fe"),
  function(x, ...){
    step_message("Get Ferroptosis data.")
    t.ferroptosisRegulators <- get_fe_data(...)
    p.ferroptosisRegulatorsDistribution <- new_pie(
      rep(names(t.ferroptosisRegulators),
        vapply(t.ferroptosisRegulators, function(x) length(unique(x$symbol)), integer(1)))
    )
    x <- plotsAdd(x, p.ferroptosisRegulatorsDistribution)
    x <- tablesAdd(x, t.ferroptosisRegulators)
    x <- methodAdd(x, "从数据库 `FerrDb V2` {cite_show('FerrdbV2UpdaZhou2023')} 获取与铁死亡相关的调控因子或铁死亡与疾病之间的关联信息 <http://www.zhounan.org/ferrdb/current/>。")
    # a list
    s.com <- try_snap(t.ferroptosisRegulators, "symbol")
    x <- snapAdd(x, "铁死亡相关调控因子统计：{s.com}")
    return(x)
  })

setMethod("map", signature = c(x = "job_fe", ref = "list"),
  function(x, ref, use = c("all", "sep")){
    use <- match.arg(use)
    if (is.null(names(ref)) && length(ref) == 1) {
      names(ref) <- "Related_targets"
    }
    ref <- lapply(ref, unique)
    alls <- lapply(x@tables$step1$t.ferroptosisRegulators, function(x) unique(x$symbol))
    if (use == "all") {
      alls <- unique(unlist(alls))
      lst <- c(ref, list(Ferroptosis_all = alls))
    } else {
      names(alls) <- paste0("Ferroptosis_", names(alls))
      lst <- c(ref, alls)
    }
    x$upset <- new_venn(lst = lst)
    x$.map_heading <- "FerrDb 与铁死亡相关基因的交集"
    return(x)
  })

get_fe_data <- function(use.symbol = T, for_gsea = F,
  path = .prefix("../data/ferroptosis_2024-12-05_small.rds"),
  add_internal_job = F, collate = F)
{
  if (collate) {
    # <http://www.zhounan.org/ferrdb/current/>
    fe_db <- list(marker = "~/Downloads/ferroptosis_marker.csv",
      driver = "~/Downloads/ferroptosis_driver.csv",
      suppressor = "~/Downloads/ferroptosis_suppressor.csv",
      unclassifier = "~/Downloads/ferroptosis_unclassified.csv",
      inducer = "~/Downloads/ferroptosis_inducer.csv",
      inhibitor = "~/Downloads/ferroptosis_inhibitor.csv",
      disease = "~/Downloads/ferroptosis_disease.csv"
    )
    fe_db <- lapply(fe_db, ftibble)
    saveRDS(fe_db, path)
  }
  data <- readRDS(path)
  if (use.symbol) {
    data <- lapply(data,
      function(x) {
        if (any(colnames(x) == "symbol")) {
          return(x)
        }
      })
    data <- lst_clear0(data)
  }
  if (for_gsea) {
    gsea <- data.table::rbindlist(data, idcol = T, fill = T)
    gsea <- dplyr::mutate(gsea, term = paste0("Ferroptosis_", .id))
    gsea <- as_tibble(dplyr::relocate(gsea, term, symbol))
    return(gsea)
  }
  if (add_internal_job) {
    job <- .job(method = "Database of `FerrDb V2` used for obtaining ferroptosis regulators",
      cite = "[@FerrdbV2UpdaZhou2023]")
    .add_internal_job(job)
  }
  data <- .set_lab(data, "Ferroptosis regulators", names(data))
  lab(data) <- "Ferroptosis regulators"
  data
}
