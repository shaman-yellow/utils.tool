# ==========================================================================
# workflow of Mendelian Randomization (AI) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mr <- setClass("job_mr",
  contains = "job",
  prototype = prototype(
    pg = "mr",
    info = c("https://mrcieu.github.io/TwoSampleMR/"),
    cite = "[@TwoSampleMRHemani2018; @MendelianRandomizationYavorska2023]",
    method = "",
    tag = "mr",
    analysis = "孟德尔随机化"
  )
)

job_mr <- function(exposure, outcome, exposure_name = "Exposure", outcome_name = "Outcome") {
  # exposure/outcome: GWAS ID (character) or data.frame
  if (is.character(exposure)) {
    exposure <- e(TwoSampleMR::extract_instruments(exposure))
  }
  if (is.character(outcome)) {
    outcome <- e(TwoSampleMR::extract_outcome_data(exposure$SNP, outcome))
  }
  x <- .job_mr(object = list(exposure_raw = exposure, outcome_raw = outcome))
  x$exposure_name <- exposure_name
  x$outcome_name <- outcome_name
  x <- methodAdd(x, "暴露数据与结局数据通过 R 包 `TwoSampleMR` 获取并整理。")
  x <- snapAdd(x, "暴露 {exposure_name} 包含 {nrow(exposure)} 个工具变量候选 SNP。")
  x <- snapAdd(x, "结局 {outcome_name} 包含 {nrow(outcome)} 个 SNP。")
  return(x)
}

setMethod("step0", signature = "job_mr",
  function(x) {
    step_message("Prepare data with `job_mr()`.")
  }
)

setMethod("step1", signature = "job_mr",
  function(x, action = 2) {
    step_message("Harmonise exposure and outcome data.")
    stop("AI Generated, not supported now!!!")
    dat <- e(TwoSampleMR::harmonise_data(
      exposure_dat = object(x)$exposure_raw,
      outcome_dat = object(x)$outcome_raw,
      action = action
    ))
    object(x)$harmonised <- dat
    t.harmonised <- set_lab_legend(
      dat,
      glue::glue("{x@sig} Harmonised data"),
      glue::glue("Harmonised 数据表|||对齐暴露与结局的效应等位基因后保留的 SNP 信息，包含效应量、标准误、等位基因频率等。")
    )
    x <- tablesAdd(x, t.harmonised)
    x <- methodAdd(x, "使用 `TwoSampleMR::harmonise_data()` 对齐暴露与结局的效应等位基因 (action = {action})。")
    x <- snapAdd(x, "Harmonise 后保留 {nrow(dat)} 个 SNP 用于 MR 分析。")
    return(x)
  }
)

setMethod("step2", signature = "job_mr",
  function(x, methods = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")) {
    step_message("Perform MR analysis.")
    dat <- object(x)$harmonised
    mr_res <- e(TwoSampleMR::mr(dat, method_list = methods))
    object(x)$mr_res <- mr_res
    t.mr_res <- set_lab_legend(
      mr_res,
      glue::glue("{x@sig} MR results"),
      glue::glue("MR 分析结果表|||各方法的效应估计值 (b)、标准误 (se)、P 值及置信区间。")
    )
    x <- tablesAdd(x, t.mr_res)
    ivw <- dplyr::filter(mr_res, method == "Inverse variance weighted")
    if (nrow(ivw)) {
      or <- exp(ivw$b)
      ci <- exp(ivw$b + c(-1, 1) * 1.96 * ivw$se)
      x <- snapAdd(x, "IVW 方法：OR = {round(or, 2)} (95% CI: {round(ci[1], 2)}–{round(ci[2], 2)})，P = {format(ivw$pval, scientific = TRUE, digits = 3)}。")
    }
    x <- methodAdd(x, "以 `TwoSampleMR::mr()` 执行 MR 分析，使用的方法包括 {bind(methods)}。")
    return(x)
  }
)

setMethod("step3", signature = "job_mr",
  function(x) {
    step_message("Sensitivity analyses.")
    dat <- object(x)$harmonised
    hetero <- e(TwoSampleMR::mr_heterogeneity(dat))
    pleio  <- e(TwoSampleMR::mr_pleiotropy_test(dat))
    leaveone <- e(TwoSampleMR::mr_leaveoneout(dat))
    object(x)$hetero <- hetero
    object(x)$pleio <- pleio
    object(x)$leaveone <- leaveone
    t.hetero <- set_lab_legend(hetero, glue::glue("{x@sig} Heterogeneity tests"), "异质性检验结果表|||Cochran's Q 检验评估 SNP 效应间的异质性。")
    t.pleio  <- set_lab_legend(pleio,  glue::glue("{x@sig} Pleiotropy test"), "水平多效性检验结果表|||MR‑Egger 截距检验评估是否存在方向性多效性。")
    x <- tablesAdd(x, t.hetero, t.pleio)
    het_ivw <- dplyr::filter(hetero, method == "Inverse variance weighted")
    if (nrow(het_ivw)) {
      x <- snapAdd(x, "异质性检验 (IVW) Cochran's Q P = {format(het_ivw$Q_pval, scientific = TRUE)}。")
    }
    if (!is.null(pleio$pval)) {
      x <- snapAdd(x, "多效性检验 (MR‑Egger intercept) P = {format(pleio$pval, scientific = TRUE)}。")
    }
    x <- methodAdd(x, "敏感性分析包括异质性检验 (`mr_heterogeneity`)、多效性检验 (`mr_pleiotropy_test`) 和留一法分析 (`mr_leaveoneout`)。")
    return(x)
  }
)

setMethod("step4", signature = "job_mr",
  function(x) {
    step_message("Generate MR plots.")
    dat <- object(x)$harmonised
    mr_res <- object(x)$mr_res
    p.scatter <- e(TwoSampleMR::mr_scatter_plot(mr_res, dat))[[1]]
    p.scatter <- set_lab_legend(
      wrap(p.scatter, 8, 6),
      glue::glue("{x@sig} MR scatter plot"),
      glue::glue("MR 散点图|||每个点代表一个 SNP，横坐标为 SNP 对暴露的效应，纵坐标为对结局的效应。各方法拟合线叠加展示。")
    )
    p.forest <- e(TwoSampleMR::mr_forest_plot(object(x)$leaveone))[[1]]
    p.forest <- set_lab_legend(
      wrap(p.forest, 8, 6),
      glue::glue("{x@sig} Leave‑one‑out forest plot"),
      glue::glue("留一法森林图|||逐次剔除一个 SNP 后重新计算 IVW 效应，用于评估单个 SNP 对整体结果的影响。")
    )
    p.funnel <- e(TwoSampleMR::mr_funnel_plot(object(x)$leaveone))[[1]]
    p.funnel <- set_lab_legend(
      wrap(p.funnel, 8, 6),
      glue::glue("{x@sig} Funnel plot"),
      glue::glue("漏斗图|||横坐标为 SNP 效应估计值，纵坐标为其标准误的倒数，对称性越好提示偏倚越小。")
    )
    x <- plotsAdd(x, p.scatter, p.forest, p.funnel)
    x <- snapAdd(x, "MR 可视化结果：散点图{aref(p.scatter)}、留一法森林图{aref(p.forest)}及漏斗图{aref(p.funnel)}。")
    return(x)
  }
)
