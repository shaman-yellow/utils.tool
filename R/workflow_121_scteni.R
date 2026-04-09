# ==========================================================================
# workflow of scteni
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scteni <- setClass("job_scteni", 
  contains = c("job"),
  prototype = prototype(
    pg = "scteni",
    info = c("https://github.com/cailab-tamu/scTenifoldKnk"),
    cite = "",
    method = "",
    tag = "scteni",
    analysis = "scTenifoldKnk 虚拟敲除"
    ))


setGeneric("asjob_scteni",
  function(x, ...) standardGeneric("asjob_scteni"))

setMethod("asjob_scteni", signature = c(x = "ANY"),
  function(x, ref){
    fea <- resolve_feature_snapAdd_onExit("x", ref)
    x <- .job_scteni(object = fea)
    return(x)
  })

setMethod("step0", signature = c(x = "job_scteni"),
  function(x){
    step_message("Prepare your data with function `job_scteni`.")
  })

setMethod("step1", signature = c(x = "job_scteni"),
  function(x){
    step_message("Do nothing")
    x <- methodAdd(x, "**scTenifoldKnk** 是一种基于单细胞转录组数据构建基因调控网络并进行 *in silico* 虚拟敲除的计算方法，其核心目的是在无需实际基因编辑实验的情况下，评估特定基因在细胞系统中的功能重要性及其对全局转录调控网络的影响。具体而言，该方法通过对目标基因进行网络层面的“移除”，并比较敲除前后基因调控网络结构及表达模式的变化，从而识别潜在的下游调控通路及受影响的关键基因模块。该分析有助于优先筛选具有重要生物学意义的候选关键基因，为后续机制研究和实验验证提供理论依据。")
    x <- methodAdd(x, "以 R 包 `scTenifoldKnk` ⟦pkgInfo('scTenifoldKnk')⟧ 模拟敲除分析，以推测关键基因下游功能影响。构建好的基因调控网络（gene regulatory network，GRN）中进行无监督虚拟敲低，通过模拟关键基因节点的缺失来模拟其敲低状态，进而筛选出因关键基因敲低而调控关系发生显著改变的基因。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_scteni"),
  function(x, use.p = c("p.adj", "p.value"), cut.p = .05, 
    cut.z = 2, recode = NULL, dir = "sct_results",
    lst_diff = .read_scteni_results.hb(dir))
  {
    step_message("Load results file and draw volcano plot.")
    use.p <- match.arg(use.p)
    lst_diff <- sapply(names(lst_diff), simplify = FALSE,
      function(name) {
        data <- lst_diff[[name]]
        data <- dplyr::filter(
          data, gene != !!name, !!rlang::sym(use.p) < !!cut.p
        )
        if (!is.null(cut.z)) {
          data <- dplyr::filter(data, abs(Z) > !!cut.z)
        }
        if (!is.null(recode)) {
          data <- dplyr::mutate(
            data, gene = dplyr::recode(gene, !!!recode)
          )
        }
        data
      })
    if (!is.null(recode)) {
      names(lst_diff) <- dplyr::recode(names(lst_diff), !!!recode)
    }
    lst_diff <- set_lab_legend(
      lst_diff,
      glue::glue("{x@sig} data of genes affected by {names(lst_diff)} knockout"),
      glue::glue("受 {names(lst_diff)} 敲除影响最显著的基因")
    )
    x <- tablesAdd(x, t.all_diff = lst_diff)
    feature(x) <- as_feature(lapply(lst_diff, function(x) x$gene), "受关键基因敲除而显著影响的基因")
    ps.volcano <- sapply(names(lst_diff), simplify = FALSE,
      function(name) {
        data <- lst_diff[[name]]
        cut.fc <- if (is.null(cut.z)) 0 else cut.z
        p <- plot_volcano(
          data, "gene", use = use.p, use.fc = "Z",
          fc = cut.fc, label.fc = "Z-score", f.nudge = .1,
          mode_fc = 1, HLs = head(data$gene, n = 10), show_legend = FALSE
        )
        set_lab_legend(
          wrap(p, 5, 6),
          glue::glue("{x@sig} genes affected by {name} knockout"),
          glue::glue("受 {name} 敲除影响最显著的基因|||横坐标为标准化Z分数，Z 的绝对值越大表示该基因受 KO 扰动越显著；纵坐标为下游基因。")
        )
      })
    x <- plotsAdd(x, ps.volcano)
    snap <- .stat_table_by_pvalue(
      dplyr::bind_rows(lst_diff, .id = "ID"),
      n = 10, split = "ID", use.p = use.p, colName = "gene", 
      target = "基因", by = "被敲除后显著影响到", needSum = FALSE
    )
    x <- snapAdd(
      x, "关键基因敲除导致调控关系发生变化，{snap}"
    )
    if (!is.null(cut.z)) {
      snap.z <- glue::glue("& |Z| &gt; {cut.z}")
    } else {
      snap.z <- ""
    }
    x <- methodAdd(x, "筛选因关键基因敲低而调控关系发生显著改变的基因（阈值为{use.p} &lt; {cut.p} {snap.z}）。")
    return(x)
  })

.read_scteni_results.hb <- function(dir, type = "3000HVG",
  pattern = paste0("scTenifoldKnk_", type, "_"))
{
  files <- list.files(dir, pattern, full.names = TRUE)
  res <- lapply(files,
    function(file) {
      readRDS(file)$diffRegulation
    })
  setNames(
    res, s(
      tools::file_path_sans_ext(basename(files)), pattern, ""
    )
  )
}
