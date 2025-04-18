# ==========================================================================
# workflow of musite
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_musite <- setClass("job_musite", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "musitePython",
    info = c("https://github.com/duolinwang/MusiteDeep_web"),
    cite = "[@MusitedeepADWang2020]",
    method = "Python tool `MusiteDeep` was used for protein post-translational modification site prediction and visualization",
    tag = "musite",
    analysis = "MusiteDeep 蛋白质转录后修饰位点预测"
    ))

setGeneric("asjob_musite", group = list("asjob_series"),
   function(x, ...) standardGeneric("asjob_musite"))

job_musite <- function(symbols)
{
  symbols <- resolve_feature_snapAdd_onExit("x", symbols, fun_generator = asjob_musite)
  x <- .job_musite(object = symbols)
  return(x)
}

setMethod("step0", signature = c(x = "job_musite"),
  function(x){
    step_message("Prepare your data with function `job_musite`.")
  })

setMethod("step1", signature = c(x = "job_musite"),
  function(x, org = c("hsa", "mmu"), save = "Seq")
  {
    step_message("Prepare protein sequences.")
    org <- match.arg(org)
    x$mart <- new_biomart(org)
    x$seqs <- get_seq.pro(object(x), x$mart)
    x$seqs_file <- write(x$seqs$fasta, save, max = NULL)
    x <- snapAdd(x, "以 `biomaRt` 获取蛋白质 ({less(object(x))}) 的序列 (Peptide)。")
    x <- methodAdd(x, "以 `biomaRt` {cite_show('MappingIdentifDurinc2009')} 获取蛋白质 ({org}) 的序列 (`biomaRt::getSequence` 获取 peptide)。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_musite"),
  function(x,
    type = c("Hydroxylysine", "Hydroxyproline", "Methylarginine",
      "Methyllysine", "N-linked_glycosylation", "N6-acetyllysine",
      "O-linked_glycosylation", "Phosphoserine_Phosphothreonine",
      "Phosphotyrosine", "Pyrrolidone_carboxylic_acid", "S-palmitoyl_cysteine",
      "SUMOylation", "Ubiquitination"
      ), cutoff = .5, respective = TRUE,
    output = timeName("musite"))
  {
    step_message("Compute PTM position.")
    if (is.null(x$output)) {
      x$output <- output
      dir.create(x$output)
    }
    x$type <- match.arg(type, several.ok = TRUE)
    message("Use PTM type: ", bind(x$type))
    models <- paste0(pg("musiteModel"), "/", paste0(x$type, collapse = ";"))
    message("Use models: ", models)
    output_file <- paste0(x$output, "/res_results.txt")
    if (!file.exists(output_file)) {
      cdRun(pg(x), " ", pg("musitePTM"),
        " -input ", x$seqs_file,
        " -output ", x$output, "/res",
        " -model-prefix \"", models, "\""
      )
    }
    lst <- sep_list(readLines(output_file), "^>", TRUE)
    t.data <- lapply(lst[-1],
      function(x) {
        data.table::fread(text = x[-1])
      })
    names(t.data) <- gs(vapply(lst[-1], function(x) x[1], character(1)), "^>", "")
    t.data <- frbind(t.data, idcol = "Sequence_name")
    colnames(t.data)[-1] <- c("Position",	"Residue", "PTMscores", "Cutoff")
    t.data <- dplyr::select(t.data, -Cutoff)
    t.data <- tidyr::separate(t.data, PTMscores, c("PTM_type", "PTM_score"), ":")
    t.data <- dplyr::mutate(t.data, PTM_score = as.double(PTM_score))
    t.data <- dplyr::relocate(t.data, Sequence_name, PTM_type)
    t.data <- dplyr::arrange(t.data, dplyr::desc(PTM_score))
    t.data <- .set_lab(t.data, sig(x), "prediction PTM of", x$type)
    t.data <- setLegend(t.data, "以 `MusiteDeep` 预测的修饰位点以及得分附表。")
    t.cutoff <- dplyr::filter(t.data, PTM_score > cutoff)
    t.cutoff <- .set_lab(t.cutoff, sig(x), "high score prediction PTM of", x$type)
    t.cutoff <- setLegend(t.cutoff, "预测的修饰位点以及得分 (Score cutoff: {cutoff}) 附表。 ")
    p.PTMsNumber <- lapply(split(t.cutoff, ~ PTM_type),
      function(data) {
        new_pie(data$Sequence_name)
    })
    p.PTMsNumber <- .set_lab(p.PTMsNumber, sig(x), names(p.PTMsNumber), "PTM numbers")
    p.PTMsNumber <- setLegend(p.PTMsNumber, "为预测到的 PTMs 数量饼图。")
    p.hasPTMs <- new_pie(object(x) %in% unique(t.cutoff$Sequence_name))
    strFasta <- Read(x$seqs$fasta)
    fun_str <- function(name, pos) {
      mapply(FUN = function(n, p) {
        substr(strFasta[[ n ]], p-2, p+2)
      }, name, pos)
    }
    data <- dplyr::mutate(t.data, str_around = fun_str(Sequence_name, Position))
    data <- tidyr::separate(data, str_around, c(NA, n(pos, 5)), "")
    data <- tidyr::pivot_longer(data,
      tidyselect::starts_with("pos", FALSE),
      names_to = "pos", values_to = "amino"
    )
    data <- dplyr::mutate(data,
      isTarget = ifelse(pos == "pos3", TRUE, FALSE),
      pos = (as.double(strx(pos, "[0-9]+")) - 6L) / 20,
      cand = paste0(Sequence_name, "_", PTM_type, "_", Position)
    )
    if (respective) {
      lst <- split(data, ~ Sequence_name)
    } else {
      lst <- list(All = data)
    }
    p.tops <- lapply(lst,
      function(data) {
        ggplot(data) +
          geom_text(aes(x = reorder(cand, PTM_score), y = pos, label = amino)) +
          geom_text(data = dplyr::filter(data, isTarget),
            aes(x = reorder(cand, PTM_score), y = min(data$pos) * 1.5, label = Position),
            color = "red", hjust = 0) +
          geom_point(data = dplyr::filter(data, isTarget),
            aes(x = reorder(cand, PTM_score), y = pos), color = "red", alpha = .3, size = 5) +
          geom_col(data = dplyr::filter(data, isTarget),
            aes(x = reorder(cand, PTM_score), y = PTM_score, fill = PTM_score)) +
          coord_flip() +
          ggforce::facet_grid_paginate(PTM_type ~ Sequence_name, scales = "free_y", space = "free") +
          scale_y_continuous(breaks = seq(0, 1, by = .2)) +
          geom_hline(yintercept = .5, linetype = 4) +
          labs(x = "PTM position", y = "PTM score") +
          guides(fill = "none") +
          theme_minimal() +
          theme(axis.text.y = element_blank())
      })
    p.tops <- .set_lab(p.tops, sig(x), names(p.tops), "PTM score")
    p.tops <- setLegend(p.tops, glue::glue("为 {names(p.tops)} 的修饰位点以及得分可视化图。"))
    lab(p.tops) <- paste0(sig(x), " ", "All PTM score")
    x@plots[[ 2 ]] <- namel(p.tops, p.PTMsNumber, p.hasPTMs)
    x@tables[[ 2 ]] <- namel(t.data, t.cutoff)
    feature(x) <- unique(t.cutoff$Sequence_name)
    x <- methodAdd(x, "以 Python 工具 `MusiteDeep` {cite_show('MusitedeepADWang2020')} 预测 {bind(x$type)} 修饰位点，设定 PTM 得分截断为 {cutoff}。")
    x <- snapAdd(x, "以 `MusiteDeep` 预测 {bind(x$type)} 修饰位点。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_musite"),
  function(x, n = 1)
  {
    step_message("Visualize the PDB structure")
    output_file <- paste0(x$output, "/ptm2Structure.json")
    cdRun(pg(x), " ", pg("musitePTM2S"),
      " -ptmInput ", x$seqs_file,
      " -ptmOutput ", x$output, "/res_results.txt",
      " -o ", x$output,
      " -maxPDB ", n
    )
    return(x)
  })

