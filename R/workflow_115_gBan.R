# ==========================================================================
# workflow of gBan
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gBan <- setClass("job_gBan", 
  contains = c("job"),
  prototype = prototype(
    pg = "gBan",
    info = c("https://github.com/HamidHadipour/GraphBAN"),
    cite = "",
    method = "",
    tag = "gBan",
    analysis = "GraphBAN 药物预测"
    ))

job_gBan <- function()
{
  .job_gBan()
}

setGeneric("asjob_gBan",
  function(x, ...) standardGeneric("asjob_gBan"))

setMethod("asjob_gBan", signature = c(x = "feature"),
  function(x){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_gBan(object = fea)
    return(x)
  })

setMethod("step0", signature = c(x = "job_gBan"),
  function(x){
    step_message("Prepare your data with function `job_gBan`.")
  })

setMethod("step1", signature = c(x = "job_gBan"),
  function(x, dir_save = paste0("GraphBAN_", x@sig), 
    batman = TRUE, zinc = FALSE, file_batman_compounds_info = getOption("file_batman_compounds_info"))
  {
    step_message("Got amino acid sequence and drug smiles.")
    x$dir_save <- dir_save
    x$mart <- new_biomart("hsa")
    x$seqs <- get_seq.pro(object(x), x$mart)
    x$file_seqs <- write(
      x$seqs$fasta, name = "peptide", dir = dir_save, max = NULL
    )
    x <- methodAdd(x, "以 `biomaRt` ({packageVersion('biomaRt')}) 获取多肽序列 (先以 Symbol 获取 'ensembl_peptide_id' 以及 'transcript_is_canonical'，从而选择经典转录本的 'ensembl_peptide_id' 获取多肽序列) 。")
    x$smiles_compounds <- list()
    if (batman) {
      if (is.null(file_batman_compounds_info)) {
        stop('is.null(file_batman_compounds_info).')
      }
      if (!file.exists(file_batman_compounds_info)) {
        stop('!file.exists(file_batman_compounds_info).')
      }
      data_batman <- ftibble(file_batman_compounds_info)
      if (!any(colnames(data_batman) == "SMILES")) {
        stop('!any(colnames(data_batman) == "SMILES").')
      }
      x$smiles_compounds$data_batman <- dplyr::distinct(data_batman, smiles = SMILES)
      x <- methodAdd(x, "获取 BATMAN-TCM (<http://bionet.ncpsb.org.cn/batman-tcm>) 化合物 SMILES 结构式。")
    }
    if (zinc) {
      data_zinc <- ftibble(get_url_data(
        "250k_rndm_zinc_drugs_clean_3.csv",
        "https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv",
        "zinc", fun_decompress = NULL
      ))
      data_zinc <- dplyr::mutate(
        data_zinc, smiles = sub("\n$", "", smiles)
      )
      x$smiles_compounds$data_zinc <- data_zinc
      x <- methodAdd(x, "获取 ZINC-250k 小分子库 (<https://www.kaggle.com/datasets/basu369victor/zinc250k>) 的化合物 SMILES 结构式。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_gBan"),
  function(x, cl = 10, mem = 1000, w.cutoff = 1000, 
    rerun = FALSE, filter_by = c("rcdk", "rdkit"))
  {
    step_message("Filter drugs.")
    compounds <- unique(unlist(lapply(x$smiles_compounds, function(x) x$smiles)))
    filter_by <- match.arg(filter_by)
    if (filter_by == "rdkit") {
      x$info_compounds <- expect_local_data(
        "tmp", "compounds_weight_rdkit", inBatches_get_compounds_weight.rdkit,
        list(smiles_list = compounds), rerun = rerun
      )
      x <- methodAdd(
        x, "以 RDKit 过滤无法解析的分子结构，剔除分子量 &gt; {w.cutoff}，含重金属 (以原子序号大于钙为条件过滤) 的化合物，生成标准化 SMILES。"
      )
    } else if (filter_by == "rcdk") {
      x$info_compounds <- expect_local_data(
        "tmp", "compounds_weight", inBatches_get_compounds_weight.rcdk,
        list(smiles_list = compounds, mem = mem, cl = cl), ignore = "cl", rerun = rerun
      )
      x <- methodAdd(
        x, "以 R 包 `rcdk` ({packageVersion('rcdk')}) 过滤无法解析的分子结构，剔除分子量 &gt; {w.cutoff}，含重金属 (以原子序号大于钙为条件过滤) 的化合物。"
      )
    }
    input_compounds <- dplyr::filter(
      x$info_compounds, MolecularWeight < w.cutoff, !HasHeavyMetal
    )
    input_compounds <- dplyr::distinct(input_compounds, smiles)
    x$input_compounds <- dplyr::mutate(input_compounds, id = seq_len(nrow(input_compounds)), .before = 1)
    combn <- expand.grid(x$input_compounds$id, x$seqs$data$hgnc_symbol)
    combn <- dplyr::rename(combn, id = 1, hgnc_symbol = 2)
    combn <- Reduce(merge, list(combn, x$input_compounds, x$seqs$data))
    combn <- dplyr::mutate(combn, Y = 1)
    combn <- dplyr::relocate(
      combn, hgnc_symbol, id, SMILES = smiles, Protein = peptide, Y
    )
    x$file_combn <- file.path(
      x$dir_save, "graphBan_input.csv"
    )
    write.csv(combn, x$file_combn, row.names = FALSE)
    x$combn <- combn
    x <- snapAdd(x, "从数据库获取到的化合物共 {length(compounds)} 个。经过滤后得到 {nrow(input_compounds)} 个唯一化合物。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_gBan"),
  function(x){
    step_message("Do nothing")
    return(x)
  })

setMethod("step4", signature = c(x = "job_gBan"),
  function(x, pattern = "graphBan_res_", cutoff = .95, 
    reRead = FALSE, method_keep = c("respective", "all"))
  {
    step_message("Collate results.")
    method_keep <- match.arg(method_keep)
    files_res <- list.files(x$dir_save, pattern, full.names = TRUE)
    if (!length(files_res)) {
      stop('length(files_res), no results file found.')
    }
    res <- pbapply::pblapply(files_res, 
      function(file) {
        res <- expect_local_data(
          "tmp", "gbanResRead", ftibble, list(files = file, select = "pred"), rerun = reRead
        )
        if (nrow(res) != nrow(x$combn)) {
          stop('nrow(res) != nrow(x$combn).')
        }
        whichKeep <- which(res$pred > cutoff)
        res <- dplyr::bind_cols(
          res[whichKeep, ], x$combn[whichKeep, ]
        )
        return(res)
      })
    names(res) <- s(basename(files_res), "\\..*$", "")
    res <- dplyr::bind_rows(res, .id = "model")
    x$res_graphBan <- res
    split_by_genes <- split(res, res$hgnc_symbol)
    s.eachModel <- vapply(
      split_by_genes, try_snap, character(1), main = "model", sub = "id"
    )
    s.eachModel <- paste0(
      paste0(
        s(names(s.eachModel), pattern, ""), " 得到候选药物 ", s.eachModel
      ), collapse = "。"
    )
    s.eachTarget <- bind(try_snap(split_by_genes, "id"))
    x <- methodAdd(x, "基于 GraphBAN（Graph-Based Attention Network）的CPI（化合物-蛋白质相互作用）预测框架，结合 BindingDB、BioSNAP 和 KIBA 等公开药物–靶点互作数据集构建基础图网络，用于模型训练与验证，输出每对“化合物-蛋白质”组合的相互作用概率。")
    if (method_keep == "all") {
      ins <- unique(ins(lst = lapply(split_by_genes, function(x) x$SMILES)))
      s.ins <- glue::glue("共同作用于各靶点的药物有 {length(ins)} 个。")
      x$smiles_keep <- ins
      x$split_by_genes <- lapply(
        split_by_genes, function(x) x[x$SMILES %in% ins, ]
      )
    } else if (method_keep == "respective") {
      x$smiles_keep <- unique(res$SMILES)
      s.ins <- ""
      x$split_by_genes <- split_by_genes
    }
    x$smiles_from_gban <- x$smiles_keep
    x$file_smiles_for_admet <- file.path(x$dir_save, "smiles_for_admet.txt")
    writeLines(x$smiles_keep, x$file_smiles_for_admet)
    x <- snapAdd(
      x, "经 GraphBAN (BindingDB, BioSNAP, KIBA 模型) 预测，筛选药物与靶点互作概率大于 {cutoff} 的化合物。{s.eachModel} (三大模型的合集)。各靶点的唯一化合物统计为 {s.eachTarget}。{s.ins}"
    )
    return(x)
  })

setMethod("step5", signature = c(x = "job_gBan"),
  function(x, file_admet = NULL, cutoff = .7)
  {
    step_message("ADMET ...")
    smiles <- x$smiles_from_gban
    if (is.null(file_admet)) {
      stop('is.null(file_admet).')
    }
    x$admet <- admet <- ftibble(file_admet)
    if (nrow(admet) != length(smiles)) {
      stop('nrow(admet) != length(smiles).')
    }
    x$admet_filter <- admet <- dplyr::filter(
      admet,
      caco2 > -5.15, hia < cutoff, PPB < 90,
      DILI < cutoff, Ames < cutoff,
      dplyr::if_all(tidyselect::matches("CYP.*inh"), ~ .x < cutoff)
    )
    x$smiles_keep <- x$smiles_from_admet <- admet$raw_smiles
    x$file_smiles_for_swiss <- file.path(x$dir_save, "smiles_for_swiss.txt")
    writeLines(x$smiles_keep, x$file_smiles_for_swiss)
    message(glue::glue("After filtered by admet: {nrow(admet)}"))
    x <- methodAdd(
      x, "以 ADMETlab 3.0 (https://admetlab3.scbdd.com/server/screening) 对交集化合物开展系统性 ADMET 分析，评估其药代动力学特征 (caco2 &gt; -5.15, hia &lt; {cutoff}, PPB &lt; 90, DILI &lt; {cutoff}, Ames &lt; {cutoff}, All CYP inhibitor &lt; {cutoff}) (参考 <https://admetlab3.scbdd.com/explanation/#/>) "
    )
    x <- snapAdd(x, "以 ADMETlab 评估其药物特性后，获得 {nrow(admet)} 个化合物。")
    return(x)
  })

setMethod("step6", signature = c(x = "job_gBan"),
  function(x, file_swiss = NULL, n_pass = 4, method = "linpinski"){
    step_message("swissADME.")
    if (is.null(file_swiss)) {
      stop('is.null(file_swiss)')
    }
    swissAdme <- ftibble(file_swiss)
    if (nrow(swissAdme) != length(x$smiles_from_admet)) {
      stop('nrow(swissAdme) != x$smiles_from_admet.')
    }
    swissAdme <- dplyr::mutate(
      swissAdme, raw_smiles = x$smiles_from_admet, .after = 1
    )
    colnames(swissAdme) <- formal_name(colnames(swissAdme))
    x$swissAdme <- swissAdme
    x <- methodAdd(x, "利用 SwissADME 平台 (<http://www.swissadme.ch/>) 评估核心候选化合物的成药性。")
    if (method == "linpinski") {
      prins <- list(
        quote(MW < 500),
        quote(X_H_bond_acceptors < 10),
        quote(X_H_bond_donors <= 5),
        quote(Consensus_Log_P <= 4.15),
        quote(TPSA < 140)
      )
      passes <- lapply(prins, 
        function(prin) {
          dplyr::filter(swissAdme, !!prin)$raw_smiles
        })
      res <- table(unlist(passes))
      x$smiles_keep <- x$smiles_from_swiss <- names(res)[ res > n_pass ]
      x <- methodAdd(x, "计算分子量 (MW)、氢键受体 (HBA) 数量、氢键供体 (HBD) 数量、LogP 值、拓扑极性表面积 (TPSA)，依据 Lipinski 五法则筛选 (标准：MW &lt; 500、HBA &lt; 10、HBD ≤ 5、LogP ≤ 4.15、TPSA &lt; 140 A²)。违反条目 ≤ {5-n_pass} 条的化合物判定为具有良好口服生物利用度及成药潜力。")
    }
    x <- snapAdd(x, "以 SwissADME 评估其药物特性后，余下 {length(x$smiles_keep)} 个化合物。")
    return(x)
  })

setMethod("step7", signature = c(x = "job_gBan"),
  function(x){
    expect_package("PubChemR", "3.0.0")
    smiles <- x$smiles_keep
    cids <- expect_local_data(
      "tmp", "pubchemr_cids", PubChemR::get_cids, list(
        identifier = smiles, namespace = "smiles"
      )
    )
    x$cids <- cids <- e(PubChemR::CIDs(cids))
    numNotGot <- sum(cids$CID == 0)
    message(glue::glue("Not got: {numNotGot}"))
    cids <- dplyr::filter(cids, CID != 0)
    x <- snapAdd(x, "以 R 包 `PubChemR` ({packageVersion('PubChemR')}) 从 PubChem 数据库 (<https://pubchem.ncbi.nlm.nih.gov/>) 检索化合物信息。将 PubChem 有记录条目的化合物视为候选药物 (n = {nrow(cids)})，并获取其对应化合物名。")
    x$synos <- try_get_syn(cids$CID)
    feature(x) <- as_feature(
      x$synos$Synonym, "候选药物", nature = "compounds"
    )
    return(x)
  })

# <https://admetlab3.scbdd.com/explanation/#/>
# 
# Ames, f30, hia (intestine), BBB, Drug-induced liver injury (DILI)
# 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0: poor (red)
# Caco-2
# > -5.15: excellent (green); otherwise: poor (red)
# PPB
# ≤ 90%: excellent (green); otherwise: poor (red)
# "CYP1A2-inh"  "CYP2C19-inh" "CYP2C9-inh"  "CYP2D6-inh"  "CYP3A4-inh"  "CYP2B6-inh"  "CYP2C8-inh"
# Category 0: Non-substrate / Non-inhibitor; Category 1: substrate / inhibitor. 


setMethod("asjob_vina", signature = c(x = "job_gBan"),
  function(x){
    # resolve_feature_snapAdd_onExit("x", feature(x))
    cpds <- map(x$synos, "CID", x$cids, "CID", "SMILES", col = "SMILES")
    layout <- dplyr::filter(x$res_graphBan, SMILES %in% cpds$SMILES)
    layout <- merge(layout, cpds, by = "SMILES", all.x = TRUE)
    layout <- dplyr::select(layout, Synonym, hgnc_symbol, CID)
    layout <- dplyr::mutate(layout, hgnc_symbol = as.character(hgnc_symbol))
    layout <- dplyr::distinct(layout)
    x <- job_vina(.layout = layout)
    x <- snapAdd(x, "对候选药物与对应靶点进行 AutoDock vina 分子对接。")
    return(x)
  })

run_remote_graphBan.huibang <- function(x,
  remote = "graphBan")
{
  remote_dir <- glue::glue(
    "~/{s(guess_project(), '^[0-9]+_', '')}"
  )
  files <- paste(
    x$file_combn, file.path(.expath, "job_templ", "graphBan.sh")
  )
  cmd_prepare <- glue::glue("mkdir {remote_dir}")
  cmd_send <- glue::glue("scp {files} {remote}:{remote_dir}")
  cmd_run <- glue::glue("ssh {remote} 'cd {remote_dir} && nohup sh graphBan.sh > task.log 2>&1 &'")
  writeLines(cmd_prepare)
  writeLines(cmd_send)
  writeLines(cmd_run)
}

inBatches_get_compounds_weight.rdkit <- function(smiles_list, python = getOption("rdkit_python"))
{
  stop(glue::glue("..."))
  # if (!is.null(python)) {
  #   e(base::Sys.setenv(RETICULATE_PYTHON = python))
  #   e(reticulate::use_python(python))
  #   e(reticulate::py_config())
  # }
  # rdkit::import("rdkit.Chem")
  # rdkit::import("rdkit.Chem.Descriptors")
  # rdkit::import("rdkit.Chem.rdMolDescriptors")
}

inBatches_get_compounds_weight.rcdk <- function(smiles_list, 
  cl = NULL, ..., mem = 1000, dir_db = .prefix("smiles_compounds_weight", "db"))
{
  # IDs: your query, col: the ID column, res: results table
  # smiles_list <- head(smiles_list, n = 3000)
  dir.create(dir_db, FALSE)
  db <- new_db(file.path(dir_db, "compounds_weight.rdata"), "smiles")
  db <- not(db, smiles_list)
  query <- db@query
  if (length(query)) {
    groups <- grouping_vec2list(query, mem, TRUE)
    message("Total group: ", length(groups))
    cli::cli_alert_info("rcdk::parse.smiles")
    res <- pbapply::pblapply(groups, cl = cl,
      function(smiles) {
        res <- try(silent = TRUE, callr::r(
          get_compounds_weight, args = c(list(smiles_list = smiles), list(...)),
          show = TRUE
          ))
        message(glue::glue("Group {attr(smiles, 'name')} finished."))
        if (!inherits(res, "try-error")) {
          res
        } else NULL
      })
    res <- frbind(res)
    db <- upd(db, res)
  }
  res <- dplyr::filter(db@db, smiles %in% !!smiles_list)
  res
}

get_compounds_weight <- function(smiles_list, HEAVY_METAL_THRESHOLD = 20, 
  molecules = NULL)
{
  # smiles_list <- head(smiles_list, n = 1000)
  if (is.null(molecules)) {
    molecules <- rcdk::parse.smiles(smiles_list)
  }
  results <- pbapply::pblapply(smiles_list,
    function(smile) {
      molecule <- molecules[[ smile ]]
      if (is.null(molecule)) {
        warning(paste("Unable to parse SMILES: ", smile))
        return()
      }
      rcdk::convert.implicit.to.explicit(molecule)
      formula <- rcdk::get.mol2formula(molecule, charge = 0)
      mol_weight <- formula@mass
      has_heavy_metal <- FALSE
      atoms <- rcdk::get.atoms(molecule)
      for (atom in atoms) {
        atomic_num <- rcdk::get.atomic.number(atom)
        if (atomic_num > HEAVY_METAL_THRESHOLD) {
          has_heavy_metal <- TRUE
          break
        }
      }
      data.frame(
        smiles = smile,
        MolecularWeight = mol_weight,
        HasHeavyMetal = has_heavy_metal,
        stringsAsFactors = FALSE
      )
    })
  do.call(dplyr::bind_rows, results)
  dplyr::as_tibble(data.table::rbindlist(results))
}

