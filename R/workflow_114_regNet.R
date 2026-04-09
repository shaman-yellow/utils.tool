# ==========================================================================
# workflow of regNet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_regNet <- setClass("job_regNet", 
  contains = c("job"),
  prototype = prototype(
    pg = "regNet",
    info = c(
      targetScanHuman = "https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi",
      mirdb = "https://mirdb.org/download.html",
      encori = "https://rnasysu.com/encori/tutorialAPI.php",
      miRNetR = "https://github.com/xia-lab/miRNetR/blob/master/vignettes/miRNetR.Rmd",
      tarbase = "https://dianalab.e-ce.uth.gr/tarbasev9/downloads",
      LncBase = "https://diana.e-ce.uth.gr/lncbasev3/home",
      npinter5 = "http://bigdata.ibp.ac.cn/npinter5/download/"
    ),
    cite = "",
    method = "",
    tag = "regNet",
    analysis = "分子调控网络分析"
    ))

setGeneric("asjob_regNet",
  function(x, ...) standardGeneric("asjob_regNet"))

setMethod("asjob_regNet", signature = c(x = "feature"),
  function(x){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_regNet(object = fea)
    return(x)
  })

setMethod("step0", signature = c(x = "job_regNet"),
  function(x){
    step_message("Prepare your data with function `job_regNet`.")
  })

setMethod("step1", signature = c(x = "job_regNet"),
  function(x, recode = NULL, tsh = TRUE, tbs = TRUE, 
    miRNetR = FALSE, mdb = FALSE, miRNetR_cache = NULL)
  {
    step_message("Get miRNA data.")
    targets <- object(x)
    if (!is.null(recode)) {
      targets <- unique(c(targets, names(recode)))
      recode <- as.list(recode)
    }
    fun_recode <- function(data, col, recode) {
      if (!is.null(recode)) {
        data[[ col ]] <- dplyr::recode(
          data[[ col ]], !!!recode, .default = data[[ col ]]
        )
      }
      data
    }
    x$all_miRNA <- list()
    if (tsh) {
      targetScanHuman <- ftibble(get_url_data(
          "Conserved_Family_Info.txt",
          # "Predicted_Targets_Info.default_predictions.txt",
          "https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Family_Info.txt.zip",
          # "https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip",
          "targetScanHuman"
          ))
      colnames(targetScanHuman) <- formal_name(colnames(targetScanHuman))
      targetScanHuman <- dplyr::filter(targetScanHuman, Gene_Symbol %in% targets)
      which_not_in_data(targetScanHuman, "Gene_Symbol", targets)
      targetScanHuman <- fun_recode(
        targetScanHuman, "Gene_Symbol", recode
      )
      targetScanHuman <- dplyr::rename(
        targetScanHuman, mirna_name = miR_Family, gene_name = Gene_Symbol
      )
      targetScanHuman <- reframe_split(
        targetScanHuman, "mirna_name", "/"
      )
      targetScanHuman <- dplyr::mutate(
        targetScanHuman, mirna_name = ifelse(
          grpl(mirna_name, "^miR"), 
          mirna_name, paste0("miR-", mirna_name)
        ),
        mirna_name = paste0("hsa-", mirna_name)
      )
      targetScanHuman <- dplyr::relocate(
        targetScanHuman, mirna_name, gene_name
      )
      targetScanHuman <- set_lab_legend(
        targetScanHuman,
        glue::glue("{x@sig} miRNA-RNA data from targetScanHuman"),
        glue::glue("从 targetScanHuman 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$targetScanHuman <- targetScanHuman
      x <- methodAdd(
        x, "TargetScanHuman (<https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi>)，"
      )
    }
    if (tbs) {
      tarbase <- ftibble(get_url_data(
          "Homo_sapiens_TarBase-v9.tsv.gz",
          "https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz",
          "tarbase", fun_decompress = NULL
          ))
      tarbase <- dplyr::filter(tarbase, gene_name %in% targets)
      which_not_in_data(tarbase, "gene_name", targets)
      tarbase <- fun_recode(tarbase, "gene_name", recode)
      tarbase <- dplyr::relocate(tarbase, mirna_name, gene_name)
      tarbase <- set_lab_legend(
        tarbase,
        glue::glue("{x@sig} miRNA-RNA data from TarBase"),
        glue::glue("从 Tarbase 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$tarbase <- tarbase
      x <- methodAdd(
        x, "TarBase (<https://dianalab.e-ce.uth.gr/tarbasev9/downloads>)，"
      )
    }
    if (FALSE && mdb) {
      mirdb <- ftibble(get_url_data(
          "miRDB_v6.0_prediction_result.txt.gz",
          "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
          "mirdb", fun_decompress = NULL
          ))
      x$all_miRNA$mirdb <- mirdb
    }
    if (miRNetR) {
      if (!is.null(miRNetR_cache)) {
        miRNetR <- ftibble(miRNetR_cache)
      } else {
        miRNetR <- callr::r(
          map_genes_in_miRNetR,
          list(genes = targets, wd = .prefix("miRNetR", "db"),
            extra = "--no-check-certificate"), show = TRUE
          )$mir.res
        miRNetR <- tibble::as_tibble(miRNetR)
      }
      which_not_in_data(miRNetR, "Target", targets)
      miRNetR <- fun_recode(miRNetR, "Target", recode)
      miRNetR <- dplyr::relocate(miRNetR, mirna_name = ID, gene_name = Target)
      miRNetR <- dplyr::mutate(
        miRNetR, mirna_name = gs(mirna_name, "hsa-mir", "hsa-miR")
      )
      miRNetR <- set_lab_legend(
        miRNetR,
        glue::glue("{x@sig} miRNA-RNA data from miRNet"),
        glue::glue("从 miRNet 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$miRNetR <- miRNetR
      x <- methodAdd(x, "miRnet (<https://www.mirnet.ca/>)，")
    }
    x <- methodAdd(x, "以上数据库用于检索以基因集 (mRNA) 为靶点 miRNA 数据。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_regNet"),
  function(x, use = "all"){
    step_message("miRNA intersection.")
    sets <- x$all_miRNA
    if (!identical(use, "all")) {
      sets <- sets[ use ]
    }
    ins <- .merge_list_by_cols(sets, by = c("mirna_name", "gene_name"))
    ins <- dplyr::distinct(ins, mirna_name, gene_name)
    message(glue::glue("Got data: {try_snap(ins, 'gene_name', 'mirna_name')}"))
    x$ins_mirna <- ins
    ins.snap <- try_snap(ins, "gene_name", "mirna_name")
    x <- snapAdd(x, "将数据库 {bind(names(sets))} 预测或记录的 mRNA 的上游 miRNA，二者 mRNA-miRNA 关系对取交集，共得到 {nrow(ins)} 对调控关系【{ins.snap}】。")
    x <- methodAdd(x, "取以上数据库的交集，确定高可信度的miRNA-mRNA调控配对。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_regNet"),
  function(x, enc = TRUE, npi = TRUE, mi = unique(x$ins_mirna$mirna_name)){
    step_message("Got lncRNA.")
    x$all_lncRNA <- list()
    if (enc) {
      encori <- get_encori_miRNA_lncRNA(mi)
      encori <- dplyr::relocate(encori, mirna_name = miRNAname, lncrna_name = geneName)
      encori <- set_lab_legend(
        encori,
        glue::glue("{x@sig} miRNA-lncRNA data from ENCORI"),
        glue::glue("从 ENCORI 获取的 miRNA-lncRNA 数据。")
      )
      x$all_lncRNA$encori <- encori
      x <- methodAdd(x, "ENCORI (<https://rna.sysu.edu.cn/encori>)，")
    }
    if (npi) {
      npinter <- ftibble(get_url_data(
          "miRNA_interaction.txt.gz",
          "http://bigdata.ibp.ac.cn/npinter5/download/file/miRNA_interaction.txt.gz",
          "npinter_miRNA_interaction", fun_decompress = NULL
          ))
      npinter <- dplyr::relocate(npinter, mirna_name = V5, lncrna_name = V2)
      if (!grpl(npinter[1, 1, drop = TRUE], "miR")) {
        stop('!grpl(npinter[1, 1, drop = TRUE], "miR")')
      }
      npinter <- dplyr::filter(npinter, mirna_name %in% mi)
      which_not_in_data(npinter, "mirna_name", mi)
      npinter <- set_lab_legend(
        npinter,
        glue::glue("{x@sig} miRNA-lncRNA data from NPinter"),
        glue::glue("从 NPinter 获取的 miRNA-lncRNA 数据。")
      )
      x$all_lncRNA$npinter <- npinter
      x <- methodAdd(x, "NPInter (<http://bigdata.ibp.ac.cn/npinter5>)，")
    }
    x <- methodAdd(x, "数据库用于获取与上述 miRNA 靶向结合的 lncRNA。")
    return(x)
  })

setMethod("step4", signature = c(x = "job_regNet"),
  function(x, use = "all"){
    step_message("lncRNA intersection")
    sets <- x$all_lncRNA
    if (!identical(use, "all")) {
      sets <- sets[ use ]
    }
    ins <- .merge_list_by_cols(sets, by = c("lncrna_name", "mirna_name"))
    ins <- dplyr::distinct(ins, lncrna_name, mirna_name)
    message(glue::glue("Got data: {try_snap(ins, 'mirna_name', 'lncrna_name')}"))
    ins.snap <- try_snap(ins, "mirna_name", "lncrna_name")
    x <- snapAdd(x, "将数据库 {bind(names(sets))} 预测或记录的 miRNA 的上游 lncRNA，二者 miRNA-lncRNA 关系对取交集，共得到 {nrow(ins)} 对调控关系【{ins.snap}】。")
    x <- methodAdd(x, "取两个数据库结果的交集筛选出高可信度的 lncRNA-miRNA 调控配对。")
    x$ins_lncrna <- ins
    return(x)
  })

setMethod("step5", signature = c(x = "job_regNet"),
  function(x, layout = "fr"){
    step_message("Network.")
    funName <- function(x) {
      setNames(x, c("from", "to"))
    }
    edges <- rbind(funName(x$ins_mirna), funName(x$ins_lncrna))
    nodes <- list(
      mRNA = x$ins_mirna$gene_name, miRNA = x$ins_mirna$mirna_name,
      lncRNA = x$ins_lncrna$lncrna_name
    )
    nodes <- as_df.lst(lapply(nodes, unique), "type", "name")[, 2:1]
    graph <- igraph::graph_from_data_frame(
      edges, directed = TRUE, vertices = nodes
    )
    layout <- ggraph::create_layout(graph, layout = layout)
    p.regNet <- ggraph(layout) + 
      geom_edge_fan(
        edge_width = .5, color = "black", show.legend = FALSE,
        end_cap = ggraph::circle(7, 'mm'),
        arrow = arrow(length = unit(1, 'mm'))) + 
      geom_node_point(aes(color = type, size = type, shape = type)) + 
      ggrepel::geom_label_repel(aes(x = x, y = y, label = name), size = 3) +
      guides(
        size = "none", shape = "none",
        color = guide_legend(override.aes = list(size = 4))
      ) +
      scale_size_manual(values = c(mRNA = 10, miRNA = 6, lncRNA = 6)) +
      scale_shape_manual(
        values = c(mRNA = 16, miRNA = 17, lncRNA = 18)
      ) +
      labs(color = "Type") +
      theme_void()
    p.regNet <- set_lab_legend(
      p.regNet,
      glue::glue("{x@sig} expression regulation networking"),
      glue::glue("lncRNA-miRNA-mRNA 表达网络分析|||图中的节点表示对应 RNA 类型，边代表相互作用。")
    )
    x <- snapAdd(x, "构建 lncRNA-miRNA-mRNA 表达网络{aref(p.regNet)}，如图所示，共 {nrow(nodes)} 个节点，{nrow(edges)} 个边。")
    x <- methodAdd(x, "以 R 包 `ggraph` ⟦pkgInfo('ggraph')⟧ 与 `ggplot2` ⟦pkgInfo('ggplot2')⟧ 对 lncRNA-miRNA-mRNA 转录后调控网络进行整合，可视化多层级分子调控网络的整体结构。")
    x <- plotsAdd(x, p.regNet)
    return(x)
  })

.merge_list_by_cols <- function(sets, by, keep = seq_along(by)) {
  if (length(sets) > 1) {
    ins <- sets[[1]][, keep]
    for (i in 2:length(sets)) {
      ins <- merge(ins, sets[[ i ]][, keep], by = by)
    }
  } else {
    ins <- sets
  }
  ins
}


get_encori_miRNA_lncRNA <- function(miRNA, dir_db = .prefix("encori", "db"))
{
  if (any(grpl(miRNA, "^miR-"))) {
    stop('any(grpl(miRNA, "^miR-")), the organisms such as hsa- should be added to prefix.')
  }
  # IDs: your query, col: the ID column, res: results table
  dir.create(dir_db, FALSE)
  db <- new_db(file.path(dir_db, "encori_miRNA_lncRNA.rdata"), ".id")
  db <- not(db, miRNA)
  query <- db@query
  if (length(query)) {
    res <- pbapply::pbsapply(query, .api_encori_miRNA_lncRNA, simplify = FALSE)
    res <- frbind(res, idcol = ".id")
    db <- upd(db, res)
  }
  res <- dplyr::filter(db@db, .id %in% !!miRNA)
}

.api_encori_miRNA_lncRNA <- function(miRNA) {
  if (length(miRNA) > 1) {
    stop('length(miRNA) > 1.')
  }
  url <- glue::glue("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA={miRNA}&interNum=1&expNum=1&cellType=all")
  data.table::fread(text = RCurl::getURL(url))
}

get_url_data <- function(expect_filename, url,
  name = formal_name(tools::file_path_sans_ext(basename(expect_filename))),
  dir = .prefix(name, "db"), fun_decompress = utils::unzip)
{
  expect_file <- file.path(dir, expect_filename)
  if (!file.exists(expect_file)) {
    dir.create(dir, FALSE)
    download_file <- file.path(
      dir, paste0(name, ".", tools::file_ext(url))
    )
    if (is.null(fun_decompress)) {
      download_file <- expect_file
    }
    utils::download.file(url, destfile = download_file)
    if (!is.null(fun_decompress)) {
      fun_decompress(download_file, exdir = dir)
    }
  }
  message(glue::glue("The `expect_file` exist: {file.exists(expect_file)}"))
  normalizePath(expect_file)
}

which_not_in_data <- function(data, col, items, prefix = NULL) {
  name <- rlang::expr_text(substitute(data))
  whichNot <- !items %in% data[[ col ]]
  if (any(whichNot)) {
    if (!is.null(prefix)) {
      message(glue::glue("{prefix} {name} not got: {bind(items[ whichNot ])}"))
    } else {
      message(glue::glue("{name} not got: {bind(items[ whichNot ])}"))
    }
  }
  return(whichNot)
}

map_genes_in_miRNetR <- function(genes, wd, method = "wget", extra = "--no-check-certificate")
{
  if (length(ls(envir = .GlobalEnv))) {
    stop("`miRNetR` is a terrible package with Poor code standards.
      Please run this function in a completely isolated subprocess!")
  }
  dir.create(wd, FALSE)
  setwd(wd)
  message(glue::glue("Set work directory in {normalizePath(wd)}.
      This is done because miRNetR will download the file to the current workspace."))
  fileWeb <- "https://www.xialab.ca/rest/sqlite/mir2gene.sqlite"
  filename <- basename(fileWeb)
  if (!file.exists(filename)) {
    download.file(fileWeb, filename, method = method, extra = extra)
  }
  miRNetR::Init.Data("mir", "gene")
  miRNetR::SetupIndListData(genes, "hsa", "gene", "symbol", "na", "na")
  nms.vec <<- c("gene")
  miRNetR::SetCurrentDataMulti()
  sqlite.path <<- "https://www.xialab.ca/rest/sqlite/"
  .on.public.web <<- FALSE
  miRNetR::QueryMultiList()
  dataSet
}

get_miRNetR.huibang <- function() {
  pak::pkg_install
}

