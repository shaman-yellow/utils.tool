# ==========================================================================
# workflow of gutmd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gutmd <- setClass("job_gutmd", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: http://bio-annotation.cn/gutMDisorder"),
    cite = "[@GutmdisorderACheng2019]",
    method = "Database `gutMDisorder` used for finding associations between gut microbiota and metabolites",
    tag = "link:meta+micro",
    analysis = "GutMDisorder 获取微生物和代谢物关联数据"
    ))

job_gutmd <- function(db = .prefix("Gut Microbe and Metabolite-human.txt", "db"))
{
  if (!file.exists(db)) {
    data <- e(RCurl::getURL("http://bio-annotation.cn/gutmgene/public/res/Gut%20Microbe%20and%20Metabolite-human.txt"))
    write_tsv(data, db)
    db_data <- data.table::fread(text = data)
  } else {
    db_data <- ftibble(db)
  }
  db_data <- dplyr::rename_all(db_data, make.names)
  db_data <- .set_lab(db_data, "", "gutMDisorder database")
  .job_gutmd(object = db_data)
}

setMethod("step0", signature = c(x = "job_gutmd"),
  function(x){
    step_message("Prepare your data with function `job_gutmd`.")
  })

setMethod("step1", signature = c(x = "job_gutmd"),
  function(x, micro_patterns = NULL, metab_patterns = NULL, label = NULL, label.auto = T)
  {
    step_message("Match microbiota or metabolite in database")
    tables <- list()
    if (!is.null(micro_patterns)) {
      res <- list()
      res$matches <- matchThats(object(x)[[1]], micro_patterns)
      if (length(res$matches)) {
        db_matched <- dplyr::filter(object(x),
          Gut.Microbiota %in% unlist(!!res$matches, use.names = F)
        )
        dic <- as_df.lst(res$matches)
        dic <- nl(dic$name, dic$type)
        db_matched <- dplyr::mutate(db_matched, Query = dplyr::recode(Gut.Microbiota, !!!dic))
        db_matched <- dplyr::relocate(db_matched, Query)
        db_matched <- .set_lab(db_matched, sig(x), "Matched data in gutMDisorder")
        res$db_matched <- db_matched
        data <- dplyr::select(db_matched, Query, Gut.Microbiota, Substrate, Metabolite)
        p.allu <- new_allu(data, axes = 1:4, label.auto = label.auto, shiny = T)
        res$p.allu <- .set_lab(wrap(p.allu), sig(x), "alluvium plot of Matched data in gutMDisorder")
        res$data <- dplyr::distinct(data)
        x$res_micro <- res
      } else {
        message("No matched Metabolites for the `micro_patterns`")
      }
    }
    if (!is.null(metab_patterns)) {
      res <- list()
      alls <- rm.no(c(object(x)$Substrate, object(x)$Metabolite))
      res$matches <- matchThats(alls, metab_patterns)
      names <- unlist(res$matches, use.names = F)
      res$db_matched <- dplyr::filter(object(x),
        Substrate %in% !!names | Metabolite %in% !!names
      )
      x$res_metab <- res
      metab <- dplyr::select(res$db_matched, Metabolite, Substrate, Gut.Microbiota, Classification)
      metab <- .set_lab(metab, label, "gutMDisorder", "Matched metabolites and their related microbiota")
      tables$metab <- namel(metab)
    }
    x@tables <- tables
    return(x)
  })

setMethod("map", signature = c(x = "job_gutmd", ref = "job_publish"),
  function(x, ref)
  {
    if (ref@cite == "[@ProteinMetabolBenson2023]") {
      message("Gut.Microbiota (previous filtered) => Metabolite (or Substrate) => Protein")
      ## Gut.Microbiota related data
      data <- x@params$res_micro$db_matched
      data <- dplyr::select(data, Gut.Microbiota,
        Substrate, Substrate.PubChem.CID,
        Metabolite, Metabolite.PubChem.CID
      )
      ref.cids <- rm.no(object(ref)$pubchem_cid)
      data <- dplyr::filter(data,
        Substrate.PubChem.CID %in% !!ref.cids | Metabolite.PubChem.CID %in% !!ref.cids
      )
      data1 <- dplyr::mutate(data, .id = Substrate.PubChem.CID, .id_from = "Substrate")
      data2 <- dplyr::mutate(data, .id = Metabolite.PubChem.CID, .id_from = "Metabolite")
      main <- dplyr::bind_rows(data1, data2)
      main <- dplyr::filter(main, !is.na(.id))
      ##############
      ##############
      ## ref data to Genes
      cids <- rm.no(c(data$Substrate.PubChem.CID, data$Metabolite.PubChem.CID))
      ref.data <- dplyr::filter(object(ref), pubchem_cid %in% !!cids, !is.na(pubchem_cid))
      ##############
      ##############
      matched <- tbmerge(main, ref.data, by.x = ".id", by.y = "pubchem_cid", allow.cartesian = T)
      matched <- dplyr::relocate(matched, .id, .id_from, Substrate, Metabolite, Gut.Microbiota)
      matched <- dplyr::select(matched, -dplyr::ends_with("PubChem.CID"))
      matched <- dplyr::mutate(matched, Target_Gene = gs(Target_Gene, "\\s*,\\s*|\\s+", " "))
      fun_format <- function(data) {
        data1 <- dplyr::filter(data, grpl(Target_Gene, " "))
        data1 <- split_lapply_rbind(data1, 1:nrow(data1),
          function(x) {
            genes <- x$Target_Gene
            x <- dplyr::select(x, -Target_Gene)
            data.frame(x, Target_Gene = strsplit(genes, " ")[[1]])
          })
        data <- dplyr::filter(data, !grpl(Target_Gene, " "))
        dplyr::bind_rows(data, data1)
      }
      matched <- fun_format(matched)
      matched <- .set_lab(matched, sig(x), "Discover relationship between Microbiota with Host genes by matching metabolites")
      x$ProteinMetabolBenson2023 <- namel(matched)
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_limma", ref = "job_gutmd"),
  function(x, ref, tops, use = "hgnc_symbol", plot_max = 1000)
  {
    if (is.null(data <- ref$ProteinMetabolBenson2023$matched)) {
      stop("is.null(ref$ProteinMetabolBenson2023$matched) == T")
    }
    sigs <- rm.no(tops[[ use ]])
    lst <- list(Microbiota_associated = data$Target_Gene, DEGs = sigs)
    p.venn <- new_venn(lst = lst)
    p.venn <- .set_lab(p.venn, sig(x), "Intersection of Microbiota associated Genes with DEGs")
    data <- dplyr::filter(data, Target_Gene %in% p.venn$ins)
    data <- dplyr::arrange(data, dplyr::desc(abs(META_Rho)))
    data <- .set_lab(data, sig(x), "Microbiota associated Genes filtered by DEGs")
    data.network <- dplyr::select(head(data, n = plot_max), Gut.Microbiota, Metabolite, Target_Gene)
    data.network <- .set_lab(data.network, sig(x), paste0("Top ", plot_max, " relationship data"))
    p.network <- plot_network.pharm(data.network,
      ax1 = "Microbiota", ax2 = "Metabolite", ax3 = "DEGs", less.label = F
    )
    p.network <- .set_lab(p.network, sig(x), paste0("Top ", plot_max, " relationship network"))
    namel(p.venn, data, p.network, data.network)
  })

