# ==========================================================================
# workflow of gmix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gmix <- setClass("job_gmix", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/gmix/wiki"),
    cite = "[@TheDisgenetKnPinero2019; @TheGenecardsSStelze2016; @PharmgkbAWorBarbar2018]",
    method = "Databses of `DisGeNet`, `GeneCards`, `PharmGKB` used for collating disease related targets"
    ))

job_gmix <- function(disease, fuzzy = NULL,
  db_PharmGKB = get_pharmgkb_data("../database/PharmGKB_relationships.tsv.gz"),
  db_DisGeNet = get_disgenet_data("../database/DisGeNet_disease_associations.tsv.gz"))
{
  x <- .job_gmix(object = disease)
  x$fuzzy <- fuzzy
  x$db_pharm <- db_PharmGKB
  x$db_dis <- db_DisGeNet
  return(x)
}

setMethod("step0", signature = c(x = "job_gmix"),
  function(x){
    step_message("Prepare your data with function `job_gmix`.")
  })

setMethod("step1", signature = c(x = "job_gmix"),
  function(x, get_pharm = T, get_dis = T, get_genecard = T){
    step_message("Search databases.")
    fun <- function(name, col) {
      grpf(unique(x[[ name ]][[ col ]]), x$fuzzy, T)
    }
    fun_show <- function(name, max = 10) {
      if (length(x[[ name ]]) > max) {
        lst <- c(head(x[[name]], n = max), "...")
      } else {
        lst <- x[[name]]
      }
      paste0(1:length(lst), ". ", lst, "\n\t")
    }
    if (get_pharm) {
      x$lst_pharm <- fun("db_pharm", "Entity1_name")
      message("get_pharm:\n\t", fun_show("lst_pharm"))
    }
    if (get_dis) {
      x$lst_dis <- fun("db_dis", "diseaseName")
      message("get_dis:\n\t", fun_show("lst_dis"))
    }
    message("[Stats] Got items:\n\t",
      "get_pharm: ", length(x$lst_pharm), "\n\t",
      "get_dis: ", length(x$lst_dis)
    )
    x$get_genecard <- get_genecard
    return(x)
  })

setMethod("step2", signature = c(x = "job_gmix"),
  function(x, use.pharm = NULL, use.dis = NULL, use.score = NULL)
  {
    step_message("Get data.")
    tables <- list()
    lst.genes <- list()
    if (!is.null(use.pharm)) {
      cli::cli_alert_info("Get PharmGKB")
      t.pharm <- dplyr::filter(x$db_pharm, Entity1_name %in% x$lst_pharm[[ use.pharm ]])
      t.pharm <- .set_lab(t.pharm, sig(x), "PharmGKB used data")
      tables <- c(tables, namel(t.pharm))
      lst.genes[[ paste0("PharmGKB: ", x$lst_pharm[[ use.pharm ]]) ]] <- t.pharm$Entity2_name
    }
    if (!is.null(use.dis)) {
      cli::cli_alert_info("Get DisGeNet")
      data <- dplyr::filter(x$db_dis, diseaseName %in% x$lst_dis[[ use.dis ]])
      if (nrow(data) > 1) {
        stop("In DisGeNet, nrow(data) > 1")
      }
      t.dis <- .get_source.DisGeNet(data$diseaseId, data$NofGenes)
      t.dis <- .set_lab(t.dis, sig(x), "DisGeNet used data")
      tables <- c(tables, namel(t.dis))
      lst.genes[[ paste0("DisGeNet: ", x$lst_dis[[ use.dis ]]) ]] <- t.dis$Gene.Symbol
    }
    if (x$get_genecard) {
      cli::cli_alert_info("Get GeneCards")
      t.genecard <- get_from_genecards(object(x), use.score)
      t.genecard <- .set_lab(t.genecard, sig(x), "GeneCards used data")
      tables <- c(tables, namel(t.genecard))
      lst.genes[[ paste0("GeneCards: ", object(x)) ]] <- t.genecard$Symbol
    }
    p.cols <- new_col(lst = lst.genes)
    p.cols <- .set_lab(p.cols, sig(x), "Overall targets number of datasets")
    x@tables[[ 2 ]] <- tables
    x@plots[[ 2 ]] <- namel(p.cols)
    x$lst.genes <- lst.genes
    return(x)
  })

setMethod("map", signature = c(x = "job_gmix", ref = "job_gmix"),
  function(x, ref, extra = NULL)
  {
    fun <- function(x) rm.no(unlist(x@params$lst.genes, use.names = F))
    gene.x <- fun(x)
    gene.y <- fun(ref)
    lst <- nl(c(sig(x), sig(ref)), list(gene.x, gene.y))
    if (is.null(extra)) {
      p <- new_venn(lst = lst)
    } else {
      if (!is(extra, "list") || is.null(names(extra))) {
        stop('`extra` Show be a "list" with names')
      }
      lst <- c(lst, extra)
      p <- new_upset(lst = lst)
    }
    namel(raw = lst, p)
  })

.get_source.DisGeNet <- function(id, length) {
  url <- paste0("https://www.disgenet.org/browser/0/1/0/", id, "/0/", length, "/source__ALL/_b./#")
  html <- RCurl::getURL(url)
  data <- as_tibble(get_table.html(html)[[1]])
  fun_format <- function(x) {
    x <- dplyr::select(x, type = 1, value = 2)
    x <- dplyr::mutate(x, type = gs(type, "^(.+?):.*", "\\1"),
      type = ifelse(grpl(type, "^[A-Z]"), type, "UniProt"),
      type = make.names(type))
    rows <- nrow(x) / 3
    x <- split_lapply_rbind(x, rep(1:rows, each = 3),
      function(x) {
        tidyr::spread(x, type, value)
      })
    x
  }
  fun_format(data)
}

get_pharmgkb_data <- function(file_relation = "../database/PharmGKB_relationships.tsv.gz") {
  data <- ftibble(file_relation)
  data <- dplyr::filter(data, (Entity1_type == "Gene" & Entity2_type == "Disease") |
    (Entity2_type == "Gene" & Entity1_type == "Disease"))
  fun_format <- function(x) {
    n <- 0L
    split_lapply_rbind(x, ~ Entity1_type,
      function(x) {
        n <<- n + 1L
        if (n == 2L) {
          colnames(x) <- colnames(x)[c(4:6, 1:3, 7:ncol(x))]
        }
        return(x)
      })
  }
  fun_format(data)
}

get_disgenet_data <- function(file_disease = "../database/DisGeNet_disease_associations.tsv.gz") {
  data <- ftibble(file_disease)
  data
}
