# ==========================================================================
# workflow of lzerd
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_lzerd <- setClass("job_lzerd", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = paste0("https://lzerd.kiharalab.org/upload/upload/", "\n",
      "http://cadd.zju.edu.cn/hawkdock/"),
    cite = "[@LzerdWebserverChrist2021; @HawkdockAWebWeng2019]",
    method = "`LZerD` and `HawkDock` webservers used for protein–protein docking"
    ))

job_lzerd <- function(symbols)
{
  .job_lzerd(object = symbols)
}

setMethod("step0", signature = c(x = "job_lzerd"),
  function(x){
    step_message("Prepare your data with function `job_lzerd`.")
  })

setMethod("step1", signature = c(x = "job_lzerd"),
  function(x){
    step_message("Download pdb files for protein docking.")
    mart <- new_biomart()
    x$anno <- filter_biomart(mart, c("hgnc_symbol", "pdb"), "hgnc_symbol",
      object(x), distinct = F)
    x$anno <- dplyr::filter(x$anno, pdb != "")
    data <- dplyr::distinct(x$anno, hgnc_symbol, .keep_all = T)
    x$layouts <- combn(data$pdb, 2, simplify = F)
    x$pdb_files <- get_pdb(unique(unlist(x$layouts)))
    x$pdb_from <- nl(data$hgnc_symbol, data$pdb, F)
    return(x)
  })

setMethod("map", signature = c(x = "job_lzerd"),
  function(x, ref, id, rev = F, tops = 1){
    ref <- match.arg(ref, c("lzerd", "hawkdock"))
    message("Create dir and naming as `id` as savepath.")
    mapped <- list()
    if (ref == "hawkdock") {
      .get_hawkdock_results(id)
      labels <- paste0(names(x$pdb_from), ": ", unname(x$pdb_from))
      if (rev) {
        labels <- rev(labels)
      }
      p.tops <- lapply(tops,
        function(top) {
          pdb_top <- paste0(id, "/top10/model.", top, ".pdb")
          save <- paste0(get_savedir("figs"), "/",
            make.names(labels[1]), "_with_", make.names(labels[2]),
            "_top", top, ".png"
          )
          if (!file.exists(save)) {
            res <- .file_fig(vis_pdb.hawkdock(pdb_top, labels[1], labels[2], save))
          } else {
            res <- .file_fig(save)
          }
          res <- .set_lab(res, sig(x), "HawkDock docking top", top)
          res
        })
      names(p.tops) <- paste0("top", tops)
      show_score <- function() {
        data <- head(ftibble(paste0(id, "/score.tsv")), n = 10)
        p <- ggplot(data) +
          geom_col(aes(x = reorder(V1, V2, decreasing = T), y = V2, fill = V2), width = .7) +
          labs(x = "Candidate", y = "Score") +
          guides(fill = "none") +
          coord_flip()
        wrap(p, 5, 4)
      }
      p.score <- show_score()
      p.score <- .set_lab(p.score, sig(x), "HawkDock ranking of all top 10 docking")
      mapped$h <- namel(p.tops, p.score)
    }
    x$mapped <- mapped
    return(x)
  })

.get_hawkdock_results <- function(id, savedir = id) {
  if (!file.exists(savedir)) {
    dir.create(savedir)
    url_score <- paste0("http://cadd.zju.edu.cn/hawkdock/statics/calculate/PDBfiles/", id, "/top10/HawkDock.dat")
    url_model <- paste0("http://cadd.zju.edu.cn/hawkdock/statics/calculate/PDBfiles/", id, "/top10/Top10-models.tar.gz")
    data <- ftibble(RCurl::getURL(url_score))
    write_tsv(data, sf <- paste0(savedir, "/score.tsv"))
    bin <- RCurl::getURLContent(url_model)
    writeBin(as.vector(bin), tf <- paste0(savedir, "/top10.tar.gz"))
    utils::untar(tf, exdir = savedir)
    namel(sf, tf)
  }
}

.show_hawkdock_results <- function(dir, file_score = "score.tsv", file_pdb = "top10/model.1.pdb")
{
  data <- head(ftibble(paste0(dir, "/", file_score)), n = 10)
  p.score <- ggplot(data) +
    geom_col(aes(x = reorder(V1, V2, decreasing = T), y = V2, fill = V2), width = .7) +
    labs(x = "Candidate", y = "Score") +
    guides(fill = "none") +
    coord_flip()
  p.score <- wrap(p.score, 5, 4)
  namel(p.score)
}

vis_pdb.hawkdock <- function(pdb, label.a, label.b,
  save = paste0(get_savedir("figs"), "/", make.names(label.a), "_with_", make.names(label.b), ".png"),
  command = "pymol")
{
  data <- ftibble(pdb, fill = T)
  coa <- dplyr::filter(data, V13 == "1.A")
  cob <- dplyr::filter(data, V13 == "1.B")
  fun_pos <- function(x) {
    c(median(x$V7), max(x$V8), max(x$V9) + diff(range(x$V9)) / 10)
  }
  pos.a <- fun_pos(coa)
  pos.b <- fun_pos(cob)
  expr <- .expr_pretty_pdb(label.a, pos.a, label.b, pos.b)
  vis_pdb(pdb, expr, save)
  return(save)
}

vis_pdb <- function(file, expr = NULL, save = NULL, command = "pymol") {
  if (!is.null(save)) {
    expr <- paste0(expr, " -g ", save)
  }
  cdRun(command, " ", file, " ", expr)
}

.expr_pretty_pdb <- function(label.a, pos.a, label.b, pos.b,
  a = "chain A", b = "chain B")
{
  pos <- function(x) paste0("[", paste0(x, collapse = ","), "]")
  paste0(
    " -d \"",
    "color pink, ", a, "; ",
    "color yellow, ", b, "; ",
    "pseudoatom pA, pos = ", pos(pos.a), "; ",
    "pseudoatom pB, pos = ", pos(pos.b), "; ",
    "set label_size, 25",
    "\" ",
    " -d \"", "label pA, '", label.a, "'", "\" ",
    "-d \"", "label pB, '", label.b, "'", "\" "
    )
}
