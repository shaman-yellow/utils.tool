# ==========================================================================
# workflow of miranda
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_miranda <- setClass("job_miranda", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@Enrigh2003; @MirbaseFromMKozoma2019]",
    method = paste0("The tool of `miRanda` used for predicting mRNA targets for microRNAs (miRNA) ",
      "and `miRBase` used for getting sequence of miRNA")
    ))

job_miranda <- function(mirna, rna)
{
  message("The `rna` should be symbol (hgnc_symbol)\n",
    "The `mirna` can be character (pattern) for match miRNA.")
  mi <- .prepare_mirna(mirna)
  x <- .job_miranda(object = list(mirna = mi$mirna, rna = gname(rna)))
  x$mi <- mi
  x
}

setMethod("step0", signature = c(x = "job_miranda"),
  function(x){
    step_message("Prepare your data with function `job_miranda`."
    )
  })

setMethod("step1", signature = c(x = "job_miranda"),
  function(x, wd = timeName("miranda"), rna_type = "coding"){
    step_message("Prepare computional files.")
    if (is.null(x$mart)) {
      mart <- new_biomart()
      x$mart <- mart
    } else {
      mart <- x$mart
    }
    dir.create(wd, F)
    x$wd <- wd
    # miRNA seq
    x$mirna_seq <- x$mi$mirdb[ object(x)$mirna ]
    x$mirna_file <- paste0(wd, "/", make.names(object(x)$mirna), ".fa")
    e(Biostrings::writeXStringSet(x$mirna_seq, x$mirna_file))
    # RNA seq
    x$rna_seq <- get_seq.rna(object(x)$rna, mart, to = rna_type)
    x$rna_file <- write(x$rna_seq$fasta, paste0(wd, "/rna"))[[1]]
    return(x)
  })

setMethod("step2", signature = c(x = "job_miranda"),
  function(x, command = "conda run -n base miranda"){
    step_message("Run miranda.")
    x$outfile <- paste0(x$wd, "/", "results.txt")
    cdRun(command,
      " ", x$mirna_file,
      " ", x$rna_file,
      " -out ", x$outfile
    )
    x$alls <- .read_miranda_results(x$outfile)
    names(x$alls$cans) <- object(x)$rna
    t.overs <- lapply(x$alls$cans,
      function(x) {
        x <- lapply(x, function(x) data.frame(score = x$score, energy = x$energy))
        data.table::rbindlist(x)
      })
    names(t.overs) <- object(x)$rna
    t.overs <- as_tibble(data.table::rbindlist(t.overs, idcol = T))
    t.overs <- dplyr::relocate(t.overs, ref = .id)
    p.tops <- lapply(x$alls$cans,
      function(lst) {
        p <- .plot_unit.miranda(lst[[1]])
      })
    p.tops <- .set_lab(p.tops, sig(x), names(p.tops), "Top match")
    p.overs <- .plot_score.miranda(t.overs)
    p.overs <- .set_lab(p.overs, sig(x), "overview of all candidates")
    x@tables[[2]] <- namel(t.overs)
    x@plots[[2]] <- namel(p.tops, p.overs)
    return(x)
  })

.prepare_mirna <- function(mirna) {
  mirdb <- .get_mature_miRNA.mirbase()
  if (grpl(mirna, "^miR")) {
    message("Use `hsa` as default species.")
    mirna <- paste0("hsa-", mirna)
  }
  foundThat <- grpf(names(mirdb), mirna, ignore.case = T)
  message("All found miRNA:\n\n", stringr::str_trunc(paste0(foundThat, collapse = "\n"), 100), "\n")
  if (length(foundThat) > 1) {
    message("Use the first.")
  }
  mirna <- foundThat[1]
  namel(mirna, foundThat, mirdb)
}

.read_miranda_results <- function(file) {
  lines <- readLines(file)
  sep.h <- grp(lines, "^Read Sequence")[1]
  head <- lines[1:sep.h]
  body <- lines[(sep.h + 1):length(lines)]
  lst <- sep_list(body, "^Performing Scan", T)
  cans <- lapply(lst,
    function(x) {
      it <- sep_list(x, "^>")
      it <- lapply(it,
        function(x) {
          if (any(grpl(x, "Query"))) x
        })
      it <- lst_clear0(it)
      lapply(it, .read_unit.miranda)
    })
  cans <- lst_clear0(cans)
  namel(head, cans)
}

.plot_score.miranda <- function(data) {
  data <- split_lapply_rbind(data, ~ref,
    function(x) {
      x <- dplyr::mutate(x, candidate = paste0("C ", 1:nrow(x)))
    })
  p <- ggplot(data) +
    geom_col(aes(x = reorder(candidate, score), y = score, fill = energy), width = .7) +
    facet_wrap(~ ref) +
    coord_flip(ylim = zoRange(data$score, 1.2)) +
    labs(x = "Score", y = "Candidates", fill = "Energy") +
    scale_fill_gradientn(colors = rev(ggsci::pal_npg()(3))) +
    theme_minimal()
  if (length(unique(data$ref)) == 1) {
    return(wrap(p, 8, 1 + nrow(data) * .5))
  }
  wrap(p)
}

.read_unit.miranda <- function(x) {
  x <- x[x != ""]
  score <- grpf(x, "Forward.*Score:")
  score <- round(as.double(gs(score, ".*Score:\\s*([0-9.]+).*", "\\1")), 2)
  unit <- "kCal/Mol"
  energy <- gs(grpf(x, "Energy:"), ".*Energy:\\s*([0-9.\\-]+).*", "\\1")
  energy <- round(as.double(energy), 2)
  p.content <- grp(x, "Query:")
  content <- x[p.content:(p.content + 2)]
  content <- strsplit(content, "")
  set <- content[[1]]
  for (i in 1:length(set)) {
    if (set[i] == "3") {
      sig <- i
      break
    }
  }
  content <- lapply(content, function(x) x[sig:length(x)] )
  names(content) <- c("miRNA", "", "Ref")
  content <- as_df.lst(content)
  content <- split_lapply_rbind(content, ~type,
    function(x) {
      x <- dplyr::mutate(x, x = 1:nrow(x))
    })
  content <- dplyr::mutate(content, type = factor(type, levels = c("Ref", "", "miRNA")))
  namel(content, score, energy, unit)
}

.plot_unit.miranda <- function(lst) {
  data <- lst$content
  data <- dplyr::mutate(data, match = ifelse(name %in% LETTERS, "Match", "..."))
  p <- ggplot(data) +
    geom_text(aes(x = x, y = type, label = name, color = match), size = 10) +
    theme_minimal() +
    scale_color_manual(values = c(Match = "red", `...` = "black")) +
    annotate("text", x = 0, y = -1,
      label = paste0("Score: ", lst$score, "\n", "Energy: ", lst$energy, " ", lst$unit),
      size = 6, hjust = 0) +
    annotate("text", x = 0, y = -3, label = "") +
    theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(size = 20)) +
    geom_blank()
  wrap(p, 3 + max(data$x) * .35, 4)
}

.get_mature_miRNA.mirbase <- function(url = "https://mirbase.org/download/mature.fa",
  save = "../mature_miRNA_miRBase.fa")
{
  if (!file.exists(save)) {
    ch <- e(RCurl::getURL(url))
    writeLines(ch, save)
    rm(ch)
  }
  # Biostrings::writeXStringSet
  db <- e(Biostrings::readRNAStringSet(save))
  db
}


