# ==========================================================================
# workflow of hybrid
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_hybrid <- setClass("job_hybrid", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c(""),
    cite = "[@FastAndEffectRehmsm2004; @CircbaseADatGlazar2014]",
    method = paste0("`RNAhybrid` used as a means for microRNA target prediction (circRNA-miRNA) and ",
      "`circBase` used for querying circRNA sequences")
    ))

job_hybrid <- function(mirna, circrna, circFa = NULL)
{
  message("The `mirna` can be character (pattern) for match miRNA.\n",
    "The `circrna` can be character (numbers) (pattern) for match circRNA."
  )
  mi <- .prepare_mirna(mirna)
  ci <- .prepare_circrna(circrna, circFa)
  x <- .job_hybrid(object = list(mirna = mi$mirna, circrna = ci$circrna))
  x$mi <- mi
  x$ci <- ci
  x
}

setMethod("step0", signature = c(x = "job_hybrid"),
  function(x){
    step_message("Prepare your data with function `job_hybrid`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_hybrid"),
  function(x, wd = timeName("hybrid")){
    step_message("Prepare computional files.")
    dir.create(wd, F)
    x$wd <- wd
    # miRNA seq
    x$mirna_seq <- x$mi$mirdb[ object(x)$mirna ]
    x$mirna_file <- paste0(wd, "/", make.names(object(x)$mirna), ".fa")
    e(Biostrings::writeXStringSet(x$mirna_seq, x$mirna_file), text = "miRNA")
    # circRNA seq
    x$circrna_seq <- x$ci$circdb[ object(x)$circrna ]
    x$circrna_file <- paste0(wd, "/", make.names(object(x)$circrna), ".fa")
    e(Biostrings::writeXStringSet(x$circrna_seq, x$circrna_file), text = "circRNA")
    return(x)
  })

setMethod("step2", signature = c(x = "job_hybrid"),
  function(x, command = "conda run -n base RNAhybrid"){
    step_message("Run RNAhybrid.")
    x$outfile <- paste0(x$wd, "/", "results.txt")
    cdRun(command,
      " -q ", x$mirna_file,
      " -t ", x$circrna_file,
      " -s 3utr_human ",
      " -b 10 ",
      " -e -5 ",
      " > ", x$outfile
    )
    x$alls <- .read_hybrid_results(x$outfile)
    t.overs <- lapply(x$alls$cans,
      function(x) {
        x <- lapply(x, function(x) data.frame(pvalue = x$pvalue, energy = x$energy))
        data.table::rbindlist(x)
      })
    t.overs <- as_tibble(data.table::rbindlist(t.overs, idcol = T))
    t.overs <- dplyr::relocate(t.overs, ref = .id)
    p.tops <- lapply(x$alls$cans,
      function(lst) {
        p <- .plot_unit.hybrid(lst[[1]])
      })
    p.tops <- .set_lab(p.tops, sig(x), names(p.tops), "Top match")
    p.overs <- .plot_energy.hybrid(t.overs)
    p.overs <- .set_lab(p.overs, sig(x), "overview of all candidates")
    x@tables[[2]] <- namel(t.overs)
    x@plots[[2]] <- namel(p.tops, p.overs)
    return(x)
  })

.plot_energy.hybrid <- function(data) {
  data <- split_lapply_rbind(data, ~ref,
    function(x) {
      x <- dplyr::mutate(x, candidate = paste0("C ", 1:nrow(x)))
    })
  p <- ggplot(data) +
    geom_col(aes(x = reorder(candidate, energy, decreasing = T), y = energy, fill = pvalue), width = .7) +
    facet_wrap(~ ref) +
    coord_flip(ylim = zoRange(data$energy, 1.2)) +
    labs(x = "Energy", y = "Candidates", fill = "P-value") +
    scale_fill_gradientn(colors = rev(ggsci::pal_npg()(3))) +
    theme_minimal()
  if (length(unique(data$ref)) == 1) {
    return(wrap(p, 8, 1 + nrow(data) * .5))
  }
  wrap(p)
}

.plot_unit.hybrid <- function(lst) {
  data <- lst$content
  p <- ggplot(data) +
    geom_text(aes(x = x, y = type, label = name, color = type), size = 10) +
    scale_color_manual(values = c(` ` = "red", `  ` = "red", Target = "black", miRNA = "black")) +
    annotate("text", x = 0, y = -1,
      label = paste0("P-value: ", lst$pvalue, "\n", "Energy: ", lst$energy, " ", lst$unit, "\n",
        "Target start position: ", lst$position),
      size = 6, hjust = 0) +
    annotate("text", x = 0, y = -3, label = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(size = 20)) +
    geom_blank()
  wrap(p, 3 + max(data$x) * .35, 4)
}

.read_hybrid_results <- function(file) {
  body <- readLines(file)
  cans <- sep_list(body, "^target:", T)
  cans <- lapply(cans, .read_unit.hybrid)
  belong <- vapply(cans, function(x) x$target, character(1))
  cans <- split(cans, make.names(belong))
  cans <- lst_clear0(cans)
  namel(cans)
}

.read_unit.hybrid <- function(x) {
  x <- x[x != ""]
  target <- stringr::str_extract(grpf(x, "^target:"), "[^ ]+$")
  pvalue <- grpf(x, "p-value:")
  pvalue <- round(as.double(gs(pvalue, "p-value:\\s*([0-9.]+).*", "\\1")), 4)
  position <- grpf(x, "^position")
  position <- as.integer(gs(position, "position\\s*([0-9.]+).*", "\\1"))
  unit <- "kCal/Mol"
  energy <- gs(grpf(x, "mfe:"), "mfe:\\s*([0-9.\\-]+).*", "\\1")
  energy <- round(as.double(energy), 2)
  p.content <- grp(x, "target 5'")
  content <- x[p.content:(p.content + 3)]
  content <- strsplit(content, "")
  set <- content[[4]]
  for (i in 1:length(set)) {
    if (set[i] == "3") {
      sig <- i
      break
    }
  }
  content <- lapply(content, function(x) x[sig:length(x)] )
  names(content) <- .type <- c("Target", "  ", " ", "miRNA")
  content <- as_df.lst(content)
  content <- split_lapply_rbind(content, ~type,
    function(x) {
      x <- dplyr::mutate(x, x = 1:nrow(x))
    })
  content <- dplyr::mutate(content, type = factor(type, levels = !!.type))
  namel(target, content, pvalue, energy, unit, position)
}

.prepare_circrna <- function(circrna, customdb = NULL) {
  if (is.null(customdb)) {
    circdb <- .get_spliced_circRNA.circBase()
    if (grpl(circrna, "^[0-9]*$")) {
      message("Match ID number.")
      circrna <- paste0("_[0]*", circrna, "[^0-9]")
    }
  } else {
    message("Use custom database (fasta) for circRNA.")
    circdb <- e(Biostrings::readRNAStringSet(customdb))
  }
  foundThat <- grpf(names(circdb), circrna, ignore.case = T)
  message("All found miRNA:\n\n", stringr::str_trunc(paste0(foundThat, collapse = "\n"), 100), "\n")
  if (length(foundThat) > 1) {
    message("Use the first.")
  }
  circrna <- foundThat[1]
  circdb <- circdb[foundThat]
  namel(circrna, foundThat, circdb)
}

.get_spliced_circRNA.circBase <- function(
  url = "http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz",
  save = "../splice_circRNA_circBase.fa.gz")
{
  if (!file.exists(save)) {
    cdRun("wget ", url, " -O ", save)
  }
  # Biostrings::writeXStringSet
  db <- e(Biostrings::readRNAStringSet(save))
  db
}
