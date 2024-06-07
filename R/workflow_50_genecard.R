# ==========================================================================
# workflow of genecard
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_genecard <- setClass("job_genecard", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://www.genecards.org/"),
    cite = "[@TheGenecardsSStelze2016]",
    method = "The Human Gene Database `GeneCards` used for disease related genes prediction",
    tag = "dis:genecards",
    analysis = "GeneCards 基因获取"
    ))

job_genecardn <- function(...) {
  .job_genecardn(list(...))
}

.job_genecardn <- setClass("job_genecardn",
  contains = c("list")
)

setValidity("job_genecardn", 
  function(object){
    if (all(vapply(object, function(x) is(x, "job_genecard"), logical(1))))
      T
    else F
  })

job_genecard <- function(disease)
{
  .job_genecard(object = disease)
}

setMethod("step0", signature = c(x = "job_genecard"),
  function(x){
    step_message("Prepare your data with function `job_genecard`.")
  })

setMethod("step1", signature = c(x = "job_genecard"),
  function(x, score = NULL, restrict = F)
  {
    step_message("Get from GeneCards website.")
    t.genecards <- get_from_genecards(object(x), score = score, restrict = restrict)
    t.genecards <- .set_lab(t.genecards, sig(x), "disease related targets from GeneCards")
    x@tables[[ 1 ]] <- namel(t.genecards)
    return(x)
  })

get_from_genecards <- function(query, score = 5, keep_drive = F, restrict = F,
  advance = F, term = c("compounds"))
{
  link <- start_drive(browser = "firefox")
  Sys.sleep(3)
  link$open()
  query_raw <- query
  if (grpl(query, " ")) {
    query <- gs(query, " ", "%20")
  }
  if (restrict) {
    query <- paste0("\"", query, "\"")
  }
  if (advance) {
    term <- match.arg(term)
    query <- paste0("%20%5B", term, "%5D%20%20(%20%20", query, "%20%20)%20")
  }
  url <- paste0('https://www.genecards.org/Search/Keyword?queryString=', query,
    '&pageSize=25000&startPage=0')
  message("Navigate to URL:\n\t", url)
  link$navigate(url)
  Sys.sleep(5)
  html <- link$getPageSource()[[1]]
  link$close()
  end_drive()
  html <- XML::htmlParse(html)
  table <- XML::readHTMLTable(html)
  table <- as_tibble(data.frame(table[[1]]))
  colnames(table) %<>% gs("\\.+", "_")
  colnames(table) %<>% gs("X_|_$", "")
  table <- dplyr::select(table, -1, -2)
  table <- dplyr::mutate(table, Score = as.double(Score))
  if (is.null(score)) {
    ss <- c(seq(1, 30, by = 1), 0L)
    chs <- vapply(ss, FUN.VALUE = double(1),
      function(x) {
        length(which(table$Score > x))
      })
    which <- menu(paste0("Score ", ss, " Gene ", chs))
    score <- ss[ which ]
  }
  table <- filter(table, Score > !!score)
  attr(table, ".score") <- score
  lich <- list(
    "The GeneCards data was obtained by querying" = query_raw,
    "Restrict (with quotes)" = restrict,
    "Filtering by Score:" = paste0("Score > ", score)
  )
  if (advance) {
    lich <- c(lich, list("Advance search:" = paste0(" [", term, "]  (  ", query_raw, "  ) ")))
  }
  attr(table, "lich") <- new_lich(lich)
  return(table)
}

setMethod("map", signature = c(x = "job_genecard", ref = "job_genecard"),
  function(x, ref, names = c(x@object, ref@object)){
    lst <- lapply(list(x, ref), function(x) x@tables$step1$t.genecards$Symbol)
    names(lst) <- names
    new_venn(lst = lst)
  })

setMethod("cal_corp", signature = c(x = "job_limma", y = "job_genecardn"),
  function(x, y, from, names = NULL, use = if (x$isTcga) "gene_name" else "hgnc_symbol",
    HLs = NULL, ..., top = NULL, tidyheatmap = T, linear = F)
  {
    genes.to <- lapply(y,
      function(x) {
        if (!is.null(top)) {
          genes <- head(x@tables$step1$t.genecards$Symbol, n = top)
        } else {
          genes <- x@tables$step1$t.genecards$Symbol         
        }
        nl(x@object, list(genes))
      })
    genes.to <- unlist(genes.to, recursive = F)
    dat.to <- as_df.lst(genes.to)
    datUni.to <- dplyr::distinct(dat.to, name, .keep_all = T)
    if (length(y) == 1 && is.null(names)) {
      names <- c("From", make.names(y[[1]]@object))
    }
    ## cal_corp ...
    lst.cor <- cal_corp(x, NULL, from, datUni.to$name, names = names, use = use, ...)
    if (tidyheatmap) {
      data <- as_tibble(lst.cor$corp)
      if (!is.null(names)) {
        from.name <- names[[1]]
        to.name <- names[[2]]
      } else {
        from.name <- "From"
        to.name <- "To"
      }
      data <- map(data, to.name, datUni.to, "name", "type", col = "Factors")
      if (!is.null(HLs)) {
        data <- dplyr::mutate(data, Type.from = ifelse(!!rlang::sym(from.name) %in% HLs, "Key", "Others"))
        # data <- dplyr::mutate(data, Type.to = ifelse(!!rlang::sym(to.name) %in% HLs, "Key", "Others"))
      }
      data <- dplyr::group_by(data, Factors, Type.from)
      tihp.cor <- new_hp.cor(data, names = c(from.name, to.name))
      lst.cor$tihp.cor <- tihp.cor
    }
    if (linear) {
      data <- dplyr::filter(data, !!rlang::sym(from.name) %in% !!HLs, pvalue < .05)
      linear.cor <- cal_corp(x, NULL, data[[ from.name ]], data[[ to.name ]], use = use, mode = "l")
      linear.cor <- map(linear.cor, to.name, datUni.to, "name", "type", col = "Factors")
      lp.cor <- vis(.corp(linear.cor), "Factors")
      lst.cor$lp.cor <- lp.cor
    }
    lst.cor
  })

setMethod("res", signature = c(x = "job_genecard"),
  function(x){
    x@tables$step1$t.genecards
  })

plot_col.genecard <- function(data, top = 10, facet = T)
{
  data <- head(data, top)
  p <- ggplot(data, aes(x = reorder(Symbol, Score), y = Score, fill = Score)) +
    geom_col(width = .5) +
    geom_text(aes(x = Symbol, y = Score + max(Score) * .01, label = Score), hjust = 0, size = 3) +
    ylim(c(0, max(data$Score) * 1.2)) +
    coord_flip() +
    labs(x = "Gene Symbol", y = "Score") +
    rstyle("theme") +
    theme(legend.position = "") +
    if (facet) {
      facet_grid(rows = ggplot2::vars(Category), scales = "free_y")
    } else {
      geom_blank()
    }
  p <- wrap(p, 7, nrow(data) * .4 + .5)
  p <- .set_lab(p, "Genecard Score visualization")
  p
}
