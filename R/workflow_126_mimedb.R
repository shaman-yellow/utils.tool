# ==========================================================================
# workflow of mimedb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mimedb <- setClass("job_mimedb", 
  contains = c("job"),
  prototype = prototype(
    pg = "mimedb",
    info = c("https://mimedb.org/downloads"),
    cite = "",
    method = "",
    tag = "mimedb",
    analysis = "MiMeDB 代谢物数据挖掘"
    ))

job_mimedb <- function(patterns)
{
  data <- ftibble(get_url_data(
    "mimedb_metabolites_v2.csv",
    "https://mimedb.org/system/downloads/2.0/mimedb_metabolites_v2.csv",
    "mimedb_metabolites_v2", dir = .prefix("mimedb", "db"), fun_decompress = NULL
  ))
  data <- dplyr::filter(
    data, grpl(
      name, patterns, TRUE
    ) | grpl(description, patterns, TRUE)
  )
  if (!nrow(data)) {
    stop('!nrow(data).')
  } else {
    message(glue::glue("Got data ({nrow(data)}): "))
    print(data)
  }
  x <- .job_mimedb()
  x$ids <- setNames(data$mime_id, data$name)
  x$data <- data
  return(x)
}

setMethod("step0", signature = c(x = "job_mimedb"),
  function(x){
    step_message("Prepare your data with function `job_mimedb`.")
  })

setMethod("step1", signature = c(x = "job_mimedb"),
  function(x, dir = "~/Downloads", which = 1L)
  {
    step_message("Collate data.")
    used <- c("Human_Proteins_and_Enzymes", "Microbial_Sources")
    files <- list.files(dir, "MiMeDB.*.html", full.names = TRUE)
    id_files <- strx(files, "MMDBc[0-9]+")
    files <- setNames(files, id_files)
    lst <- sapply(x$ids, simplify = FALSE,
      function(id) {
        lst <- .format_mimedb_html_data(files[id])
        lst <- lst[ names(lst) %in% used ]
        lst$Human_Proteins_and_Enzymes <- dplyr::select(
          lst$Human_Proteins_and_Enzymes, 
          `Protein Type`, `Protein Name`, `Uniprot ID`
        )
        lst$Human_Proteins_and_Enzymes <- dplyr::arrange(
          lst$Human_Proteins_and_Enzymes, `Protein Type`
        )
        lst$Human_Proteins_and_Enzymes <- dplyr::group_by(
          lst$Human_Proteins_and_Enzymes, `Protein Type`
        )
        lst$Microbial_Sources <- dplyr::select(
          lst$Microbial_Sources, 
          `Superkingdom/Kingdom`, Phylum, `Total species`, `Hosts and Body Sites`
        )
        lst$Microbial_Sources <- dplyr::arrange(
          lst$Microbial_Sources, `Superkingdom/Kingdom`, Phylum
        )
        lst$Microbial_Sources <- dplyr::group_by(
          lst$Microbial_Sources, `Superkingdom/Kingdom`
        )
        return(lst)
      })
    names(lst) <- formal_name(names(x$ids))
    ts.proteins <- lapply(lst, function(x) x$Human_Proteins_and_Enzymes)
    ts.proteins <- set_lab_legend(
      ts.proteins,
      glue::glue("{x@sig} {names(ts.proteins)} Human_Proteins_and_Enzymes"),
      glue::glue("人类中相关蛋白与酶")
    )
    ts.microbe <- lapply(lst, function(x) x$Microbial_Sources)
    ts.microbe <- set_lab_legend(
      ts.microbe,
      glue::glue("{x@sig} {names(ts.microbe)} Microbial_Sources"),
      glue::glue("代谢物来源微生物")
    )
    x <- snapAdd(x, "搜集整理 {bind(names(x$ids))} 在 MiMeDB 数据库中的微生物来源和相关基因或蛋白数据记录。")
    if (length(ts.microbe) == 1) {
      num <- length(unique(ts.proteins[[1]]$`Protein Name`))
      x <- snapAdd(x, "表格{aref(ts.proteins)}为 {names(x$ids)} 在人体中的相关蛋白，共包含 {num} 种不同的蛋白。")
      x <- snapAdd(x, "表格{aref(ts.microbe)}为 {names(x$ids)} 在人体中的来源微生物，共包含 {nrow(ts.microbe[[1]])} 个条目。")
    }
    x <- methodAdd(x, "MiMeDB (<https://mimedb.org/>) 可用于系统性挖掘代谢物相关的蛋白及其微生物来源，其核心在于整合“代谢物–微生物–酶/蛋白–宿主效应”的多组学关联信息。该数据库不仅提供代谢物的结构与分类信息，还明确记录其来源微生物、参与代谢反应的酶（蛋白）以及相关基因信息，从而支持从代谢物反向追溯其生物来源及作用机制的分析。")
    x <- tablesAdd(x, ts.microbe, ts.proteins)
    return(x)
  })

.format_mimedb_html_data <- function(file) {
  lst <- get_table.html(readLines(file))
  lst <- lapply(lst, 
    function(data) {
      data <- as_tibble(data)
      if (colnames(data)[1] == "V1") {
        if (all(!is.na(names <- unlist(data[1, ])))) {
          data <- setNames(data[-1, ], names)
        }
      }
      data
    })
  colFirst <- vapply(lst, function(x) colnames(x)[1], character(1))
  names <- dplyr::recode(
    colFirst, "Superkingdom/Kingdom" = "Microbial Sources",
    "Health Outcome/Bioactivity" = "Disease",
    "Protein ID" = "Human Proteins and Enzymes"
  )
  setNames(lst, formal_name(names))
}

