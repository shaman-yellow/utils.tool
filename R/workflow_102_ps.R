# ==========================================================================
# workflow of ps
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ps <- setClass("job_ps", 
  contains = c("job"),
  prototype = prototype(
    pg = "ps",
    info = c("http://db.phasep.pro/"),
    cite = "[@PhaSepDB_a_dat_You_K_2020]",
    method = "",
    tag = "ps",
    analysis = "PhaSepDB 相分离相关蛋白"
    ))

job_ps <- function(dir = .prefix("llps", "db"), distinct = TRUE)
{
  dir.create(dir, FALSE)
  url <- "http://db.phasep.pro/static/db/database/phasepdbv2_1_llps.xlsx"
  file <- file.path(dir, basename(url))
  if (!file.exists(file)) {
    utils::download.file(url, file)
  }
  object <- fxlsx(file)
  object <- dplyr::mutate(
    object, location = ifelse(
      location == "_", "Unknown", location
    ), .before = 1
  )
  if (distinct) {
    object <- dplyr::distinct(object, gene_name, .keep_all = TRUE)
  }
  object <- set_lab_legend(object, "LLPS all data", "为 PhaSepDB LLPS 所有数据附表。")
  p.location <- new_pie(object$location, title = "Location")
  p.location <- set_lab_legend(p.location, "LLPS Location", "为液液相分离 (LLPS) 蛋白的位置分布饼图。")
  p.org <- new_pie(
    object$organism, title = "Organism", legend = FALSE, legend_ncol = NULL
  )
  p.org <- set_lab_legend(p.org, "LLPS organism", "为液液相分离 (LLPS) 蛋白的来源分布饼图。")
  x <- .job_ps(object = object)
  x <- methodAdd(x, "获取 PhaSepDB ({cite_show('PhaSepDB_a_dat_You_K_2020')}) LLPS 蛋白数据。")
  x$p.location <- p.location
  x$p.org <- p.org
  return(x)
}

setMethod("step0", signature = c(x = "job_ps"),
  function(x){
    step_message("Prepare your data with function `job_ps`.")
  })

setMethod("step1", signature = c(x = "job_ps"),
  function(x, org = "Homo sapiens"){
    step_message("Quality control (QC).")
    org <- match.arg(org)
    db <- dplyr::filter(object(x), organism == org)
    x$db <- set_lab_legend(db, paste("LLPS", org), glue::glue("为 LLPS 数据 {org} 来源附表。"))
    feature(x) <- unique(unlist(strsplit(x$db$Gene_names, " ")))
    x <- snapAdd(x, "使用来源为 {org} 的条目。(注：仅根据基因名，合并这些条目，共包含 Symbol {length(x$.feature)} 个。) ")
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_ps"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
