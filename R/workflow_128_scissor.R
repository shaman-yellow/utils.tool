# ==========================================================================
# workflow of scissor  (AI) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scissor <- setClass("job_scissor",
  contains = "job",
  prototype = prototype(
    pg = "scissor",
    info = c(
      "https://github.com/sunduanchen/Scissor",
      "https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html"
    ),
    cite = "",
    method = "",
    tag = "scissor",
    analysis = "Scissor 表型关联单细胞亚群鉴定"
  )
)

setGeneric("do_scissor",
  function(x, ref, ...) standardGeneric("do_scissor"))

setMethod("do_scissor", signature = c(x = "job_seurat", ref = "job_deseq2"),
  function(x, ref, group = "group", family = "binomial", assay = SeuratObject::DefaultAssay(object(x)))
  {
    if (x@step < 3L) {
      stop('x@step < 3L.')
    }
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    levels <- NULL
    if (family == "binomial") {
      levels <- rev(.guess_compare_deseq2(ref))
      meta_bulk <- as.integer(factor(object(ref)@colData[[group]], levels = levels)) - 1L
    } else {
      stop('family != "binomial", not yet ready.')
    }
    expr_bulk <- ref$vst@assays@data[[1L]]
    object <- object(x)
    assayObj <- object[[ assay ]]
    if (is(assayObj, "Assay5") || is(assayObj, "SCTAssay")) {
      message(glue::glue("Detected 'Assay5' or 'SCTAssay', convert to 'Assay'."))
      assayObj <- as(assayObj, "Assay")
      object[[ "RNA" ]] <- assayObj
    }
    if (assay != "RNA") {
      # Scissor only extracts "RNA" assay internally!!!
      nameSnn <- glue::glue("{assay}_snn")
      object@graphs[[ "RNA_snn" ]] <- object@graphs[[ nameSnn ]]
      SeuratObject::DefaultAssay(object) <- "RNA"
      object[[ assay ]] <- NULL
    }
    pr <- params(x)
    x <- .job_scissor(object = list(bulk_dataset = expr_bulk, sc_dataset = object, phenotype = meta_bulk))
    x <- methodAdd(x, "Scissor 是一种用于整合单细胞转录组数据与群体水平表型信息（如临床性状或 bulk 转录组数据）的分析方法，其主要目的是在单细胞分辨率下识别与特定表型显著相关的细胞亚群。该方法通过构建单细胞与 bulk 样本之间的表达相关性，并结合回归模型（如惩罚回归）筛选出与目标表型正相关或负相关的细胞群体（分别称为 Scissor+ 和 Scissor− 细胞）。在此基础上，可进一步对相关细胞群进行差异分析与功能注释。")
    x@params <- append(x@params, pr)
    x$family <- family
    x$tag <- levels
    return(x)
  }
)

setMethod("step1", signature = "job_scissor",
  function(x, mode = c("normal", "small", "middle"), rerun = FALSE,
    dir_cache = create_job_cache_dir(x), ...)
  {
    step_message("Run Scissor algorithm.")
    file_save <- file.path(dir_cache, "scissor_inputs.rdata")
    mode <- match.arg(mode)
    if (mode == "small") {
      alphas <- c(5e-3, 5e-4, 1e-4, 5e-5)
    } else if (mode == "middle") {
      alphas <- c(3e-3, 1e-3, 7e-4, 9e-4)
    } else if (mode == "normal") {
      alphas <- c(.005, .05, .1, .2)
    }
    fun_scissor <- function(...) {
      if (file.exists(file_save)) {
        Load_file <- file_save
      } else {
        Load_file <- NULL
      }
      options(future.globals.maxSize = Inf)
      old_plan <- future::plan()
      future::plan(future::multicore, workers = length(alphas))
      on.exit(future::plan(old_plan))
      res_lst <- e(future.apply::future_lapply(
        alphas,
        function(alpha) {
          Scissor::Scissor(
            object(x)$bulk_dataset, object(x)$sc_dataset, object(x)$phenotype,
            Save_file = file_save, family = x$family, 
            Load_file = Load_file,
            tag = x$tag, alpha = alpha
          )
        }
      ))
      setNames(res_lst, paste0("alpha_", alphas))
    }
    x$res_scissor <- expect_local_data(
      "tmp", "scissor", fun_scissor, rerun = rerun,
      list(
        colnames(object(x)$bulk_dataset),
        colnames(object(x)$sc_dataset),
        object(x)$phenotype,
        alphas
      )
    )
    return(x)
  }
)

