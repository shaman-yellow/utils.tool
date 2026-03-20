# ==========================================================================
# workflow of reactome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_reactome <- setClass("job_reactome", 
  contains = c("job"),
  prototype = prototype(
    pg = "reactome",
    info = c("https://bioconductor.posit.co/packages/devel/bioc/vignettes/ReactomeGSA/inst/doc/analysing-scRNAseq.html"),
    cite = "",
    method = "",
    tag = "reactome",
    analysis = "ReactomeGSA 富集分析"
    ))

setGeneric("asjob_reactome",
  function(x, ...) standardGeneric("asjob_reactome"))

setMethod("asjob_reactome", signature = c(x = "job_seurat"),
  function(x, ref = NULL, split = "orig.ident", group.by = x$group.by)
  {
    object <- object(x)
    metadata <- as_tibble(object@meta.data, idcol = "cell")
    e(Seurat::Idents(object) <- group.by)
    message(glue::glue("Clusters for enrichment: {bind(levels(Seurat::Idents(object)))}"))
    if (is.null(ref)) {
      fea <- as_feature(levels(Seurat::Idents(object)), "关键细胞", "cell")
    }
    snapAdd_onExit("x", "将{snap(fea)}进行通路富集分析。")
    if (!is.null(split)) {
      if (!is.null(object@meta.data[[split]])) {
        message(glue::glue("Split seurat by group of '{split}'."))
        splits <- unique(object@meta.data[[split]])
        lst <- sapply(splits, simplify = FALSE,
          function(name) {
            message(glue::glue("Split get {split}: {name}"))
            e(SeuratObject:::subset.Seurat(object, !!rlang::sym(split) == !!name))
          })
      } else {
        stop('!is.null(object@meta.data[[split]]).')
      }
    } else {
      lst <- list(All = object)
    }
    cli::cli_alert_info("ReactomeGSA::analyse_sc_clusters")
    fun_analysis <- function(name, ...) {
      ReactomeGSA::analyse_sc_clusters(lst[[name]], verbose = TRUE)
    }
    objects <- pbapply::pbsapply(names(lst), simplify = FALSE,
      function(name) {
        # To ensure the uniqueness of the cache, a non functional ID was added
        expect_local_data(
          "tmp", "reactome", fun_analysis, list(
            name = name, id = sig(x)
          )
        )
      }
    )
    x <- .job_reactome(object = objects)
    x$metadata <- metadata
    x$split <- split
    x$from <- "Seurat"
    x <- methodAdd(x, "反应组基因集分析（Reactome Gene Set Analysis，ReactomeGSA）是 Reactome 知识库旗下的多组学通路分析工具，可通过 Camera、PADOG 等方法将单细胞数据转化为通路水平信息，通过计算细胞集群平均表达量实现通路富集分析，识别细胞功能特征，解析细胞异质性。")
    x <- methodAdd(x, "使用 R 包 `ReactomeGSA` ({packageVersion('ReactomeGSA')}) 基于基因表达量对 Reactome 通路进行“活性评分”。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_reactome"),
  function(x){
    step_message("Prepare your data with function `job_reactome`.")
  })

setMethod("step1", signature = c(x = "job_reactome"),
  function(x, gather.by = "group")
  {
    step_message("Gather for test.")
    if (!is.null(x$from) && x$from == "Seurat") {
      if (length(unique(x$metadata[[gather.by]])) != 2) {
        stop('length(unique(x$metadata[[gather.by]])) != 2, only comparison of two group allowed.')
      }
      paths <- lapply(object(x),
        function(obj) {
          ReactomeGSA::get_result(obj, "pathways", "Seurat")
        })
      genes <- lapply(object(x), 
        function(obj) {
          ReactomeGSA::get_result(obj, "fold_changes", "Seurat")
        })
      resTest <- lapply(setNames(list(paths, genes), c("pathways", "genes")),
        function(lst) {
          nrows <- vapply(lst, nrow, integer(1))
          message(
            glue::glue("The results number of each {x$split}: {bind(less(nrows, 10))}")
          )
          data <- dplyr::bind_rows(lst, .id = x$split)
          whichIdCols <- which(vapply(data, is.character, logical(1)))
          if (!all(whichIdCols %in% 1:3)) {
            stop('!all(whichIdCols %in% 1:3).')
          }
          idCols <- colnames(data)[whichIdCols]
          if (!all(idCols %in% c("Pathway", "Name", "Identifier", x$split))) {
            stop('!all(idCols %in% c("Pathway", "Name", "Identifier", x$split)).')
          }
          data <- dplyr::select(data, -dplyr::contains("Pathway"))
          data <- dplyr::rename(data, Name = 2)
          if (!identical(colnames(data)[1:2], c(x$split, "Name"))) {
            stop('!identical(colnames(data)[1:2], c(x$split, "Name")).')
          }
          data <- tidyr::pivot_longer(
            data, dplyr::where(is.double), names_to = "type", values_to = "value"
          )
          data <- map(
            data, x$split, x$metadata, x$split, gather.by, col = gather.by
          )
          # drop the 'sample' columns, for pivot_wider gathered to generate 'list'
          data <- dplyr::select(data, -!!rlang::sym(x$split))
          data <- tidyr::pivot_wider(data, names_from = !!rlang::sym(gather.by), values_from = value)
          stop_debug(namel(data))
          data <- dplyr::group_by(data, Name, type)
          groups <- unique(x$metadata[[gather.by]])
          data <- dplyr::summarise(
            data, test = wilcox.test(.data[[groups[1]]], .data[[groups[2]]])
          )
          data
        })
    }
    return(x)
  })

