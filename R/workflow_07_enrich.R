# ==========================================================================
# workflow of enrich
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_enrich <- setClass("job_enrich", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("...")
    ))

job_enrich <- function(ids, annotation, from = "hgnc_symbol", to = "entrezgene_id")
{
  if (is.null(names(ids)))
    stop("is.null(names(ids))")
  if (missing(annotation))
    stop("missing(annotation)")
  maps <- lapply(ids,
    function(id) {
      unique(filter(annotation, !!rlang::sym(from) %in% !!id)[[ to ]])
    })
  en <- .job_enrich(object = maps)
  en@params$raw <- ids
  en@params$annotation <- annotation
  en
}

setMethod("step0", signature = c(x = "job_enrich"),
  function(x){
    step_message("Prepare your data with function `job_enrich`.
      "
    )
  })

# Biological Process, Molecular Function, and Cellular Component groups
setMethod("step1", signature = c(x = "job_enrich"),
  function(x, orgDb = 'org.Hs.eg.db', cl = 4, maxShow = 10){
    step_message("Use clusterProfiler for enrichment.
      "
    )
    cli::cli_alert_info("clusterProfiler::enrichKEGG")
    res.kegg <- multi_enrichKEGG(object(x))
    p.kegg <- vis_enrich.kegg(res.kegg, maxShow = maxShow)
    p.kegg <- lapply(p.kegg, function(x) wrap(x, 8, 4 * (maxShow / 10)))
    fun <- function(sets) {
      lapply(sets,
        function(set) {
          from_ids <- x@params$annotation$entrezgene_id
          to_names <- x@params$annotation$hgnc_symbol
          to_names[ match(set, from_ids) ]
        })
    }
    res.kegg <- lapply(res.kegg, mutate, geneName_list = fun(geneID_list))
    cli::cli_alert_info("clusterProfiler::enrichGO")
    res.go <- multi_enrichGO(object(x), orgDb = orgDb, cl = cl)
    p.go <- vis_enrich.go(res.go, maxShow = maxShow)
    x@tables[[ 1 ]] <- namel(res.kegg, res.go)
    x@plots[[ 1 ]] <- namel(p.kegg, p.go)
    x@params$check_go <- check_enrichGO(res.go)
    return(x)
  })

setMethod("step2", signature = c(x = "job_enrich"),
  function(x){

  })
