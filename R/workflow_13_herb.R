# ==========================================================================
# workflow of herb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_herb <- setClass("job_herb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@HerbAHighThFang2021]"
    ))

job_herb <- function(herbs, db = get_herb_data())
{
  x <- .job_herb(object = db)
  x@params$herbs <- herbs
  herbs_info <- filter(object(x)$herb, Herb_cn_name %in% !!params(x)$herbs)
  x@params$herbs_info <- herbs_info
  print(herbs_info)
  x
}

setMethod("step0", signature = c(x = "job_herb"),
  function(x){
    step_message("Prepare your data with function `job_herb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_herb"),
  function(x){
    step_message("Dowload compounds of herbs.
      "
    )
    link <- start_drive(browser = "firefox")
    Sys.sleep(3)
    tryCatch(
      download_herbCompounds(link, params(x)$herbs_info$Herb_),
      finally = end_drive()
    )
    if (file.exists("herbs_ingredient")) {
      unlink("herbs_ingredient", T, T)
    }
    herbs_compounds <- moveToDir_herbs(params(x)$herbs_info$Herb_)
    x@tables[[ 1 ]] <- namel(herbs_compounds)
    return(x)
  })

setMethod("step2", signature = c(x = "job_herb"),
  function(x, group_number = 50, group_sleep = 3)
  {
    step_message("Dowload targets of compounds")
    link <- start_drive(browser = "firefox")
    Sys.sleep(3)
    all_ids <- unique(x@tables$step1$herbs_compounds$Ingredient.id)
    all_ids <- grouping_vec2list(all_ids, group_number, T)
    lst <- pbapply::pblapply(all_ids,
      function(ids) {
        to <- paste0("compounds_target", "/", attr(ids, "name"))
        if (!file.exists(to)) {
          download_compoundTargets(link, ids)
        }
        data <- moveToDir_herbs(ids, to = to, .id = "Ingredient_id")
        cli::cli_alert_info("Sleep between groups ...")
        Sys.sleep(group_sleep)
        data
      })
    end_drive()
    compounds_targets <- as_tibble(data.table::rbindlist(lst, fill = T))
    x@tables[[ 2 ]] <- namel(compounds_targets)
    return(x)
  })

setMethod("step3", signature = c(x = "job_herb"),
  function(x, disease = NULL, mart_dataset = "hsapiens_gene_ensembl")
  {
    step_message("Get disease targets from genecards and annotate targets with biomaRt.
      As well, plot the intersection of herbs targets.
      "
    )
    if (is.null(x@params$mart)) {
      mart <- new_biomart(mart_dataset)
      x@params$mart <- mart
    } else {
      mart <- x@params$mart
    }
    ## annotate compounds targets
    compounds_targets <- x@tables$step2$compounds_targets
    targets <- unique(compounds_targets$Target.name)
    targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol", targets)
    x@params$targets_annotation <- targets_annotation
    compounds_targets_annotation <- tbmerge(
      compounds_targets, targets_annotation,
      by.x = "Target.name", by.y = "hgnc_symbol", all.x = T
    )
    ## create data: herbs_target
    herbs_compounds <- x@tables$step1$herbs_compounds
    herbs_targets <- tbmerge(
      herbs_compounds, compounds_targets,
      by.x = "Ingredient.id", by.y = "Ingredient_id", all.x = T,
      allow.cartesian = T
    )
    x@params$herbs_targets <- herbs_targets
    easyRead <- map(x@params$herbs_targets, "herb_id", object(x)$herb, "Herb_", "Herb_pinyin_name")
    x@params$easyRead <- easyRead
    ## plot upset of herbs targets
    sets <- distinct(herbs_targets, herb_id, Target.name)
    sets <- split(sets$Target.name, sets$herb_id)
    fun <- function(sets) x@params$herbs_info$Herb_pinyin_name[match(names(sets), x@params$herbs_info$Herb_)]
    names(sets) <- fun(sets)
    p.herbs_targets <- new_upset(lst = sets)
    ## plot upset of herbs ingredient
    sets <- split(herbs_compounds$Ingredient.id, herbs_compounds$herb_id)
    names(sets) <- fun(sets)
    p.herbs_compounds <- new_upset(lst = sets)
    ## set plots and tables
    plots <- namel(p.herbs_targets, p.herbs_compounds)
    tables <- namel(compounds_targets_annotation)
    if (!is.null(disease)) {
      disease_targets <- get_from_genecards(disease)
      x@params$disease_targets <- disease_targets
      disease_targets_annotation <- filter_biomart(mart, general_attrs(), "hgnc_symbol",
        disease_targets$Symbol)
      tables <- c(tables, namel(disease_targets_annotation))
      p.targets_disease <- new_upset(
        targets = targets,
        targets_annotation = targets_annotation$hgnc_symbol,
        disease_targets = disease_targets$Symbol,
        disease_targets_annotation = disease_targets_annotation$hgnc_symbol
      )
      plots <- c(plots, namel(p.targets_disease))
    }
    x@plots[[ 3 ]] <- plots
    x@tables[[ 3 ]] <- tables
    x@params$ppi_used <- dplyr::filter(
      targets_annotation,
      hgnc_symbol %in% dplyr::all_of(disease_targets_annotation$hgnc_symbol)
    )
    return(x)
  })



