# ==========================================================================
# workflow of tfbs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_tfbs <- setClass("job_tfbs", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://tfbsdb.systemsbiology.net/"),
    cite = "[@CausalMechanisPlaisi2016]",
    method = "The `Transcription Factor Target Gene Database` (<https://tfbsdb.systemsbiology.net/>) was used for discovering relationship between transcription factors and genes. ",
    tag = "tf",
    analysis = "Transcription Factor Target Gene Database"
    ))

job_tfbs <- function(genes)
{
  .job_tfbs(object = rm.no(genes))
}

setMethod("step0", signature = c(x = "job_tfbs"),
  function(x){
    step_message("Prepare your data with function `job_tfbs`.")
  })

setMethod("step1", signature = c(x = "job_tfbs"),
  function(x, db = .prefix("tfbs/tfbs.rds", "db"), db_anno = .prefix("tfbs/anno.rds", "db"),
    file_anno = .prefix("tfbs/humanTFMotifEntrezMappings.xlsx", "db"), cl = NULL)
  {
    step_message("Obtain Transcription Factor.")
    db <- new_db(db, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      res <- pbapply::pbsapply(db@query, simplify = FALSE,
        function(gene) {
          url <- paste0("https://tfbsdb.systemsbiology.net/searchgene?searchterm=", gene)
          res <- RCurl::getURL(url)
          res <- try(get_table.html(res), TRUE)
          if (inherits(res, "try-error")) {
            message("\nCan not get results of ", gene, ", continue to next.")
            return()
          }
          res <- lapply(res,
            function(x) {
              if (is(x, "df"))
                as_tibble(x)
            })
          res$tfbs
        })
      res <- frbind(res, idcol = TRUE, fill = TRUE)
      db <- upd(db, res)
    }
    res <- dplyr::filter(db@db, .id %in% object(x))
    if (TRUE) {
      if (!file.exists(file_anno)) {
        cont <- e(RCurl::getURLContent(
            "https://tfbsdb.systemsbiology.net/static/downloads/humanTFMotifEntrezMappings.xlsx"
            ))
        writeBin(as.vector(cont), file_anno)
      }
      anno <- fxlsx(file_anno)
      res <- map(res, "Motif", anno, "Motif.Name", "Gene", col = "Symbol")
      if (any(is.na(res$Symbol))) {
        ids <- unique(res$Motif[ is.na(res$Symbol) ])
        db <- new_db(db_anno, ".id")
        db <- not(db, ids)
        if (length(db@query)) {
          message("Need query: ", length(db@query), " (total: ", length(unique(res$Motif)), ")")
          anno <- pbapply::pbsapply(db@query, simplify = FALSE, cl = cl,
            function(id) {
              url <- paste0("https://tfbsdb.systemsbiology.net/searchtf?searchterm=", id)
              res <- try_get_url(url)
              lst <- suppressMessages(try(get_table.html(res, which = 1:2), TRUE))
              if (inherits(lst, "try-error")) {
                message("Skip: ", id, " as error occurred.")
                return(data.frame(Symbol = NA_character_))
              }
              res <- lapply(lst,
                function(x) {
                  if (is(x, "df"))
                    as_tibble(x)
                })
              des <- strx(res[[2]][1,1]$V1, "^[^\n]+")
              Sys.sleep(1)
              data.frame(Symbol = strx(des, "^[^ ]+"), des = des)
            })
          anno <- frbind(anno, idcol = TRUE, fill = TRUE)
          db <- upd(db, anno)
        }
        anno <- dplyr::filter(db@db, .id %in% ids)
        res$Symbol[ is.na(res$Symbol) ] <- anno$Symbol[match(res$Motif[ is.na(res$Symbol) ], anno$.id)]
      }
      if (TRUE) {
        res$Symbol[ is.na(res$Symbol) ] <- toupper(strx(res$Motif[is.na(res$Symbol)], "[^_]{2,}"))
      }
    }
    res <- dplyr::relocate(res, target = .id, TF_symbol = Symbol)
    res <- .set_lab(res, sig(x), "Transcription Factor binding sites")
    x@tables[[ 1 ]] <- namel(res)
    return(x)
  })

setMethod("map", signature = c(x = "job_tfbs", ref = "character"),
  function(x, ref, axes = 1:3, ...)
  {
    data <- dplyr::select(x@tables$step1$res, target, TF_symbol, Motif, PValue)
    data <- dplyr::mutate(data, PValue = as.double(PValue))
    data <- split_lapply_rbind(data, seq_len(nrow(data)), args = list(fill = TRUE),
      function(x) {
        if (grpl(x$TF_symbol, "/")) {
          symbols <- strsplit(x$TF_symbol, "/")[[1]]
          x <- do.call(dplyr::bind_rows, rep(list(x), length(symbols)))
          x$TF_symbol <- symbols
        }
        x
      })
    p.venn <- new_venn(DB_TFs = data$TF_symbol, Candidates = ref)
    p.venn <- .set_lab(p.venn, sig(x), "intersection of TFs with queried Candidates")
    data <- dplyr::filter(data, TF_symbol %in% ref)
    p.allu <- new_allu(data, axes = axes, col.fill = 4, ...) +
      theme_minimal() +
      labs(fill = "P-value") +
      theme(axis.text.y = element_blank(),
        axis.title = element_blank()) +
      geom_blank()
    p.allu <- .set_lab(wrap(p.allu, 10, 7), sig(x), "The genes and related TFs")
    x$mapped <- namel(p.venn, p.allu, data)
    return(x)
  })
