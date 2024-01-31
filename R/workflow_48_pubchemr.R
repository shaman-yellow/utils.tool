# ==========================================================================
# workflow of pubchemr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_pubchemr <- setClass("job_pubchemr", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/satijalab/pubchemr/wiki"),
    cite = "",
    method = "R package `PubChemR` used for querying compounds information"
    ))

job_pubchemr <- function(cids)
{
  cids <- cids[ !is.na(cids) ]
  x <- .job_pubchemr(object = cids)
  x
}

setMethod("step0", signature = c(x = "job_pubchemr"),
  function(x){
    step_message("Prepare your data with function `job_pubchemr`.")
  })

setMethod("step1", signature = c(x = "job_pubchemr"),
  function(x){
    step_message("Get synonyms and smiles of compounds.")
    synonyms <- e(PubChemR::get_synonyms(object(x)))
    props <- e(PubChemR::get_properties(c("IsomericSMILES", "InChIKey"), object(x)))
    fun_format <- function(x, get) {
      names <- vapply(x, function(x) x[[ "CID" ]], numeric(1))
      x <- lapply(x,
        function(x) {
          if (is(x, "list")) {
            x[[ get ]]
          } else NULL
        })
      names(x) <- names
      return(x)
    }
    x$smiles <- fun_format(props, "IsomericSMILES")
    x$inks <- fun_format(props, "InChIKey")
    x$synonyms <- fun_format(synonyms, "Synonym")
    return(x)
  })

setMethod("filter", signature = c(x = "job_pubchemr"),
  function(x, log){
    if (length(object(x)) != length(log)) {
      stop("length(object(x)) != length(log)")
    }
    fun <- function(x) x[ log ]
    object(x) <- fun(object(x))
    x$smiles <- fun(x$smiles)
    x$inks <- fun(x$inks)
    x$synonyms <- fun(x$synonyms)
    return(x)
  })

setMethod("map", signature = c(x = "job_pubchemr", ref = "character"),
  function(x, ref, extra = NULL, only = T){
    names <- c(unique(ref), extra)
    names <- tolower(names)
    finds <- lapply(x$synonyms,
      function(x) {
        x <- tolower(x)
        if (only) {
          unique(names[ names %in% x ])
        } else {
          names[ names %in% x ]
        }
      })
    lst <- lst_clear0(finds)
    if (only) {
      unlist(lst)
    } else {
      lst
    }
  })

setMethod("map", signature = c(x = "job_tcmsp", ref = "job_pubchemr"),
  function(x, ref, name_recode = NULL, filter = T, merge = T)
  {
    if (x@step < 2L) {
      stop("x@step < 2L")
    }
    if (!is.null(name_recode)) {
      fun_recode <- function(data) {
        data <- dplyr::mutate(data,
          `Molecule Name` = dplyr::recode(`Molecule Name`, !!!as.list(name_recode),
            .default = `Molecule Name`)
        )
        data
      }
      x@tables$step2$compounds_targets %<>% fun_recode
      x@tables$step2$ingredients %<>% fun_recode
    }
    if (filter) {
      finds <- map(ref, x@tables$step2$ingredients$`Molecule Name`)
      x <- filter(x, tolower(`Molecule Name`) %in% rm.no(finds))
    }
    if (merge) {
      finds <- map(ref, x@tables$step2$ingredients$`Molecule Name`)
      lst <- nl(unname(finds), names(finds))
      x@tables$step2$ingredients <- dplyr::mutate(x@tables$step2$ingredients,
        cid = dplyr::recode(tolower(`Molecule Name`), !!!lst, .default = NA_character_),
        smiles = dplyr::recode(cid, !!!ref$smiles, .default = NA_character_))
    }
    return(x)
  })

try_get_cids.name <- function(name) {
  db <- PubChemR::get_cids(name)
  db <- dplyr::distinct(db, Identifier, .keep_all = T)
  db <- dplyr::mutate(db, CID = ifelse(grpl(CID, "^No"), NA_integer_, as.integer(CID)))
  db
}


