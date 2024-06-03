# ==========================================================================
# workflow of unitmp
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_unitmp <- setClass("job_unitmp", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "unitmp",
    info = c("https://www.unitmp.org/"),
    cite = "[@UnitmpUnifiedDobson2024; @TheHumanTransDobson2015]",
    method = "The UNIfied database of TransMembrane Proteins (UniTmp) was used for transmembrane protein information retrieving",
    tag = "unitmp"
    ))

job_unitmp <- function(mode = c("htp_all", "htp_exists"), db = .prefix("unitmp", "db"))
{
  mode <- match.arg(mode)
  url <- switch(mode, htp_all = "https://htp.unitmp.org/data/HTP/data/d.2.1/sets/htp_all.xml",
    htp_exists = "https://htp.unitmp.org/data/HTP/data/d.2.1/sets/Exists.xml"
  )
  file <- paste0(db, "/", get_filename(url))
  if (!file.exists(file)) {
    message("Download file from URL:", url)
    xmlStr <- RCurl::getURL(url)
    writeLines(xmlStr, file)
  }
  xml <- XML::xmlParse(file)
  attrs <- XML::xpathApply(xml, "//HtpItem", XML::xmlAttrs)
  values <- XML::xpathApply(xml, "//HtpItem/Name", XML::xmlValue)
  if (length(attrs) != length(values)) {
    stop("length(attrs) != length(values)")
  }
  dat <- dplyr::bind_rows(attrs)
  dat$Protein_name <- unlist(values)
  dat <- .set_lab(dat, "UniTmp data of ", mode)
  ## get Gene Name
  message("Try Mapping Gene Names.")
  mapdb <- paste0(db, "/EntryName2gene.rds")
  ldb <- new_db(mapdb, "From")
  ldb <- not(ldb, dat$id)
  if (length(ldb@query)) {
    mapped <- e(UniProt.ws::mapUniProt(
      from = "UniProtKB_AC-ID",
      to = "Gene_Name",
      columns = c("accession", "id"),
      query = ldb@query
    ))
    ldb <- upd(ldb, mapped, ldb@query)
  }
  mapped <- dplyr::filter(ldb@db, From %in% dat$id)
  dat <- map(dat, "id", mapped, "From", "To", col = "Gene_Name")
  .job_unitmp(object = dat)
}

setMethod("step0", signature = c(x = "job_unitmp"),
  function(x){
    step_message("Prepare your data with function `job_unitmp`.")
  })

setMethod("step1", signature = c(x = "job_unitmp"),
  function(x){
    step_message("Quality control (QC).")
    return(x)
  })
