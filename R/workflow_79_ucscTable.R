# ==========================================================================
# workflow of ucsc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ucscTable <- setClass("job_ucscTable", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "ucscTable",
    info = c("https://genome.ucsc.edu/cgi-bin/hgTables"),
    cite = "[@RtracklayerAnLawren2009]",
    method = "R package `rtracklayer` used for UCSC data query",
    tag = "ucsc",
    analysis = "rtracklayer UCSC Table 数据获取"
    ))

# DNAseI hypersensitive sites

job_ucscTable <- function(genome)
{
  x <- .job_ucscTable(object = genome)
  query <- e(rtracklayer::ucscTableQuery(genome))
  x$alls <- e(rtracklayer::tableNames(query))
  return(x)
}

setMethod("step0", signature = c(x = "job_ucscTable"),
  function(x){
    step_message("Prepare your data with function `job_ucscTable`.")
  })

setMethod("step1", signature = c(x = "job_ucscTable"),
  function(x){
    step_message("")
    return(x)
  })

get_dnase_data <- function(genome, chr, range,
  table = "wgEncodeRegDnaseClustered")
{
  query <- rtracklayer::ucscTableQuery(genome, table = table,
    rtracklayer::GRangesForUCSCGenome(genome, chr, range)
  )
  rtracklayer::getTable(query)
}

get_cpgIsland_data.ucsc <- function(genome, table = "cpgIslandExt", ...)
{
  get_data.ucsc(genome, table, ...)
}

get_data.ucsc <- function(genome, table, db_dir = .prefix("ucscTable", "db"), addInternalJob = T)
{
  dir.create(db_dir, F)
  file <- file.path(db_dir, paste0(genome, "_", table, ".csv"))
  if (!file.exists(file)) {
    query <- e(rtracklayer::ucscTableQuery(genome, table = table))
    data <- as_tibble(e(rtracklayer::getTable(query)))
    data.table::fwrite(data, file)
  } else {
    data <- ftibble(file)
  }
  if (addInternalJob) {
    .add_internal_job(.job_ucscTable())
  }
  data
}

# example
# 
# list available tables:
# dat <- rtracklayer::ucscTableQuery("hg38")
# alls <- rtracklayer::tableNames(dat)
#
# Download data:
# dat <- rtracklayer::ucscTableQuery("hg38", table = "wgEncodeRegDnaseClustered")
# datSet <- rtracklayer::getTable(dat)
#
#
# rtracklayer::browserSession
# GenomeInfoDb::genome
# rtracklayer::ucscTableQuery
# rtracklayer::getTable
# rtracklayer::tableNames
#

# range <- GRangesForUCSCGenome("mm9", "chr12", IRanges(57795963, 57815592))
# dat <- rtracklayer::ucscTableQuery("hg38", table = "wgEncodeRegDnaseClustered", range = range)
