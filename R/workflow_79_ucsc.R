# ==========================================================================
# workflow of ucsc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ucsc <- setClass("job_ucsc", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "ucsc",
    info = c("https://genome.ucsc.edu/cgi-bin/hgTables"),
    cite = "[@RtracklayerAnLawren2009]",
    method = "R package `rtracklayer` used for UCSC data query"
    ))

job_ucsc <- function()
{
  .job_ucsc()
}

setMethod("step0", signature = c(x = "job_ucsc"),
  function(x){
    step_message("Prepare your data with function `job_ucsc`.")
  })

setMethod("step1", signature = c(x = "job_ucsc"),
  function(x){
    step_message("")
    return(x)
  })

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
