# ==========================================================================
# workflow of methy
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_methy <- setClass("job_methy", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "methy",
    info = c("https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html",
      "https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html",
      "https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html#43_DataTrack"
      ),
    cite = "[@MissmethylAnPhipso2016; @MinfiAFlexibAryee2014; @VisualizingGenHahne2016; @CallingDifferePeters2021]",
    method = "R package `missMethyl`, `minfi`, `Gviz`, `DMRcate` were used for analysing methylation array data",
    tag = "methyl:array"
    ))

job_methy <- function()
{
  .job_methy()
}

setMethod("step0", signature = c(x = "job_methy"),
  function(x){
    step_message("Prepare your data with function `job_methy`.")
  })

setMethod("step1", signature = c(x = "job_methy"),
  function(x){
    step_message("")
    return(x)
  })
