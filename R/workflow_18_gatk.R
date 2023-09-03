# ==========================================================================
# workflow of gatk
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gatk <- setClass("job_gatk", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://gatk.broadinstitute.org/hc/en-us/sections/360007226631-Tutorials")
    ))

setGeneric("asjob_gatk", 
  function(x, ...) standardGeneric("asjob_gatk"))

setMethod("asjob_gatk", signature = c(x = "job_fastp"),
  function(x, wd){
    if (any(duplicated(x@params$metadata$SampleName)))
      stop("any(duplicated(x@params$metadata$SampleName)) == T")
    x <- .job_gatk(object = x@params$metadata)
    x@params$wd <- wd
    if (!file.exists(wd))
      dir.create(wd)
    x
  })

setMethod("step0", signature = c(x = "job_gatk"),
  function(x){
    step_message("Prepare your data with function `asjob_gatk`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_gatk"),
  function(x, ref = "../hg38.fa", picard = "~/operation/picard.jar",
    java = "/usr/lib/jvm/java-18-openjdk-amd64/bin/java")
  {
    step_message("Quality control (QC).
      This do:
      "
    )
    path <- get_path(ref)
    if (!file.exists(ref)) {
      cdRun("wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz",
        path = path)
      cdRun("gunzip hg38.fa.gz", path = path)
      cdRun("bwa index -a bwtsw hg38.fa", path = path)
      cdRun("samtools faidx hg38.fa", path = path)
      cdRun(java, " -jar ", picard,
        " CreateSequenceDictionary R=hg38.fa O=hg38.dict",
        path = path)
    }
    x@params$java <- java
    x@params$ref <- normalizePath(ref)
    x@params$picard <- picard
    return(x)
  })

setMethod("step2", signature = c(x = "job_gatk"),
  function(x, workers = 9, mem = 28, use.minimap2 = F){
    step_message("Alignment to reference gnome.
      Then sorted the bam file.
      ...
      "
    )
    cli::cli_alert_info("bwa mem")
    pbapply::pbapply(object(x), 1,
        function(vec) {
          vec <- as.list(vec)
          output <- paste0(vec$SampleName, ".bam")
          if (!is_workflow_object_exists(output)) {
            if (!use.minimap2) {
              cdRun("bwa mem",
                " -t ", workers, " -M -Y ",
                " -R '@RG\\tID:", vec$SampleName, "\\tSM:", vec$SampleName, "\\tLB:WES\\tPL:Illumina'",
                " ", x@params$ref, 
                " ", vec[[ "forward-absolute-filepath" ]],
                " ", vec[[ "reverse-absolute-filepath" ]],
                " | samtools view -Sb - > ", output,
                path = x@params$wd
              )
            } else {
              stop("emm... the code of using 'minimap2' is in the future.")
            }
          }
        })
    cli::cli_alert_info("samtools sort")
    pbapply::pbapply(object(x), 1,
      function(vec) {
        vec <- as.list(vec)
        output <- paste0(vec$SampleName, ".sorted.bam")
        if (!is_workflow_object_exists(output)) {
          cdRun("samtools sort ",
            " -@ ", workers,
            " -o ", output,
            " ", vec$SampleName, ".bam",
            path = x@params$wd
          )
        }
      })
    cli::cli_alert_info("picard MarkDuplicates")
    pbapply::pbapply(object(x), 1,
      function(vec) {
        vec <- as.list(vec)
        output <- paste0(vec$SampleName, ".sorted.markdup.bam")
        if (!is_workflow_object_exists(output)) {
          cdRun(x@params$java, " -jar ",
            " ", x@params$picard, " MarkDuplicates",
            " -Xmx", mem, "g ",
            " -I=", vec$SampleName, ".sorted.bam",
            " -M=", vec$SampleName, ".sorted.markdup.txt",
            " -O=", output,
            " --REMOVE_DUPLICATES=true",
            path = x@params$wd
          )
        }
      })
    cli::cli_alert_info("picard BuildBamIndex")
    pbapply::pbapply(object(x), 1,
      function(vec) {
        vec <- as.list(vec)
        output <- paste0(vec$SampleName, ".bai")
        if (!is_workflow_object_exists(output)) {
          cdRun(x@params$java, " -jar ",
            " ", x@params$picard, " BuildBamIndex",
            " -Xmx", mem, "g ",
            " -I=", vec$SampleName, ".sorted.markdup.bam",
            " -O=", output,
            path = x@params$wd
          )
        }
      })
    return(x)
  })

setMethod("step3", signature = c(x = "job_gatk"),
  function(x, workers = 9){
    step_message("...")
    return(x)
  })

setMethod("is_workflow_object_exists", signature = c(object = "character"),
  function(object, x){
    if (missing(x)) {
      x <- get("x", envir = parent.frame(2))
    }
    path <- x@params$wd
    file <- paste0(path, "/", object)
    if (file.exists(file)) {
      T
    } else {
      F
    }
  })
