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
    info = paste0("Tutorial: https://gatk.broadinstitute.org/hc/en-us/sections/360007226631-Tutorials",
      "\nhttps://github.com/biod/sambamba",
      "\nhttps://github.com/ExaScience/elprep",
      "\nhttps://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz",
      "\nhttps://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0"
    )
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
  function(x, ref = "../hg38.fa",
    knownSites = c("../resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf",
      "../resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"),
    vqsr_resource = knownSites,
    gatk = "~/operation/gatk4/gatk",
    picard = "~/operation/picard.jar",
    java = "/usr/lib/jvm/java-18-openjdk-amd64/bin/java",
    tmpdir = "~/disk_sdb1") 
  {
    step_message("Prepare reference data (human) for mapping.")
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
    x@params$knownSites <- knownSites
    x@params$java <- java
    x@params$ref <- normalizePath(ref)
    x@params$picard <- picard
    x@params$gatk <- gatk
    x@params$tmpdir <- tmpdir
    return(x)
  })

setMethod("step2", signature = c(x = "job_gatk"),
  function(x, workers = 9, only.bwa = T, use.sambamba = T){
    step_message("Alignment to reference gnome.
      Then remote duplicated.
      "
    )
    WorkflowGatk.bwa_mem(x)
    if (!only.bwa) {
      if (!use.sambamba) {
        WorkflowGatk.samtools_sort(x)
        WorkflowGatk.picard_MarkDuplicates(x)
        WorkflowGatk.picard_BuildBamIndex(x)
      } else {
        WorkflowGatk.sambamba_markdup(x, workers)
      }
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_gatk"),
  function(x, elprep = "conda run -n base elprep",
    bcftools = "conda run -n base bcftools",
    batch = F, mem = 28, workers = 9)
  {
    step_message("Use elprep for do almost everything.")
    if (!is.null(elprep)) {
      x@params$elprep <- elprep
      x <- WorkflowGatk.elprep_prepare(x)
      if (batch) {
        x <- WorkflowGatk.elprep_bqsr_haplotypecaller_batch(x, workers)
      } else {
        x@params$bcftools <- bcftools
        WorkflowGatk.elprep_bqsr_haplotypecaller(x, workers)
        x <- WorkflowGatk.bcftools_merge(x, workers)
      }
      x@params$use.elprep <- T
    } else {
      stop("...")
    }
    return(x)
  })

setMethod("step4", signature = c(x = "job_gatk"),
  function(x){
    step_message("GATK VQSR")
    if (is.null(x@params$use.elprep)) {
      if (is.null(x@params$vcf)) {
        x@params$vcf <- "all.vcf.gz"
      }
      x <- WorkflowGatk.gatk_prepare_refVcf(x)
      WorkflowGatk.gatk_VariantRecalibrator_snps(x)
      WorkflowGatk.gatk_VariantRecalibrator_indels(x)
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_gatk"),
  function(x, annovar = "~/operation/annovar/annotate_variation.pl", ref = "hg38",
    ref_used = c("refGene", "cytoBand", "exac03", "avsnp147", "dbnsfp30a")[1],
    ref_operation = c("gx", "r", "f", "f", "f")[1])
  {
    step_message("Use annovar for annotation of vcf.")
    path_annovar <- get_path(annovar)
    ref_path <- paste0(path_annovar, "/humandb")
    ref_pattern <- paste0("^", ref)
    if (length(list.files(ref_path, ref_pattern)) < length(ref_used)) {
      pbapply::pblapply(ref_used,
        function(name) {
          if (!file.exists(paste0(ref_path, "/", ref, "_", name, ".txt"))) {
            cdRun(annovar, " -buildver ", ref, " -downdb",
              ifelse(name == "cytoBand", NULL, " -webfrom annovar "),
              " ", name, " humandb",
              path = path_annovar)
          }
        }
      )
    }
    if (is.null(x@params$vcf_for_anno)) {
      x@params$vcf_for_anno <- "all.vcf.gz"
    }
    vcf <- x@params$vcf_for_anno
    output <- gs(vcf, "\\.vcf[^/]*$", ".input")
    if (!is_workflow_object_exists(output)) {
      cdRun(paste0(path_annovar, "/", "convert2annovar.pl"),
        " -format vcf4 -allsample -withfreq ", vcf,
        " > ", output,
        path = x@params$wd
      )
    }
    inputfile <- output
    cdRun(paste0(path_annovar, "/", "table_annovar.pl"),
      " ", inputfile, " ",
      ref_path, " ", " -buildver ", ref,
      " --outfile res_", get_realname(inputfile),
      " -remove -protocol ", paste0(ref_used, collapse = ","),
      " -operation ", paste0(ref_operation, collapse = ","),
      path = x@params$wd)
    x@params$vcf <- vcf
    return(x)
  })

setMethod("is_workflow_object_exists", signature = c(object = "character"),
  function(object, x, path = NULL){
    if (missing(x)) {
      x <- get("x", envir = parent.frame(2))
    }
    if (is.null(path)) {
      path <- x@params$wd
    }
    file <- paste0(path, "/", object)
    if (file.exists(file)) {
      T
    } else {
      F
    }
  })

WorkflowGatk.bwa_mem <- function(x) {
  cli::cli_alert_info("bwa mem")
  use.minimap2 <- F
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
}

WorkflowGatk.samtools_sort <- function(x) {
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
}

WorkflowGatk.sambamba_markdup <- function(x, workers) {
  cli::cli_alert_info("sambamba markdup")
  pbapply::pbapply(object(x), 1,
    function(vec) {
      vec <- as.list(vec)
      output <- paste0(vec$SampleName, ".markdup.bam")
      if (!is_workflow_object_exists(output)) {
        cdRun("sambamba markdup", " -r ", " -p ",
          " -t ", workers,
          " --overflow-list-size 600000",
          " ", vec$SampleName, ".bam",
          " ", output,
          path = x@params$wd
        )
      }
    })
}

WorkflowGatk.picard_MarkDuplicates <- function(x) {
  cli::cli_alert_info("picard MarkDuplicates")
  pbapply::pbapply(object(x), 1,
    function(vec) {
      vec <- as.list(vec)
      output <- paste0(vec$SampleName, ".markdup.bam")
      if (!is_workflow_object_exists(output)) {
        cdRun(x@params$java, " -jar ",
          " ", x@params$picard, " MarkDuplicates",
          " I=", vec$SampleName, ".sorted.bam",
          " O=", output,
          " M=", vec$SampleName, ".markdup.txt",
          " REMOVE_DUPLICATES=true",
          path = x@params$wd
        )
      }
    })
}

WorkflowGatk.picard_BuildBamIndex <- function(x) {
  cli::cli_alert_info("picard BuildBamIndex")
  pbapply::pbapply(object(x), 1,
    function(vec) {
      vec <- as.list(vec)
      output <- paste0(vec$SampleName, ".bai")
      if (!is_workflow_object_exists(output)) {
        cdRun(x@params$java, " -jar ",
          " ", x@params$picard, " BuildBamIndex",
          " -I=", vec$SampleName, ".markdup.bam",
          " -O=", output,
          path = x@params$wd
        )
      }
    })
}

WorkflowGatk.elprep_prepare <- function(x) {
  cli::cli_alert_info("elprep fasta-to-eflasta")
  if (is.null(x@params$ref_elfasta)) {
    ref <- get_filename(x@params$ref)
    ref_elfasta <- paste0(get_realname(ref), ".elfasta")
    path <- get_path(x@params$ref)
    if (!is_workflow_object_exists(ref_elfasta, path = path)) {
      cdRun(x@params$elprep, " fasta-to-elfasta", " ", ref, " ", ref_elfasta,
        path = path)
    }
    x@params$ref_elfasta <- normalizePath(paste0(path, "/", ref_elfasta))
  }
  cli::cli_alert_info("elprep vcf-to-elsites")
  if (is.null(x@params$elsites)) {
    elsites <- pbapply::pblapply(x@params$knownSites,
      function(file) {
        path <- get_path(file)
        filename <- get_filename(file)
        real <- gs(filename, "\\.[a-z]+$", "")
        elsites <- paste0(real, ".elsites")
        if (!is_workflow_object_exists(elsites, path = path)) {
          cdRun(x@params$elprep,
            " vcf-to-elsites ", filename, " ", elsites,
            path = path)
        }
        normalizePath(paste0(path, "/", elsites))
      })
    x@params$elsites <- unlist(elsites)
  }
  return(x)
}

WorkflowGatk.elprep_bqsr_haplotypecaller <- function(x, workers) {
  cli::cli_alert_info("elprep sfm")
  pbapply::pbapply(object(x), 1,
    function(vec) {
      vec <- as.list(vec)
      output <- paste0(vec$SampleName, ".vcf.gz")
      if (!is_workflow_object_exists(output)) {
        cdRun(x@params$elprep,
          " sfm", " ", vec$SampleName, ".bam", " ",
          " ", vec$SampleName, ".elprep.bam",
          " --nr-of-threads ", workers,
          " --tmp-path ", x@params$tmpdir,
          " --mark-duplicates",
          " --mark-optical-duplicates ", vec$SampleName, ".metrics",
          " --sorting-order coordinate",
          " --bqsr ", vec$SampleName, ".recal",
          " --known-sites ", paste0(x@params$elsites, collapse = ","),
          " --reference ", x@params$ref_elfasta,
          " --haplotypecaller ", output,
          path = x@params$wd
        )
      }
    })
}

WorkflowGatk.elprep_bqsr_haplotypecaller_batch <- function(x, workers) {
  output <- paste0("all.vcf.gz")
  if (!is_workflow_object_exists("elprep")) {
    cdRun("mkdir elprep", path = x@params$wd)
    cdRun("mv ", paste0(object(x)$SampleName, ".bam", collapse = " "),
      " -t elprep ", path = x@params$wd)
  }
  cli::cli_alert_info("elprep sfm")
  if (!is_workflow_object_exists(output)) {
    cdRun(x@params$elprep,
      " sfm", " elprep ",
      " all.elprep.bam",
      " --nr-of-threads ", workers,
      " --tmp-path ", x@params$tmpdir,
      " --mark-duplicates",
      " --mark-optical-duplicates ", "all.metrics",
      " --sorting-order coordinate",
      " --bqsr all.recal",
      " --known-sites ", paste0(x@params$elsites, collapse = ","),
      " --reference ", x@params$ref_elfasta,
      " --haplotypecaller ", output,
      path = x@params$wd
    )
  }
  x@params$vcf <- "all.vcf.gz"
  return(x)
}

WorkflowGatk.gatk_prepare_refVcf <- function(x) {
  if (is.null(x@params$vqsr_resource)) {
    x@params$vqsr_resource <- x@params$knownSites
  }
  lapply(x@params$vqsr_resource,
    function(file) {
      path <- get_path(file)
      filename <- get_filename(file)
      if (!is_workflow_object_exists(paste0(filename, ".idx"), path = path)) {
        cdRun(x@params$gatk, " IndexFeatureFile",
          " -I ", filename, path = path)
      }
    })
  x@params$vqsr_resource %<>% normalizePath()
  return(x)
}

create_gatk_vcf_index <- function(vcf, x) {
  if (!is_workflow_object_exists(paste0(vcf, ".tbi"))) {
    cdRun(x@params$gatk, " IndexFeatureFile",
      " -I ", vcf, path = x@params$wd)
  }
}

WorkflowGatk.gatk_VariantRecalibrator_snps <- function(x) {
  resources <- x@params$vqsr_resource
  snp_res <- resources[ grepl("\\.snps\\.", resources) ]
  snp_res <- paste0(
    " -resource:", get_realname(snp_res), ",",
    "known=false,training=true,truth=true,prior=",
    13 - 1:length(snp_res),
    " ", snp_res, " "
  )
  create_gatk_vcf_index(x@params$vcf, x)
  if (!is_workflow_object_exists("tmp_snp.tranches")) {
    E(cdRun(x@params$gatk, " VariantRecalibrator",
        " -R ", x@params$ref,
        " -V ", x@params$vcf,
        paste0(snp_res, collapse = " "),
        " -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum",
        " -mode SNP",
        " -O tmp_snp.recal ",
        " --tranches-file tmp_snp.tranches ",
        " --rscript-file tmp_snp.plots.R",
        path = x@params$wd
        ))
  }
  if (!is_workflow_object_exists("all_snps.VQSR.vcf.gz")) {
    E(cdRun(x@params$gatk, " ApplyVQSR",
        " -R ", x@params$ref,
        " -V ", x@params$vcf,
        " --ts-filter-level 99.0 --tranches-file tmp_snp.tranches ",
        " --recal-file tmp_snp.recal ",
        " -mode SNP",
        " -O all_snps.VQSR.vcf.gz",
        path = x@params$wd
        ))
  }
}

WorkflowGatk.gatk_VariantRecalibrator_indels <- function(x) {
  resources <- x@params$vqsr_resource
  indel_res <- resources[ grepl("\\.indels\\.", resources) ]
  indel_res <- paste0(
    " -resource:", get_realname(indel_res), ",",
    "known=true,training=true,truth=true,prior=",
    13 - 1:length(indel_res),
    " ", indel_res, " "
  )
  vcf <- "all_snps.VQSR.vcf.gz"
  create_gatk_vcf_index(vcf, x)
  if (!is_workflow_object_exists("tmp_indel.tranches")) {
    E(cdRun(x@params$gatk, " VariantRecalibrator",
        " -R ", x@params$ref,
        " -V ", vcf,
        paste0(indel_res, collapse = " "),
        " -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum",
        " -mode INDEL", " --max-gaussians 6",
        " -O tmp_indel.recal ",
        " --tranches-file tmp_indel.tranches ",
        " --rscript-file tmp_indel.plots.R",
        path = x@params$wd
        ))
  }
  if (!is_workflow_object_exists("all_snps_indels.VQSR.vcf.gz")) {
    E(cdRun(x@params$gatk, " ApplyVQSR",
        " -R ", x@params$ref,
        " -V ", vcf,
        " --ts-filter-level 99.0 --tranches-file tmp_indel.tranches ",
        " --recal-file tmp_indel.recal ",
        " -mode INDEL",
        " -O all_snps_indels.VQSR.vcf.gz",
        path = x@params$wd
        ))
  }
}

WorkflowGatk.bcftools_merge <- function(x) {
  pbapply::pblapply(object(x),
    function(vec) {
      vec <- as.list(vec)
      output <- paste0(vec$SampleName, ".vcf.gz.csi")
      if (!is_workflow_object_exists(output)) { 
        cdRun(x@params$bcftools,
          " index ", vec$SampleName, ".vcf.gz",
          path = x@params$wd)
      }
    })
  vcfs <- paste0(paste0(object(x)$SampleName, ".vcf.gz"), collapse = " ")
  output <- "all.vcf.gz"
  if (!is_workflow_object_exists(output)) {
    cdRun(x@params$bcftools, " merge ", vcfs, " -O z -o ", output)
  }
  x@params$vcf <- output
  return(x)
}
