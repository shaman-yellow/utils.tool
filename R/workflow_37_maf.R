# ==========================================================================
# workflow of maf
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_maf <- setClass("job_maf", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://bioconductor.org/packages/release/bioc/html/maftools.html"),
    cite = "[@MaftoolsEfficMayako2018]",
    method = "R package 'maftools' used for analyzing and visualizing Mutation Annotation Format (MAF) files",
    tag = "maf",
    analysis = "Maftools 变异注释结果分析"
    ))

setGeneric("asjob_maf", 
  function(x, ...) standardGeneric("asjob_maf"))

setMethod("asjob_maf", signature = c(x = "job_tcga"),
  function(x, use = "follow_up", keep_consensus = T){
    n.all_raw <- colSum(object(x)$Tumor_Sample_Barcode)
    object <- e(maftools::read.maf(object(x), isTCGA = T))
    n.all <- colSum(maftools::getClinicalData(object)$Tumor_Sample_Barcode)
    if (!is.null(x$queries$clinical)) {
      clinical <- e(TCGAbiolinks::GDCprepare_clinic(x$queries$clinical, use))
      if (use == "follow_up" & keep_consensus) {
        ## check duplicated
        data <- dplyr::distinct(clinical, bcr_patient_barcode, vital_status)
        stats <- table(data$bcr_patient_barcode)
        clinical <- dplyr::filter(clinical, bcr_patient_barcode %in% names(stats[stats == 1]))
        clinical <- dplyr::distinct(clinical, bcr_patient_barcode, .keep_all = T)
        clinical <- dplyr::filter(clinical, vital_status != "")
        object.bar <- e(maftools::getClinicalData(object)$Tumor_Sample_Barcode)
        clinical <- dplyr::filter(clinical, bcr_patient_barcode %in% object.bar)
        notInCli <- object.bar[ !object.bar %in% clinical$bcr_patient_barcode ]
        object <- e(maftools::filterMaf(object, tsb = notInCli))
        message("In clinical: ", n.cli <- colSum(clinical$bcr_patient_barcode), ".\n",
          "In maf object: ", n.sam <- colSum(maftools::getClinicalData(object)$Tumor_Sample_Barcode))
        format.cli <- dplyr::rename(
          dplyr::mutate(clinical, vital_status = dplyr::recode(vital_status, Alive = 1, Dead = 0)),
          Tumor_Sample_Barcode = bcr_patient_barcode)
      }
    }
    x <- .job_maf(object = object)
    x$clinical <- as_tibble(format.cli)
    x$info <- namel(n.all_raw, n.all, n.sam, n.cli)
    x$info.lich <- new_lich(list(`All dowloaded data number` = n.all_raw,
        `TCGAbiolinks filtered` = n.all,
        `With clinical data (vital_status)` = n.cli,
        `Finaly used sample number` = n.sam
        ))
    return(x)
  })

# mafSurvival(object(mf), "TP53",
#   clinicalData = clinicalData,
#   time = "days_to_last_followup", Status = "vital_status"
# )

setMethod("step0", signature = c(x = "job_maf"),
  function(x){
    step_message("Prepare your data with function `asjob_maf`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_maf"),
  function(x, genes = NULL, top = 30){
    step_message("Visualize the overview of data.")
    e(maftools::plotmafSummary(object(x), addStat = 'median', titvRaw = T))
    p.summary <- wrap(recordPlot())
    p.summary <- .set_lab(p.summary, sig(x), "summary of mutation")
    e(maftools::titv(object(x)))
    p.snp_class <- wrap(recordPlot())
    e(maftools::oncoplot(object(x), top = top, genes = genes))
    p.oncoplot <- wrap(recordPlot())
    p.oncoplot <- .set_lab(p.oncoplot, sig(x), "oncoplot of top genes")
    x@plots[[ 1 ]] <- namel(p.summary, p.snp_class, p.oncoplot)
    return(x)
  })


