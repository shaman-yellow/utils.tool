# ==========================================================================
# workflow of prr
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_prr <- setClass("job_prr",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "prr",
    info = c("https://github.com/paulgeeleher/pRRophetic2; https://osf.io/5xvsg; https://github.com/paulgeeleher/pRRophetic2/tree/master/pRRophetic/vignettes; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4167990/"),
    cite = "[@PrropheticAnGeeleh2014]",
    method = "R Package `pRRophetic` was used for Prediction of Clinical Chemotherapeutic Response"
    ))

job_prr <- function(data, drug = prr_drug())
{
  if (is(data, "tbl_df")) {
    if (is.character(data[[1]])) {
      message("Convert `data` to 'matrix', use first column as rownames.")
      rownames <- data[[1]]
      data <- as.matrix(data[, -1])
      rownames(data) <- gname(rownames)
    }
  } else if (!is(data, "matrix")) {
    stop('!is(data, "matrix")')
  }
  setupPRR()
  drug <- match.arg(drug)
  x <- .job_prr(object = data)
  x$drug <- drug
  x$p.drugQQ <- wrap(as_grob(expression(e(pRRophetic::pRRopheticQQplot(drug))), environment()))
  x$p.drugQQ <- .set_lab(x$p.drugQQ, sig(x), "QQ plot for distribution of the transformed IC50 data")
  return(x)
}

setGeneric("asjob_prr", 
  function(x, ...) standardGeneric("asjob_prr"))

setMethod("asjob_prr", signature = c(x = "job_tcga"),
  function(x, drug = prr_drug())
  {
    drug <- match.arg(drug)
    if (missing(drug)) {
      stop("missing(drug)")
    }
    if (x@step != 2) {
      stop("x@step != 2")
    }
    if (!all(c("clinical", "RNA", "protein") %in% names(object(x)))) {
      stop("RNA, protein data should covered in `x` (job_tcga)")
    }
    message("get RNA data")
    x <- suppressMessages(step3(x, query = "RNA"))
    message("for normalization")
    x <- suppressMessages(asjob_limma(x))
    x <- suppressMessages(step1(x, norm_vis = F))
    ## forma
    rownames <- gname(object(x)$genes$gene_name)
    useWhich <- !is.na(rownames) & !duplicated(rownames)
    data <- object(x)$E[useWhich, ]
    rownames(data) <- rownames[ useWhich ]
    x <- job_prr(data, drug)
    return(x)
  })

setMethod("asjob_prr", signature = c(x = "job_limma"),
  function(x, drug = prr_drug())
  {
    drug <- match.arg(drug)
    if (is.null(x$isTcga)) {
      stop("is.null(x$isTcga)")
    }
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    rownames <- gname(object(x)$genes$gene_name)
    useWhich <- !is.na(rownames) & !duplicated(rownames)
    data <- object(x)$E[useWhich, ]
    rownames(data) <- rownames[ useWhich ]
    x <- job_prr(data, drug)
    return(x)
  })

setMethod("step0", signature = c(x = "job_prr"),
  function(x){
    step_message("Prepare your data with function `job_prr`.")
  })


# Sensitivity, Specificity, PPV, and NPV for Predictive Biomarkers
# PMID: 26109105
# Non-response (NR) tissues with Response (R) tissues
#
# The positive predictive value (PPV) of the test can be defined as the
# probability that the survival of a marker-positive patient will be longer if
# the patient receives treatment T than if the patient receives control treatment
# C
#
# The negative predictive value (NPV) is the probability that a
# biomarker-negative patient will not have longer survival on T rather than C. 
# 

setMethod("step1", signature = c(x = "job_prr"),
  function(x){
    step_message("Cross validation on a training set to estimate prediction accuracy.")
    x$cvOut <- cvOut <- e(pRRophetic::pRRopheticCV(x$drug, cvFold = 5, testExprData = object(x)))
    x$summary_cvOut <- summary(cvOut)
    p.cvOut <- wrap(as_grob(expression(plot(cvOut)), environment()))
    p.cvOut <- .set_lab(p.cvOut, sig(x), "estimate prediction accuracy")
    x@plots[[ 1 ]] <- namel(p.cvOut)
    return(x)
  })

setMethod("step2", signature = c(x = "job_prr"),
  function(x, k = 3){
    step_message("Predicted the IC50, and perform k-means clustering.")
    predictedPtype <- e(pRRophetic::pRRopheticPredict(object(x), x$drug, selection = 1))
    group <- kmeans(predictedPtype, k)
    t.predict <- tibble::tibble(
      sample = names(predictedPtype),
      sensitivity = unname(predictedPtype),
      kmeans_group = group$cluster
    )
    t.predict <- .set_lab(t.predict, sig(x), "predicted drug sensitivity")
    attr(t.predict, "lich") <- new_lich(
      list("k-means clustering" = paste0("Centers = ", k))
    )
    x@tables[[ 2 ]] <- namel(t.predict)
    return(x)
  })


setMethod("map", signature = c(x = "job_tcga", ref = "job_prr"),
  function(x, ref, use = "protein")
  {
    use <- match.arg(use)
    if (x@step != 2) {
      stop("x@step != 2")
    }
    if (ref@step != 2) {
      stop("ref@step != 2")
    }
    meta.prr <- ref@tables$step2$t.predict
    if (use == "protein") {
      if (!any(names(object(x)) == "protein")) {
        stop('!any(names(object(x)) == "protein")')
      }
    }
    if (!any(names(object(x)) == "clinical")) {
      stop('!any(names(object(x)) == "clinical")')
    }
    x <- suppressMessages(step3(x, query = "protein", clinical.info = "drug"))
    data <- object(x)
    metadata <- x$metadata
    namel(data, metadata, meta.prr)
  })

setGeneric("do_limma", 
  function(x, ref, ...) standardGeneric("do_limma"))

setMethod("do_limma", signature = c(x = "job_tcga", ref = "job_prr"),
  function(x, ref, use = list("resistance" = max, "non_resistance" = min))
  {
    lst <- map(x, ref)
    metadata <- lst$meta.prr
    counts <- dplyr::select(lst$data, 1, dplyr::starts_with("TCGA"))
    ## the `genes` is protein info
    genes <- dplyr::select(lst$data, AGID:peptide_target)
    ## format
    fstr <- function(x) substr(x, 1, 12)
    colnames(counts)[ -1 ] %<>% fstr()
    metadata <- dplyr::mutate(metadata, sample = fstr(sample))
    ## metadata preparing
    metadata <- .assign_resistance(metadata, use)
    metadata <- dplyr::relocate(metadata, sample, group)
    ## remove not exists sample of protein data
    metadata <- dplyr::filter(metadata, sample %in% !!colnames(counts)[-1])
    ## the `genes` is protein info
    x <- job_limma_normed(counts, metadata, genes)
    x$extra_metadata <- dplyr::slice(lst$metadata, match(bcr_patient_barcode, metadata$sample))
    return(x)
  })

setMethod("map", signature = c(x = "job_limma", ref = "job_prr"),
  function(x, ref, use = list("resistance" = max, "non_resistance" = min))
  {
    if (x@step != 1L) {
      stop("x@step != 1L")
    }
    metadata <- ref@tables$step2$t.predict
    metadata <- .assign_resistance(metadata, use, filter.na = F)
    metadata <- dplyr::mutate(metadata, group = ifelse(is.na(group), "Others", group))
    object(x)$targets <- map(object(x)$targets, "sample", metadata, "sample", "group",
      col = "predicted_resistance")
    group <- object(x)$targets$predicted_resistance
    print(table(group))
    message("Set group as `predicted_resistance`")
    x$design <- mx(~0 + group)
    object(x)$targets <- dplyr::mutate(object(x)$targets, group = predicted_resistance)
    return(x)
  })

setMethod("map", signature = c(x = "job_limma", ref = "job_limma"),
  function(x, ref, from_tcga = T)
  {
    message("This method mapped the metadata of `ref` into `x`")
    if (x@step > 0) {
      stop("x@step > 0")
    }
    meta <- ref$metadata
    if (is.null(meta)) {
      meta <- ref$.metadata
    }
    if (from_tcga) {
      fstr <- function(x) substr(x, 1, 12)
      meta <- dplyr::mutate(meta, sample = fstr(sample))
      if (is.null(colnames(object(x)))) {
        colnames(object(x)) <- object(x)$samples$sample
      }
      colnames(object(x)) %<>% fstr()
      object(x)$samples <- dplyr::mutate(object(x)$samples, sample = fstr(sample))
    }
    ## subset
    object(x) <- object(x)[, match(meta$sample, colnames(object(x))) ]
    object(x)$samples <- dplyr::mutate(object(x)$samples,
      group = meta$group[ match(sample, meta$sample) ])
    message("After match, identical: ",
      as.character(identical(object(x)$samples$sample, meta$sample)))
    return(x)
  })

setupPRR <- function() {
  require("pRRophetic")
  replaceFunInPackage("calcPhenotype", .calcPhenotype, "pRRophetic")
  replaceFunInPackage("summarizeGenesByMean", .summarizeGenesByMean, "pRRophetic")
}

.calcPhenotype <- function (trainingExprData, trainingPtype, testExprData, batchCorrect = "eb", 
  powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2, 
  minNumSamples = 10, selection = -1, printOutput = TRUE, removeLowVaringGenesFrom = "homogenizeData") 
{
  if (!is(testExprData, "matrix"))
    stop("ERROR: \"testExprData\" must be a matrix.")
  if (!is(trainingExprData, "matrix"))
    stop("ERROR: \"trainingExprData\" must be a matrix.")
  if (!is(trainingPtype, "numeric"))
    stop("ERROR: \"trainingPtype\" must be a numeric vector.")
  if (ncol(trainingExprData) != length(trainingPtype)) 
    stop("The training phenotype must be of the same length as the number of columns of the training expressin matrix.")
  if ((ncol(trainingExprData) < minNumSamples) || (ncol(testExprData) < 
      minNumSamples)) {
    stop(paste("There are less than", minNumSamples, "samples in your test or training set. It is strongly recommended that you use larger numbers of samples in order to (a) correct for batch effects and (b) fit a reliable model. To supress this message, change the \"minNumSamples\" parameter to this function."))
  }
  homData <- homogenizeData(testExprData, trainingExprData, 
    batchCorrect = batchCorrect, selection = selection, printOutput = printOutput)
  if (!(removeLowVaringGenesFrom %in% c("homogenizeData", "rawData"))) {
    stop("\"removeLowVaringGenesFrom\" must be one of \"homogenizeData\", \"rawData\"")
  }
  keepRows <- seq(1:nrow(homData$train))
  if (removeLowVaryingGenes > 0 && removeLowVaryingGenes < 
    1) {
    if (removeLowVaringGenesFrom == "homogenizeData") {
      keepRows <- pRRophetic:::doVariableSelection(cbind(homData$test, 
          homData$train), removeLowVaryingGenes = removeLowVaryingGenes)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) 
        cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."))
    }
    else if (removeLowVaringGenesFrom == "rawData") {
      evaluabeGenes <- rownames(homData$test)
      keepRowsTrain <- pRRophetic:::doVariableSelection(trainingExprData[evaluabeGenes, 
        ], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRowsTest <- pRRophetic:::doVariableSelection(testExprData[evaluabeGenes, 
        ], removeLowVaryingGenes = removeLowVaryingGenes)
      keepRows <- intersect(keepRowsTrain, keepRowsTest)
      numberGenesRemoved <- nrow(homData$test) - length(keepRows)
      if (printOutput) 
        cat(paste("\n", numberGenesRemoved, "low variabilty genes filtered."))
    }
  }
  offset = 0
  if (powerTransformPhenotype) {
    if (min(trainingPtype) < 0) {
      offset <- -min(trainingPtype) + 1
      trainingPtype <- trainingPtype + offset
    }
    transForm <- car::powerTransform(trainingPtype)[[6]]
    trainingPtype <- trainingPtype^transForm
  }
  if (printOutput) 
    cat("\nFitting Ridge Regression model... ")
  trainFrame <- data.frame(Resp = trainingPtype, t(homData$train[keepRows, 
      ]))
  rrModel <- ridge::linearRidge(Resp ~ ., data = trainFrame)
  if (printOutput) 
    cat("Done\n\nCalculating predicted phenotype...")
  totBeta <- sum(abs(coef(rrModel)))
  eachBeta <- abs(coef(rrModel))
  eachContribution <- eachBeta/totBeta
  if (is(homData$test, "numeric")) {
    n <- names(homData$test)
    homData$test <- matrix(homData$test, ncol = 1)
    rownames(homData$test) <- n
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata = rbind(testFrame, 
        testFrame))[1]
  }
  else {
    testFrame <- data.frame(t(homData$test[keepRows, ]))
    preds <- predict(rrModel, newdata = testFrame)
  }
  if (powerTransformPhenotype) {
    preds <- preds^(1/transForm)
    preds <- preds - offset
  }
  if (printOutput) 
    cat("Done\n\n")
  return(preds)
}

.summarizeGenesByMean <- function (exprMat) 
{
    geneIds <- rownames(exprMat)
    t <- table(geneIds)
    allNumDups <- unique(t)
    allNumDups <- allNumDups[-which(allNumDups == 1)]
    exprMatUnique <- exprMat[which(geneIds %in% names(t[t == 
        1])), ]
    gnamesUnique <- geneIds[which(geneIds %in% names(t[t == 1]))]
    for (numDups in allNumDups) {
        geneList <- names(which(t == numDups))
        for (i in 1:length(geneList)) {
            exprMatUnique <- rbind(exprMatUnique, colMeans(exprMat[which(geneIds == 
                geneList[i]), ]))
            gnamesUnique <- c(gnamesUnique, geneList[i])
        }
    }
    if (is(exprMatUnique, "numeric")) {
        exprMatUnique <- matrix(exprMatUnique, ncol = 1)
    }
    rownames(exprMatUnique) <- gnamesUnique
    return(exprMatUnique)
}

prr_drug <- function() {
  c("A.443654", "A.770041", "ABT.263", "ABT.888", 
    "AG.014699", "AICAR", "AKT.inhibitor.VIII", "AMG.706", 
    "AP.24534", "AS601245", "ATRA", "AUY922", "Axitinib", 
    "AZ628", "AZD.0530", "AZD.2281", "AZD6244", "AZD6482", 
    "AZD7762", "AZD8055", "BAY.61.3606", "Bexarotene", "BI.2536", 
    "BIBW2992", "Bicalutamide", "BI.D1870", "BIRB.0796", 
    "Bleomycin", "BMS.509744", "BMS.536924", "BMS.708163", 
    "BMS.754807", "Bortezomib", "Bosutinib", "Bryostatin.1", 
    "BX.795", "Camptothecin", "CCT007093", "CCT018159", "CEP.701", 
    "CGP.082996", "CGP.60474", "CHIR.99021", "CI.1040", "Cisplatin", 
    "CMK", "Cyclopamine", "Cytarabine", "Dasatinib", "DMOG", 
    "Docetaxel", "Doxorubicin", "EHT.1864", "Elesclomol", 
    "Embelin", "Epothilone.B", "Erlotinib", "Etoposide", 
    "FH535", "FTI.277", "GDC.0449", "GDC0941", "Gefitinib", 
    "Gemcitabine", "GNF.2", "GSK269962A", "GSK.650394", "GW.441756", 
    "GW843682X", "Imatinib", "IPA.3", "JNJ.26854165", "JNK.9L", 
    "JNK.Inhibitor.VIII", "JW.7.52.1", "KIN001.135", "KU.55933", 
    "Lapatinib", "Lenalidomide", "LFM.A13", "Metformin", 
    "Methotrexate", "MG.132", "Midostaurin", "Mitomycin.C", 
    "MK.2206", "MS.275", "Nilotinib", "NSC.87877", "NU.7441", 
    "Nutlin.3a", "NVP.BEZ235", "NVP.TAE684", "Obatoclax.Mesylate", 
    "OSI.906", "PAC.1", "Paclitaxel", "Parthenolide", "Pazopanib", 
    "PD.0325901", "PD.0332991", "PD.173074", "PF.02341066", 
    "PF.4708671", "PF.562271", "PHA.665752", "PLX4720", "Pyrimethamine", 
    "QS11", "Rapamycin", "RDEA119", "RO.3306", "Roscovitine", 
    "Salubrinal", "SB.216763", "SB590885", "Shikonin", "SL.0101.1", 
    "Sorafenib", "S.Trityl.L.cysteine", "Sunitinib", "Temsirolimus", 
    "Thapsigargin", "Tipifarnib", "TW.37", "Vinblastine", 
    "Vinorelbine", "Vorinostat", "VX.680", "VX.702", "WH.4.023", 
    "WO2009093972", "WZ.1.84", "X17.AAG", "X681640", "XMD8.85", 
    "Z.LLNle.CHO", "ZM.447439")
}

.assign_resistance <- function(metadata, use, filter.na = T)
{
  data <- dplyr::group_by(metadata, kmeans_group)
  data <- dplyr::summarise(data, mean = mean(sensitivity))
  Resis <- data$kmeans_group[ data$mean == use$resistance(data$mean) ]
  Non_resis <- data$kmeans_group[ data$mean == use$non_resistance(data$mean) ]
  metadata <- dplyr::mutate(metadata,
    group = ifelse(kmeans_group == Resis, "Resistance",
      ifelse(kmeans_group == Non_resis, "Non_resistance", NA_character_))
  )
  if (filter.na) {
    metadata <- dplyr::filter(metadata, !is.na(group))
  }
  metadata
}
