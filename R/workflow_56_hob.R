# ==========================================================================
# workflow of hob
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_hob <- setClass("job_hob", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://github.com/whymin/HOB"),
    cite = "[@HobpreAccuratWeiM2022]",
    method = "Python tool of `HOB` was used for prediction of human oral bioavailability",
    tag = "hob",
    analysis = "HOB 人类口服利用度预测"
    ))

job_hob <- function(smiles)
{
  if (is.data.frame(smiles)) {
    message("The input is data.frame. Use first column as names, next column as smiles.")
    smiles <- nl(smiles[[1]], smiles[[2]], F)
  }
  smiles <- smiles[ !duplicated(smiles) ]
  .job_hob(object = smiles)
}

setGeneric("asjob_hob", group = list("asjob_series"),
  function(x, ...) standardGeneric("asjob_hob"))

setMethod("asjob_hob", signature = c(x = "job_pubchemr"),
  function(x){
    .job_hob(object = unlist(x$smiles))
  })

setMethod("step0", signature = c(x = "job_hob"),
  function(x){
    step_message("Prepare your data with function `job_hob`.")
  })

setMethod("step1", signature = c(x = "job_hob"),
  function(x, cutoff = c("20", "50"), wd = timeName("hob"),
    command = paste(pg("hobPython"), pg("hobPredict")),
    model = pg("hobModel"), extra = pg("hobExtra"), dir = NULL)
  {
    step_message("Run HOB for Prediction.")
    if (is.null(dir)) {
      if (!file.exists(strx(command, "[^ ]+$"))) {
        stop('file.exists(strx(command, "[^ ]+$")) == F')
      }
      if (!file.exists(model)) {
        stop("file.exists(model) == F")
      }
      x$wd <- wd
      dir.create(wd, F)
      file.copy(extra, wd)
      writeLines(object(x), paste0(wd, "/smiles.txt"))
      cutoff <- match.arg(cutoff)
      cdRun(command, " ", model, " smiles.txt ", cutoff, path = x$wd)
    } else {
      if (dir.exists(dir))
        wd <- dir
    }
    file_res <- paste0(wd, "/pred_result.csv")
    t.hob <- ftibble(file_res)
    if (!is.null(names(object(x)))) {
      if (.is.cid(names(object(x)))) {
        names <- names(object(x))
        lab.x <- "Compound CID"
      } else {
        names <- stringr::str_trunc(names(object(x)), 30)
        lab.x <- "Compound"
      }
    } else {
      names <- paste0("Compound ", seq_along(object(x)))
      lab.x <- ""
    }
    t.hob <- dplyr::mutate(t.hob, name = !!names, isOK = ifelse(prediction, T, F))
    if (F) {
      p.hob <- ggplot(t.hob) +
        geom_col(aes(x = reorder(name, `probability(+)`), y = `probability(+)`, fill = `probability(+)`),
          width = .7) +
        coord_flip() +
        geom_hline(yintercept = 50, linetype = "dashed", size = 0.7, color = "grey") +
        labs(x = lab.x, y = "Probability") +
        ggtitle("HOB (20%) prediction") +
        guides(fill = "none") +
        rstyle("theme") +
        theme(plot.title = element_text(size = 10))
      if (nrow(t.hob) > 100) {
        p.hob <- p.hob + facet_wrap(~ isOK, scales = "free")
        p.hob <- wrap(p.hob, 7, 2 + .07 * nrow(t.hob))
      } else {
        p.hob <- wrap(p.hob, 7, 2 + .15 * nrow(t.hob))
      }
    } else {
      p.hob <- new_pie(t.hob$isOK, title = "HOB (20%) Prediction")
    }
    x <- methodAdd(x,
      "以 Python 工具 `HOB` {cite_show('HobpreAccuratWeiM2022')} 预测化合物人类口服利用度 ({cutoff}%)。")
    x <- methodAdd(x, "化合物以 PubChem CID 获取结构式 (SMILES) 随后通过 `HOB` 程序预测口服利用度 {cutoff}% 是否达标。")
    x <- snapAdd(x, "`HOB` 预测结果，所有用于预测的化合物 (含有结构式信息的) {nrow(t.hob)} 个，达到口服利用度标准的有 {lenth(which(t.hob$isOK))} 个。")
    p.hob <- .set_lab(p.hob, sig(x), "HOB 20 prediction")
    t.hob <- .set_lab(t.hob, sig(x), "data of HOB 20 prediction")
    x@tables[[ 1 ]] <- namel(t.hob)
    x@plots[[ 1 ]] <- namel(p.hob)
    return(x)
  })

setMethod("map", signature = c(x = "job_pubchemr", ref = "job_hob"),
  function(x, ref){
    object(x)[as.logical(res(ref)$prediction)]
  })

setMethod("res", signature = c(x = "job_hob"),
  function(x){
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    x@tables$step1$t.hob
  })

.is.cid <- function(x) {
  all(grpl(x, "^[0-9]+$"))
}
