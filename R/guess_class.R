# ==========================================================================
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.general_prefix_of_class <- function() {
  prefix <- c(srn = "job_seurat5n",
    sr = "job_seurat",
    cc = "job_cellchat",
    lm = "job_limma",
    des = "job_deseq2",
    ven = "job_venn",
    vn = "job_vina",
    mn = "job_monocle2",
    rt = "job_reactome",
    mb = "job_mebocost",
    sc = "job_scenic",
    hd = "job_hdwgcna",
    ssr = "job_scissor"
  )
  if (any(duplicated(names(prefix)))) {
    stop('any(duplicated(names(prefix))), in `.general_prefix_of_class`')
  }
  prefix
}

