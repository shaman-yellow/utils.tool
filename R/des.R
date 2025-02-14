# ==========================================================================
# description for conception (colnames)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.get_des <- function(ref) {
  db <- .get_db_des()
  lst <- db[ names(db) %in% ref ]
  vapply(names(lst), function(x) lst[[ x ]], character(1))
}

.get_db_des <- function(fresh = FALSE) {
  if (is.null(maybe <- getOption("des")) | fresh) {
    des <- .parse_md()
    options(des = des)
    return(des)
  }
  maybe
}

.parse_md <- function(dir = .prefix(), pattern = "^help.*md$|description.*md$",
  excludes = c("id", ".id"))
{
  files <- list.files(dir, pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    return(list())
  }
  lst <- lapply(files,
    function(file) {
      lines <- readLines(file)
      if (grpl(name <- get_realname(file), "description")) {
        lines <- lines[ grpl(lines, "^[0-9\\-].*:") ]
        lines <- gs(lines, "^-\\s*|`", "")
      } else if (grpl(name, "help")) {
        lines <- sep_list(lines, "^$")
        lines <- unlist(lapply(lines, function(x) paste0(x, collapse = " ")))
        lines <- gs(lines, "^\\s*|\\s*$", "")
        lines <- gs(lines, "\\s{2,}", " ")
      }
      lines <- strsplit(lines, ":")
      lines
    })
  lst <- unlist(lst, recursive = FALSE)
  lst <- lapply(lst, function(x) if (length(x) > 1) x else NULL)
  lst <- lst_clear0(lst)
  names(lst) <- NULL
  lst <- nl(vapply(lst, function(x) x[1], character(1)),
    vapply(lst, function(x) x[2], character(1)))
  lst <- lst[ !duplicated(names(lst)) ]
  lst <- lst[ !names(lst) %in% excludes ]
  lst <- lapply(lst, function(x) gs(x, "([{}])", "\\\\\\1"))
  lst
}

.tag_anno <- function() {
  c(`pharm` = "网络药理学 (__dep__)",
    `herb:HERB` = "挖掘 HERB 数据库(__dep__)",
    `herb:BATMAN` = "挖掘 BATMAN (网络药理学) 数据库",
    `dis:genecards` = "获取疾病或条件相关的基因集:Genecards",
    `ppi:stringdb` = "构建 PPI 网络",
    `enrich:clusterProfiler` = "KEGG、GO富集分析",
    `docking:vina` = "全自动批量分子对接 (从蛋白质、分子名字到对接结果)",
    `scrna:anno` = "单细胞数据分析注释细胞群",
    `scrna:pseudo:monocle` = "单细胞拟时分析",
    `crd:crdscore` = "是否昼夜节律紊乱预测:rna (__dep__)",
    `scrna:cellchat` = "CellChat 单细胞细胞通讯",
    `db:TCGA` = "TCGA 数据挖掘和常规分析 (下载和整理)",
    `surv` = "生存分析",
    `wgcna` = "WGCNA 共表达分析",
    `filter:lasso` = "构建 LASSO-COX 预后模型",
    `rna:raw` = "从头开始的 RNA-seq 数据分析:从 fastq 开始质控、注释等",
    `rna:diff` = "limma, edgeR 差异分析",
    `gene:anno` = "biomaRt 基因注释、信息获取",
    `rbp` = "RBP 预测:RNA 结合蛋白与 RNA 的结合",
    `db:herb` = "HERB (网络药理学) 数据库",
    `16s:qiime2+mp` = "常规 16s RNA 肠道菌数据分析:alpha多样性+beta多样性+biomarker筛选等",
    `link:micro+meta+0` = "肠道菌数据和代谢物数据联合分析:通过公共数据挖掘建立联系 (__dep__)",
    `raw:16s` = "从公共数据库挖掘可用的 16s RNA 肠道菌群数据 (__dep__)",
    `raw:meta` = "从公共数据库挖掘可用的代谢物数据 (来自 GNPS 的原始数据 (较少)，或文献的补充材料) (__dep__)",
    `dock:vina` = "全自动批量分子对接",
    `dock:pep` = "分子对接肽与蛋白",
    `enrich:fella` = "代谢物富集分析",
    `wes:raw` = "从头开始的 WES 数据分析:从 fastq 开始质控、注释等",
    `wes:filter` = "WES 变异筛选",
    `annova` = "变异注释 (__dep__)",
    `raw:geo` = "从 GEO 挖掘可用数据:RNA-seq",
    `raw:tcga` = "TCGA RNA-seq 数据挖掘和常规分析",
    `enrich:pathview` = "Pathview 可视化富集通路",
    `rpackage` = "制作特定功能的 R 包",
    `raw:pdb` = "蛋白质 pdb 数据 (PDB，AlphaFold 互为补充) 获取用于可视化",
    `raw:protein` = "蛋白质信息获取",
    `scrna:cancer` = "单细胞数据鉴定癌细胞 (copyKAT)",
    `scrna:pseudo` = "单细胞拟时分析",
    `scrna:integrate:risc` = "不同来源的 scRNA-seq 数据整合:RISC整合消除批次效应 (__dep__)",
    `enrich` = "富集分析 (__dep__)",
    `raw:rna` = "从头开始的 RNA-seq 数据分析:从原始数据 fastq 开始质控、注释等",
    `diff` = "差异分析 (__dep__)",
    `link:meta+micro+0` = "肠道菌数据和代谢物数据联合分析:通过公共数据挖掘建立联系 (文献的补充材料) (__dep__)",
    `meta:fella` = "代谢物富集分析",
    `meta` = "代谢物数据分析",
    `link:protein+meta` = "蛋白质或基因和代谢物数据联合分析:通过公共数据挖掘建立联系 (文献的补充材料) (__dep__)",
    `get:fe` = "铁死亡相关基因获取 (数据库)",
    `enrich:gsea` = "GSEA 富集分析",
    `scrna:strna` = "空间转录组常规分析 (大致同单细胞数据) ",
    `scrna:anno:scsa` = "SCSA、ChatGPT-4 细胞注释",
    `scrna:monocle` = "Monocle3 拟时分析",
    `scrna:sub:pseudo` = "以拟时分析区分细胞亚群",
    `ppi` = "STRINGdb 构建 PPI 网络",
    `scrna:integrate` = "不同来源的 scRNA-seq 数据整合 (__dep__)",
    `rnaedit:edqtl` = "RNA 编辑位点信息挖掘:GTEx数据库edQTL挖掘 (__dep__)",
    `lite:bib` = "bibliometrix 常规文献计量分析和可视化",
    `link:micro+meta+1` = "肠道菌群和代谢物常规关联分析 (如 Pearson) :需来自于同一批样本 (__dep__)",
    `tcga:snp` = "TCGA 变异数据 (Variant) 挖掘和常规分析",
    `scrna:anno:sub` = "单细胞数据细胞注释到具体细胞亚群:如 Treg 细胞",
    `link:rna+rna` = "结合不同 (如疾病) 条件下 RNA-seq DEGs 筛选基因",
    `bindrna:circ+mi` = "circRNA 结合 miRNA 位点预测",
    `bindrna:mi+m` = "miRNA 结合 mRNA 位点预测",
    `db:tcmsp` = "TCMSP (网络药理学) 数据库数据挖掘",
    `crawl` = "编写爬虫工具批量获取网页数据",
    `target:swiss` = "SwissTargetPrediction 预测化合物靶点 (限制访问)",
    `link:meta+micro` = "肠道菌数据和代谢物数据联合分析:通过公共数据挖掘建立联系 (文献的补充材料)",
    `link:genes+meta` = "蛋白质或基因和代谢物数据联合分析:通过公共数据挖掘建立联系 (文献的补充材料)",
    `link:genes+meta+micro` = "蛋白质或基因、肠道菌群、代谢物数据联合分析:通过公共数据 (文献的补充材料)",
    `link:genes(rna)+eqtl` = "基因 (RNA) 和 eQTL 关联性挖掘:通过公共数据 (文献的补充材料)",
    `link:eqtl+snp` = "eQTL 和 SNP 关联性挖掘:通过公共数据 (__dep__)",
    `link:snp+meta` = "SNP 和代谢物之间的关联性挖掘:通过公共数据 (文献的补充材料) (__dep__)",
    `dock:protein+protein` = "ClusPro 蛋白质与蛋白质之间的对接模拟",
    `db:plant` = "从 Plantaedb 数据库获取植物成分信息 (作为中药数据库缺失记录的补充) ",
    `fe` = "铁死亡相关基因挖掘",
    `hob` = "20% HOB 口服生物利用度预测",
    `m6a:site` = "m6A 编辑位点获取",
    `cpd:info` = "化合物信息获取:同义名、smiles、inchikey等",
    `target` = "数据库获取化合物靶点",
    `target:sp` = "预测化合物靶点",
    `dis:pharm` = "疾病相关基因集:PharmGKB 数据库挖掘 (__dep__)",
    `dis:dis` = "疾病相关基因集:DisGeNet 数据库挖掘 (__dep__)",
    `tf` = "调控该基因的相关转录因子 (TF) 数据获取",
    `cpd:class` = "化合物系统分类注释 (ClassyFire) ",
    `link:meta+genes+0` = "蛋白质或基因和代谢物数据联合分析:通过公共数据挖掘建立联系 (__dep__)",
    `lcms` = "从头开始的 LC-MS/MS 数据处理:数据格式转换、峰检测、化合物鉴定等",
    `lcms:sirius` = "LC-MS/MS 数据以尖端预测工具鉴定更多化合物 (费时)",
    `lcms:deep` = "LC-MS/MS 数据深入分析:差异代谢物筛选、分子网络可视化、富集分析等",
    `ocr` = "批量 OCR 识别图片中文字并整理结果",
    `fastq` = "fastq 文件预处理、对齐参考基因",
    `metagenome` = "宏基因组数据分析",
    `fea` = "特征筛选 (__dep__)",
    `gwas` = "GWAS 数据获取和分析",
    `mix` = "多数据库获取疾病靶点 (__dep__)",
    `mass` = "质谱数据分析:原始数据预处理、特征峰解卷积等",
    `unitmp` = "跨膜蛋白数据获取",
    `prr` = "药物临床化疗反应预测、耐药性预测",
    `musite` = "蛋白转录后修饰位点预测 (PTMs)",
    `epi` = "表观遗传调控因子数据获取 (数据库)",
    `ctd` = "化合物和疾病关联数据获取 (CTD 数据库)",
    `dl` = "药物 Drug-likeness 预测 (__dep__)",
    `scrna:flux` = "单细胞数据预测代谢通量",
    `methyl:array` = "甲基化芯片数据分析",
    `methyl:dmr` = "差异甲基化分析、可视化",
    `cpg` = "CpG islands 数据获取",
    `ucsc` = "UCSC table 数据获取、基因位点注释",
    `cardinal` = "空间代谢组数据分析:Cardinal",
    `maf` = "变异注释分析和可视化:maftools",
    `bt2` = "基因组比对:Bowtie 2",
    `mfuzz` = "对基因聚类分析 (时序分析) mfuzz",
    `pubmed` = "PubMed 文献搜索",
    `gds` = "GEO 数据库搜索、调研，获取可用数据",
    `pathview` = "KEGG 通路可视化和表达量映射",
    `xena` = "UCSCXenaTools 癌症相关数据获取 (相较于 TCGA 补充了 GTEx 和 Target 数据库，适合 Cancer VS Normal)",
    `fusion` = "TWAS全转录组关联研究",
    `venn` = "取交集",
    `vep` = "变异注释 (SNP 获取 rsID 等注释)",
    `seurat5n` = "单细胞数据分析",
    `lnctard` = "LncRNA调控靶点"
  )
}
