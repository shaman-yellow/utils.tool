# ==========================================================================
# workflow of risc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_risc <- setClass("job_risc", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://github.com/bioinfoDZ/RISC"),
    cite = "[@RobustIntegratLiuY2021]",
    method = "R package `RISC` used for scRNA-seq data integration"
    ))

setGeneric("asjob_risc", 
  function(x, ...) standardGeneric("asjob_risc"))

setMethod("asjob_risc", signature = c(x = "list"),
  function(x, filter, group.by){
    if (is.null(names(x))) {
      stop("is.null(names(x))")
    }
    lapply(x,
      function(x) {
        if (!is(x, "job_seurat")) {
          stop("is(x, 'job_seurat') == F")
        }
      })
    if (!is.null(filter)) {
      n <- 0L
      x <- lapply(x,
        function(x) {
          n <<- n + 1L
          x <- getsub(x, cells = grp(x@object@meta.data[[ group.by ]], filter))
          message("Object ", n, " after filtered: ", paste0(ids(x, group.by), collapse = ", "))
          return(x)
        })
    }
    message(crayon::red("Make sure all Gene names in the same ID levels. eg, all as hgnc_symbol."))
    x <- e(lapply(x,
        function(x) {
          use <- if (is(x, "job_seuratSp")) "Spatial" else "RNA"
          counts <- object(x)@assays[[ use ]]@counts
          rownames(counts) <- gs(rownames(counts), "\\.[0-9]*", "")
          counts <- counts[!duplicated(rownames(counts)), ]
          genes <- data.frame(index = 1:nrow(counts))
          rownames(genes) <- rownames(counts)
          cells <- object(x)@meta.data
          if (nrow(cells) != ncol(counts)) {
            counts <- counts[, colnames(counts) %in% rownames(cells)]
          }
          RISC::readsc(counts, cells, genes, is.filter = F)
        }))
    x <- .job_risc(object = x)
    x$group.by <- group.by
    return(x)
  })

setMethod("step0", signature = c(x = "job_risc"),
  function(x){
    step_message("Prepare your data with function `asjob_risc`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_risc"),
  function(x, workers = 5, minPC = 10, maxPC = 20){
    step_message("Preprocess before integration.
      red{{Note: in default, the data Quality was considered been done
      in procedure of Seurat reading, so was not conducted herein again.}}"
    )
    x$workers <- workers
    if (is.null(x$isNormed)) {
      object(x) <- e(lapply(object(x),
          function(object) {
            RISC::scFilter(object)
          }))
      object(x) <- e(lapply(object(x),
          function(object) {
            RISC::scNormalize(object, ncore = x$workers)
          }))
      object(x) <- e(lapply(object(x),
          function(object) {
            RISC::scDisperse(object)
          }))
      x$isNormed <- T
    }
    genes.ins <- ins(lst = lapply(object(x),
        function(object) object@rowdata$Symbol
        ))
    x$genes.ins <- genes.ins
    p.reference <- try(e(RISC::InPlot(object(x), var.gene = genes.ins,
          Std.cut = 0.99, ncore = x$workers, minPC = minPC, nPC = maxPC)), silent = T)
    if (inherits(p.reference, "try-error")) {
      warning("RISC::InPlot stopped.")
    } else {
      p.reference <- wrap(p.reference, 8, 9)
    }
    p.reference <- .set_lab(p.reference, sig(x), "select reference dataset for integration")
    x@plots[[ 1 ]] <- namel(p.reference)
    return(x)
  })

setMethod("step2", signature = c(x = "job_risc"),
  function(x, ref, pattern, colors = pgc(pattern, ids(x)), seed = 123)
  {
    step_message("Integration and transformed as Seurat object.")
    if (is.null(x$isIntegrate)) {
      if (!is.character(ref)) {
        stop("is.character(ref)")
      }
      ref.which <- which(names(object(x)) == ref)
      if (length(ref.which) > 1) {
        stop("length(ref.which) > 1")
      }
      object(x)[c(1, ref.which)] <- object(x)[c(ref.which, 1)]
      set.seed(seed)
      object(x) <- e(RISC::scMultiIntegrate(object(x),
          ncore = x$workers, add.Id = names(object(x)),
          var.gene = x$genes.ins))
      x$isIntegrate <- T
    }
    if (is.null(x$isUmap)) {
      object(x) <- e(RISC::scUMAP(object(x), npc = 20, use = "PLS"))
      x$isUmap <- T
    }
    p.umap <- vis(x, group.by = x$group.by, colors = colors)
    x$palette <- colors
    x@plots[[ 2 ]] <- namel(p.umap)
    return(x)
  })

setMethod("step3", signature = c(x = "job_risc"),
  function(x, group.by = x@params$group.by){
    step_message("Cell clustering and Find markers genes.")
    object(x) <- e(RISC::scCluster(object(x), slot = "cell.pls", neighbor = 3, npc = 20))
    clusters <- object(x)@coldata[[ group.by ]]
    if (!is.factor(clusters)) {
      clusters <- as.factor(clusters)
    }
    object(x)@coldata$Cluster <- clusters
    all_markers <- dplyr::as_tibble(e(RISC::AllMarker(object(x), ncore = x$workers)))
    x@tables[[ 3 ]] <- namel(all_markers)
    ## for other package
    ## object(x)@assay$logcount
    return(x)
  })

setMethod("step4", signature = c(x = "job_risc"),
  function(x, contrasts, group.by = x$group.by, p.adjust = .01, log2fc = 1){
    step_message("Test for DEGs.")
    if (is.data.frame(contrasts)) {
      contrasts <- apply(contrasts, 1, c, simplify = F)
    }
    res <- e(lapply(contrasts,
        function(con) {
          fun <- function(pt) {
            rownames(object(x)@coldata)[grpl(ids(x, group.by, F), pt)]
          }
          cell.ctrl <- fun(con[2])
          cell.sam <- fun(con[1])
          RISC::scDEG(object(x), cell.ctrl = cell.ctrl,
            cell.sam = cell.sam, ncore = x$workers, Padj = p.adjust, log2FC = log2fc)
        }))
    res <- lapply(res, dplyr::as_tibble)
    names(res) <- vapply(contrasts, function(x) paste0(x[1], "_vs_", x[2]), character(1))
    res <- dplyr::as_tibble(data.table::rbindlist(res, idcol = T))
    res <- dplyr::rename(res, contrast = .id)
    # res <- dplyr::filter(res, p_val_adj < .05)
    x@tables[[ 4 ]] <- list(contrasts = res)
    return(x)
  })

setMethod("focus", signature = c(x = "job_risc"),
  function(x, features, cluster = 3, group.by = x$group.by, colors = x$palette){
    step_message("Heatmap for DEGs.")
    ann_col <- if (is.null(colors)) NULL else list(Group = colors)
    RISC::Heat(object(x), genes = features, colFactor = group.by,
      ann_col = ann_col, gene.cluster = cluster)
  })

setMethod("vis", signature = c(x = "job_risc"),
  function(x, group.by = x$group.by, pt.size = .7, colors = color_set()) {
    e(RISC::DimPlot(object(x), colFactor = group.by,
        size = pt.size, Colors = colors))
  })

setMethod("ids", signature = c(x = "job_risc"),
  function(x, id = x@params$group.by, unique = T){
    ids <- object(x)@coldata[[ id ]]
    if (unique)
      ids <- unique(ids)
    ids
  })

pgc <- pattern_gradientColor <- function(pattern, names,
  colors = ggsci::pal_rickandmorty()(10))
{
  if (length(pattern) > length(colors)) {
    message("`colors` provided not enough.")
    colors <- color_set()
  }
  names <- as.character(unique(names))
  colors <- colors[ 1:length(pattern) ]
  n <- 0L
  palette <- lapply(pattern,
    function(pt) {
      n <<- n + 1L
      x <- sort(names[ grpl(names, pt) ])
      colors <- colorRampPalette(c("white", colors[n]))(length(x) + 2)[-(1:2)]
      names(colors) <- x
      return(colors)
    })
  unlist(palette)
}

setMethod("asjob_seurat", signature = c(x = "job_risc"),
  function(x, name = "integrated"){
    group.by <- x$group.by
    palette <- x$palette
    data <- do.call(cbind, object(x)@assay$logcount)
    data <- e(SeuratObject::CreateAssayObject(data = data))
    meta.data <- object(x)@coldata
    object <- e(SeuratObject::CreateSeuratObject(data, name, meta.data = meta.data))
    em.pca <- object(x)@DimReduction$cell.pls
    em.umap <- object(x)@DimReduction$cell.umap
    fun <- function(x, prefix) {
      colnames(x) <- gs(colnames(x), "^[a-zA-Z]*", prefix)
      return(x)
    }
    em.pca <- fun(em.pca, "PC_")
    em.umap <- fun(em.umap, "UMAP_")
    reduc.pca <- e(SeuratObject::CreateDimReducObject(em.pca, assay = name, key = "PC"))
    reduc.umap <- e(SeuratObject::CreateDimReducObject(em.umap, assay = name, key = "UMAP"))
    object@reductions$pca <- reduc.pca
    object@reductions$umap <- reduc.umap
    object@active.assay <- name
    x <- .job_seurat(object = object)
    x$group.by <- group.by
    x$palette <- palette
    x@step <- 3L
    sig(x) <- name
    return(x)
  })

setMethod("asjob_monocle", signature = c(x = "job_risc"),
  function(x){
    x <- asjob_seurat(x)
    x <- asjob_monocle(x)
    return(x)
  })

setMethod("mutate", signature = c(x = "job_risc"),
  function(x, ...){
    object(x)@coldata <- dplyr::mutate(object(x)@coldata, ...)
    return(x)
  })

setMethod("asjob_gsea", signature = c(x = "job_risc"),
  function(x, contrast.pattern = NULL, marker.list = x@tables$step4$contrasts)
  {
    if (!is.null(contrast.pattern)) {
      topTable <- dplyr::filter(marker.list, grpl(contrast, !!contrast.pattern))
    } else {
      topTable <- marker.list
    }
    job_gsea(dplyr::relocate(topTable, hgnc_symbol = Symbol, logFC = log2FC))
  })

# setMethod("asjob_monocle", signature = c(x = "job_risc"),
  # function(x, group.by = x@params$group.by, dims = 50){
  #   if (x@step < 3L) {
  #     stop("x@step < 3L")
  #   }
  #   group.by <- group.by
  #   palette <- x$palette
  #   counts <- do.call(cbind, object(x)@assay$logcount)
  #   gene_metadata <- object(x)@rowdata
  #   gene_metadata$gene_short_name <- rownames(gene_metadata)
  #   metadata <- object(x)@coldata
  #   object <- e(monocle3::new_cell_data_set(counts,
  #       cell_metadata = metadata, gene_metadata = gene_metadata))
  #   object <- e(monocle3::preprocess_cds(object, norm_method = "none", num_dim = dims))
  #   object <- e(monocle3::reduce_dimension(object))
  #   x <- .job_monocle(object = object)
  #   x$group.by <- group.by
  #   x$palette <- palette
  #   return(x)
  # })
