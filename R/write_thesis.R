# ==========================================================================
# tools for fastly write formatting thesis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inclu.fig <- function(image, land = FALSE, saveDir = "thesis_fig", dpi = 300,
  scale = if (land) {
    list(width = 10, height = 6.2)
  } else {
    list(width = 6.2, height = 10)
  }, name = NULL)
{
  if (!file.exists(saveDir))
    dir.create(saveDir)
  upper <- dirname(image)
  file <- paste0(
    digest::digest(image, "crc32c"), "_", basename(image)
  )
  ## backup for figure
  if (grepl("\\.pdf$", file)) {
    savename <- file.path(saveDir, sub("\\.pdf$", ".png", file))
    if (!file.exists(savename)) {
      pdf_convert(image, filenames = savename, dpi = dpi, pages = 1)
      need_trim <- TRUE
    } else need_trim <- FALSE
  } else {
    savename <- file.path(saveDir, file)
    if (!file.exists(savename)) {
      file.copy(image, savename)
      need_trim <- TRUE
    } else need_trim <- FALSE
  }
  ## trim the border
  if (!need_trim) {
    gc()
    pic_trim(savename)
  }
  ## record image info
  image_info <- getOption("Inclu_image_info", list())
  if (is.null(image_info[[ savename ]])) {
    img <- magick::image_read(savename)
    info <- magick::image_info(img)
    magick::image_destroy(img)
    image_info[[ savename ]] <- list(width = info$width, height = info$height)
    options(Inclu_image_info = image_info)
  } else {
    info <- image_info[[ savename ]]
  }
  ## set figure display
  ratio <- info$width / info$height
  if ((scale$width / scale$height) >= ratio) {
    ## height as reference
    height <- scale$height
    width <- height * ratio
  } else {
    ## width as reference
    width <- scale$width
    height <- width / ratio
  }
  knitr::opts_current$set(fig.height = height, fig.width = width)
  if (knitr::pandoc_to("docx")) {
    content <- officer::external_img(savename, width, height)
    writeLines(c("", assis_docx_img(content), ""))
    if (missing(name)) {
      label <- knitr::opts_current$get("autor_label")
    } else {
      label <- name
    }
    if (is.null(label)) {
      stop(glue::glue('is.null(autor_label), Chunk: {knitr::opts_current$get("label")}'))
    }
    if (!grpl(label, "^unnamed-chunk")) {
      run_num <- officer::run_autonum(
        seq_id = knitr::opts_chunk$get("fig.lp"),
        pre_label = knitr::opts_chunk$get("fig.cap.pre"),
        post_label = knitr::opts_chunk$get("fig.cap.sep"),
        bkm = label,
        prop = officer::fp_text_lite(bold = TRUE)
      )
      caption <- officer::block_caption(as_caption(label), "Image Caption", autonum = run_num)
      content <- officer:::to_wml_block_caption_pandoc(caption)
      writeLines(c("", content, ""))
    }
  } else {
    knitr::include_graphics(savename)
  }
}

prepare.fig <- function(image, saveDir = "thesis_fig", dpi = 300)
{
  if (!file.exists(saveDir))
    dir.create(saveDir)
  upper <- dirname(image)
  if (is.na(upper))
    upper <- "."
  file <- basename(image)
  ## backup for figure
  if (grepl("\\.pdf$", file)) {
    savename <- paste0(saveDir, "/", sub("\\.pdf$", ".png", file))
    if (!file.exists(savename)) {
      pdf_convert(image, filenames = savename, dpi = dpi, pages = 1)
      need_trim <- TRUE
    } else need_trim <- FALSE
  } else {
    savename <- paste0(saveDir, "/", file)
    if (!file.exists(savename)) {
      file.copy(image, savename)
      need_trim <- TRUE
    } else need_trim <- FALSE
  }
  ## trim the border
  if (need_trim) {
    gc()
    pic_trim(savename)
  }
}

inclu.capt <- function(img, saveDir = "thesis_fig") {
  if (!file.exists(saveDir)) {
    dir.create(saveDir)
  }
  filename <- basename(img)
  savename <- paste0(saveDir, "/", filename)
  if (!file.exists(savename)) {
    png(savename, 1000, 1000, res = 150)
    grid::grid.raster(png::readPNG(img))
    dev.off()
  }
  inclu.fig(savename)
}

.env <- topenv(environment())

pic_trim <- function(file){
  img <- magick::image_read(file)
  img <- magick::image_trim(img)
  magick::image_write(img, file)
  magick::image_destroy(img)
}

pretty_flex2 <- function(data,
  caption = "This is table", footer = NULL,  bold_header = FALSE,
  weight = sc(colnames(data), rep(1, length(data))), width = 6,
  family = "Times New Roman", e.family = "SimSun", font.size = 10.5,
  caption_style = "Table Caption", form_body = FALSE, form_header = TRUE)
{
  args <- as.list(environment())
  do.call(pretty_flex, args)
}

pretty_flex <- function(data,
  caption = "This is table", footer = "This is footer",  bold_header = FALSE,
  weight = sc(colnames(data), rep(1, length(data))), width = 6,
  family = "Times New Roman", e.family = "SimSun", font.size = 10.5,
  caption_style = "Table Caption", form_body = TRUE, form_header = TRUE)
{
  if (form_body) {
    data <- dplyr::mutate_if(
      data, function(x) ifelse(is.character(x) | is.factor(x), TRUE, FALSE), form
    )
  }
  if (form_header) {
    colnames(data) <- form(colnames(data))
    if (!is.null(names(weight)))
      names(weight) <- form(names(weight))
  }
  weight <- as.list(weight)
  weight <- vapply(colnames(data),
    function(name) {
      if (is.null(weight[[ name ]])) 1 else weight[[ name ]]
    }, double(1), USE.NAMES = FALSE)
  w.width <- weight * (width / sum(weight))
  data <- flextable::flextable(data)
  data <- flextable::valign(data, valign = "top", part = "body")
  if (bold_header) {
    data <- flextable::bold(data, part = "header", bold = TRUE)
  }
  if (!is.null(footer)) {
    if (is.character(footer)) {
      footer <- paste0(footer, collapse = "\n")
    }
    data <- flextable::add_footer_lines(data, footer)
  }
  data <- flextable::width(data, width = w.width)
  data <- flextable::font(
    data, fontname = family, eastasia.family = e.family, part = 'all'
  )
  data <- flextable::fontsize(data, size = font.size, part = 'all')
  data <- flextable::set_caption(data, caption, word_stylename = caption_style)
  data
}

flex_footer <- function(x, props = fp_text(font.size = 10.5, font.family = "SimSun")) {
  flextable::as_chunk(x, props = props)
}

match_fun2 <- function(txt, pattern = "^[0-9_.A-Za-z]*(?= <-)"){
  res <- stringr::str_extract(txt, pattern)
  res <- res[ !vapply(res, is.na, logical(1)) ]
}

match_method2 <- function(txt){
  match_fun2(txt, "(?<=setMethod\\(\")[^\"]*(?=\")")
}

match_class <- function(txt){
  class <- match_fun2(txt, "(?<=setClass\\(\").*(?=\")")
  sapply(class,
    function(x) {
      if (!isVirtualClass(x))
        slotNames(new(x))
    }, simplify = FALSE)
}

as_code_list <- function(names, value = rep("", length(names)), prefix = "c"){
  writeLines(paste0(prefix, "("))
  lapply(seq_along(names),
    function(n) {
      end <- if (n == length(names)) "\"" else "\","
      writeLines(paste0("  `", names[[ n ]], "` = \"", value[[ n ]], end))
    })
  writeLines(")")
}

as_code_list2 <- function(lst){
  lapply(seq_along(lst),
    function(n) {
      writeLines(paste0("## ", names(lst[n])))
      as_code_list(rep(names(lst[n]), length(lst[[ n ]])),
        lst[[ n ]])
      writeLines("")
    })
}

as_pander.flex <- function(obj, page.width = 12) {
  require(pander)
  data <- obj$body$dataset
  caption <- obj$caption$value
  pandoc.table(data, caption)
}

as_kable.flex <- function(obj, page.width = 14, landscape.width = 19) {
  require(kableExtra)
  data <- obj$body$dataset
  if (length(n <- grep("InChIKey", colnames(data))) > 0) {
    data[[ n ]] <- gsub("([A-Z]{7})([A-Z]{1,})", "\\1- \\2", data[[ n ]])
  }
  caption <- obj$caption$value
  kable <- kable(
    data, "latex", booktabs = TRUE,
    longtable = TRUE, caption = caption
  )
  if (ncol(data) > 8) {
    page.width <- landscape.width
  }
  widths <- obj$body$colwidths
  widths <- as.list(
    paste0(widths / sum(widths) * page.width, "cm")
  )
  for (i in seq_along(widths)) {
    kable <- column_spec(kable, i, width = widths[[ i ]])
  }
  kable <- kable_styling(
    kable, latex_options = c("hold_position", "repeat_header"),
    font_size = 10.5
  )
  if (length(footnotes <- obj$footer$dataset[[ 1 ]]) > 0 ) {
    footnotes <- strsplit(footnotes, split = "\n")[[1]]
    kable <- kableExtra::footnote(
      kable, number = footnotes, fixed_small_size = TRUE
    )
  }
  kable
}

as_gt.flex <- function(obj, page.width = 650) {
  require(gt)
  data <- obj$body$dataset
  colnames <- colnames(data)
  widths <- obj$body$colwidths
  widths <- as.list(widths / sum(widths) * page.width)
  widths <- lapply(seq_along(widths),
    function(n) {
      eval(parse(text = paste0("`", colnames[n], "`", " ~ ",
            "px(", widths[n], ")")))
    })
  caption <- obj$caption$value
  gt <- pretty_table(
    data, title = NULL, footnote = NULL,
    filename = NULL, caption = md(paste0("**", caption, "**")), widths = widths
  )
  if (length(footnotes <- obj$footer$dataset[[ 1 ]]) > 0 ) {
    if (grepl(":|：", footnotes)) {
      footnotes <- strsplit(footnotes, split = "\n")[[1]]
      locates <- stringr::str_extract(footnotes, "^.*?(?=:|：)")
      footnotes <- gsub("^.*?[:：]", "", footnotes)
      for (i in seq_along(footnotes)) {
        gt <- tab_footnote(gt, footnote = footnotes[i],
          locations = cells_column_labels(columns = dplyr::starts_with(locates[i])))
      }
    } else {
      gt <- tab_footnote(gt, footnote = footnotes)
    }
  }
  gt <- tab_options(gt, table.font.size = px(10.5), footnotes.font.size = px(8))
  if (needTex()) {
    chunk_label <- knitr::opts_current$get("tab.id")
    gt <- as_latex_with_caption(gt, chunk_label, caption)
  }
  gt
}

as_latex_with_caption <- function(gt, chunk_label, caption) {
  gt <- gt::as_latex(gt)
  caption <- paste0("\\caption{\\label{tab:", chunk_label, "}", caption, "}\\\\")
  latex <- strsplit(gt[1], split = "\n")[[1]]
  pos <- grep("\\\\begin\\{longtable\\}", latex)
  # pos.end <- grep("\\\\end\\{longtable\\}", latex)
  if (!length(pos) == 1)
    stop("length(pos) != 1")
  # if (!length(pos.end) == 1)
    # stop("length(pos.end) != 1")
  latex[pos] <- paste0(latex[pos], "\n", caption, "\n")
  # latex[pos] <- paste0("\\resizebox{\\linewidth}{!}{\n", latex[pos], "\n", caption, "\n")
  # latex[pos.end] <- paste0(latex[pos.end], "\n", "}")
  latex <- paste(latex, collapse = "\n")
  gt[1] <- latex
  return(gt)
}

thesis_as_latex.word <- function(md, flex_as_kable = "(^flex_.*)") {
  if (length(md) > 1) {
    message("Guess `md` has been read by readLines.")
  } else {
    md <- readLines(md)
    md <- gsub("<!---BLOCK_LANDSCAPE_START--->", "\\\\landstart", md)
    md <- gsub("<!---BLOCK_LANDSCAPE_STOP--->", "\\\\landend", md)
    md <- gsub("<!---BLOCK_TOC--->", "\\\\tableofcontents", md)
    md <- gsub("<!---BLOCK_TOC\\{seq_id: 'fig'\\}--->", "\\\\listoffigures", md)
    md <- gsub("<!---BLOCK_TOC\\{seq_id: 'tab'\\}--->", "\\\\listoftables", md)
    md <- gsub("(\\*\\*.*?[:：]\\*\\*)", "\\\\noindent\n\\1", md)
    if (!is.null(flex_as_kable)) {
      md <- gsub(flex_as_kable, "as_kable.flex(\\1)", md)
    }
    md <- gsub("^(```\\{r.*)tab.id", "\\1label", md)
    # md <- gsub("tab.id\\s*=\\s*\"(.*?)\"", "\\1", md)
    md <- gsub("&emsp;&emsp;\\s*", "", md)
  }
}

insert_pdf.tex <- function(lines, pos, pdf, set_pagenumber = TRUE) {
  if (set_pagenumber) {
    command <- "\\includepdfset{pagecommand={\\thispagestyle{fancy}}}"
    lines[pos] <- paste0(lines[pos], "\n", command)
  }
  pages <- paste0("1-", pdftools::pdf_info(pdf)$pages)
  command <- paste0("\\includepdfmerge{", pdf, ",", pages, "}\n")
  lines[pos] <- paste0(lines[pos], "\n", command)
  return(lines)
}

insert_tocCont.tex <- function(lines, pos, name, style = "section") {
  command <- paste0("\n\\newpage\n", "\\addcontentsline{toc}{", style, "}{", name, "}")
  lines[pos] <- paste0(lines[pos], "\n", command)
  return(lines)
}

# ==========================================================================
# insert text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (requireNamespace("officer", quietly = TRUE)) {
  ftext <- officer::ftext
  fp_text <- officer::fp_text
  fpar <- officer::fpar
  fp_par <- officer::fp_par
  st.index <- fp_text(font.size = 14, font.family = "SimHei")
}

# ==========================================================================
# custom render
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

custom_render <- function(file, theme = 'thesis', fix = fix_spell,
  fun_render = rmarkdown::render, heading_mark = "#") {
  ## initialize
  global <- environment()
  clear <- function(ns) {
    lapply(ns, function(n) assign(paste0("n", n), 0, envir = global))
  }
  clear(1:9)
  attr(n1, 'part') <- 0
  ## numbering heading
  md <- fix(readLines(file))
  h.pos <- grep(paste0("^", heading_mark), md)
  chunk.pos <- get_chunkPos(md)
  h.pos <- h.pos[ !h.pos %in% chunk.pos ]
  toc <- lapply(md[h.pos],
    function(ch) {
      level <- length(stringr::str_extract_all(ch, heading_mark)[[ 1 ]])
      fun <- match.fun(paste0("rheader", ifelse(level > 4, 4, level), ".", theme))
      nch <- fun(ch, level = level, global = global, clear = clear)
      if (grepl("@@part", ch)) level <- 0
      list(heading = nch, level = level)
    })
  md[h.pos] <- vapply(toc, function(x) x$heading, character(1))
  ## cross-reference
  ids <- stringr::str_extract(md, "(?<=\\{@@id:)[a-zA-Z._0-9]*(?=\\})")
  ids.pos <- (seq_along(md))[ !is.na(ids) ]
  if (length(ids.pos) > 0) {
    ids.relToc <- vapply(ids.pos, FUN.VALUE = character(1),
      function(pos) {
        as_ref.ch(search_toc(pos, toc = toc, toc.pos = h.pos))
      })
    ids.db <- as.list(sc(ids[ !is.na(ids) ], ids.relToc))
    md[ ids.pos ] <- gsub("\\{@@id:[a-zA-Z._0-9]*\\}", "", md[ ids.pos ])
    refs <- stringr::str_extract_all(md, "(?<=\\{@@ref:)[a-zA-Z._0-9]*(?=\\})")
    refs <- unique(unlist(refs))
    refs <- refs[ !is.na(refs) ]
    if (length(refs) > 0) {
      refs.content <- vapply(refs, FUN.VALUE = character(1),
        function(ref) {
          content <- ids.db[[ ref ]]
          if (is.null(content)) {
            stop("The reference id: ", ref, " not found.")
          }
          content
        })
      for (i in seq_along(refs)) {
        md <- gsub(paste0("{@@ref:", refs[i], "}"), refs.content[i], md, fixed = TRUE)
      }
    }
  }
  ## output
  filename <- basename(file)
  filepath <- dirname(file)
  writeLines(md, nfile <- paste0(filepath, "/", "_temp_", filename))
  if (!is.null(fun_render))
    fun_render(nfile)
}

get_chunkPos <- function(md) {
  chunk.pos <- grep("^```", md)
  unlist(lapply(seq(1, length(chunk.pos), by = 2),
      function(n) {
        chunk.pos[n]:chunk.pos[ n + 1 ]
      }))
}

## chinese \u4e00-\u9fa5
## double [^\x00-\xff] 

pch <- function() "[\u4e00-\u9fa5]"

fix_spell <- function(md) {
  chunk.pos <- get_chunkPos(md)
  yaml.pos <- grep("^---", md)
  yaml.pos <- 1:(yaml.pos[2])
  pos <- (seq_along(md))
  pos <- pos[ !pos %in% c(chunk.pos, yaml.pos) ]
  md[ pos ] <- gsub("\"",  "\'", md[ pos ])
  md[ pos ] <- fix_quote(md[ pos ])
  md
}

fix_quote <- function(lines, left = "‘", right = "’", extra = "：、，。+（）") {
  pattern <- paste0("[\u4e00-\u9fa5", extra, "]")
  lines <- gsub(paste0("(?<=", pattern, ")\\s*\'\\s*(?=[a-zA-Z])"), left, lines, perl = TRUE)
  lines <- gsub(paste0("(?<=[a-zA-Z0-9])\\s*\'\\s*(?=", pattern, ")"), right, lines, perl = TRUE)
  lines
}

fix_quote.latex <- function(lines) {
  lines <- gsub('‘', '`', lines)
  lines <- gsub('’', '\'', lines)
  lines <- gsub("\'([a-zA-Z_.\\(\\)0-9]{1,})\'", "`\\1\'", lines)
  lines <- gsub("`([a-zA-Z_.\\(\\)0-9]{1,})`", "`\\1\'", lines)
  lines
}

as_ref.ch <- function(rel.toc, sep = " &gt; ") {
  if (is.list(rel.toc)) {
    rel.toc <- vapply(rel.toc, function(x) x$heading, character(1))
  }
  heading <- gsub("^#*\\s*|\\s*\\{-\\}", "", rel.toc)
  paste0(heading, collapse = sep)
}

search_toc <- function(pos, toc, toc.pos)
{
  previous <- toc[toc.pos <= pos]
  pre.levels <- vapply(previous, function(x) x$level, double(1))
  top.level <- tail(pre.levels, n = 1) + 1
  is.rel <- rev(vapply(rev(pre.levels), FUN.VALUE = logical(1),
    function(level) {
      if (level < top.level) {
        top.level <<- level
        return(TRUE)
      } else FALSE
    }))
  previous[ is.rel ]
}

rheader1.thesis <- function(x, global, clear, ...,
  fun.part = function() {
    x <- sub("# @@part ", paste0("# ", "第", chn(attr(n, 'part')), "部分 "), x)
    paste0(x, " {-}")
  },
  fun.normal = function() {
    paste0(sub("# ", paste0("# ", chn(n), "、"), x), "{-}")
  })
{
  if (grepl('\\{-\\}', x)) {
    return(x)
  } else if (grepl('@@part ', x)) {
    n <- 0
    n1 <- get("n1", envir = global)
    attr(n, 'part') <- attr(n1, 'part') + 1
    clear(1:9)
    assign("n1", n, envir = global)
    fun.part()
  } else {
    n <- get("n1", envir = global) + 1
    clear(1:9)
    assign("n1", n, envir = global)
    fun.normal()
  }
}

rheader2.thesis <- function(x, global, clear, ...,
  fun = function() {
    sub("## ", paste0("## ", "（", chn(n), "）"), sub("$", " {-}", x))
  })
{
  n <- get("n2", envir = global) + 1
  clear(2:9)
  assign("n2", n, envir = global)
  fun()
}

rheader3.thesis <- function(x, global, clear, ...,
  fun = function() {
    sub("### ", paste0("### ", n, ". "), sub("$", " {-}", x))
  }) 
{
  n <- get("n3", envir = global) + 1
  clear(3:9)
  assign("n3", n, envir = global)
  fun()
}

rheader4.thesis <- function(x, global, clear, level = 4, ...) {
  n <- get(paste0("n", level), envir = global) + 1
  clear(level:9)
  assign(paste0("n", level), n, envir = global)
  num <- vapply(3:level, function(n) get(paste0("n", n), envir = global), double(1))
  num <- paste0(num, collapse = ".")
  sub("#* ",
    paste0(paste0(rep("#", level), collapse = ""), " ", num, " "),
    sub("$", " {-}", x))
}

chn <- function(n) {
  switch(n,
    "1" = "一", "2" = "二", "3" = "三", "4" = "四", "5" = "五",
    "6" = "六", "7" = "七", "8" = "八", "9" = "九", "10" = "十"
  )
}

# ==========================================================================
# custom docx tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

custom_docx_document <- function(...){
  args <- list(...)
  expr <- args$reference_docx
  if (grepl("^`r ", expr)) {
    path <- eval(parse(text = gsub("^`r |`$", "", expr)))
    args$reference_docx <- path
  }
  args <- list(
    reference_docx = args$reference_docx,
    keep_md = TRUE,
    df_print = "tibble",
    tables = list(
      style = "Table", width = 1, topcaption = TRUE, tab.lp = "tab:",
      caption = list(
        style = "Table Caption", pre = "表", sep = " ", tnd = 0, tns = "-",
        fp_text = structure(list(
          font.size = NA, bold = TRUE, italic = NA,
          underlined = NA, color = NA, font.family = NA, bold.cs = FALSE,
          font.size.cs = NA, vertical.align = 'baseline', shading.color = NA,
          hansi.family = NA, eastasia.family = NA, cs.family = NA,
          lang.val = NA, lang.eastasia = NA, lang.bidi = NA
        ), structure = "fp_text")
      )),
    plots = list(
      style = "Normal", align = "center", fig.lp = "fig:", topcaption = FALSE,
      caption = list(
        style = "Image Caption", pre = "图", sep = " ", tnd = 0, tns = "-",
        fp_text = structure(list(
          font.size = NA, bold = FALSE, italic = NA,
          underlined = NA, color = NA, font.family = NA, bold.cs = FALSE,
          font.size.cs = NA, vertical.align = 'baseline', shading.color = NA,
          hansi.family = NA, eastasia.family = NA, cs.family = NA,
          lang.val = NA, lang.eastasia = NA, lang.bidi = NA
        ), structure = "fp_text")
        )),
    page_margins = list(
      bottom = 1, top = 1, right = 1.25, left = 1.25, header = 0.5,
      footer = 0.5, gutter = 0.5)
  )
  do.call(officedown::rdocx_document, c(list('bookdown::word_document2'), args))
}

custom_docx_document2 <- function(...){
  args <- list(...)
  expr <- args$reference_docx
  if (grepl("^`r ", expr)) {
    path <- eval(parse(text = gsub("^`r |`$", "", expr)))
    args$reference_docx <- path
  }
  args <- list(
    reference_docx = args$reference_docx,
    keep_md = TRUE,
    df_print = "tibble",
    tables = list(
      style = "Table", width = 1, topcaption = TRUE, tab.lp = "tab:",
      caption = list(
        style = "Table Caption", pre = "Tab. ", sep = " ", tnd = 0, tns = "-",
        fp_text = structure(list(
          font.size = NA, bold = TRUE, italic = NA,
          underlined = NA, color = NA, font.family = NA, bold.cs = FALSE,
          font.size.cs = NA, vertical.align = 'baseline', shading.color = NA,
          hansi.family = NA, eastasia.family = NA, cs.family = NA,
          lang.val = NA, lang.eastasia = NA, lang.bidi = NA
        ), structure = "fp_text")
      )),
    plots = list(
      style = "Normal", align = "center", fig.lp = "fig:", topcaption = FALSE,
      caption = list(
        style = "Image Caption", pre = "Fig. ", sep = " ", tnd = 0, tns = "-",
        fp_text = structure(list(
          font.size = NA, bold = TRUE, italic = NA,
          underlined = NA, color = NA, font.family = NA, bold.cs = FALSE,
          font.size.cs = NA, vertical.align = 'baseline', shading.color = NA,
          hansi.family = NA, eastasia.family = NA, cs.family = NA,
          lang.val = NA, lang.eastasia = NA, lang.bidi = NA
        ), structure = "fp_text")
        )),
    page_margins = list(
      bottom = 1, top = 1, right = 1.25, left = 1.25, header = 0.5,
      footer = 0.5, gutter = 0.5)
  )
  output_formats <- do.call(officedown::rdocx_document, c(list('bookdown::word_document2'), args))
  original_knit <- output_formats$post_knit
  ## replace the post_knit function, add post-knit ...
  new_fun <- function(metadata, input_file, runtime, ...) {
    original_knit(metadata, input_file, runtime, ...)
    output_file <- officedown:::file_with_meta_ext(input_file, "knit", "md")
    if (!file.exists(output_file)) {
      stop('!file.exists(output_file), I can not found the post-knit md file.')
    }
    lines <- readLines(output_file)
    all_placeHolder <- stringr::str_extract_all(lines, "<!-- autor_legend:[a-zA-Z0-9-]*? -->")
    n <- 0L
    autor_legend_env <- getOption("autor_legend_env")
    lines <- vapply(lines,
      function(line) {
        n <<- n + 1L
        if (length(alls <- all_placeHolder[[ n ]])) {
          replace <- ""
          label <- ""
          for (i in alls) {
            label <- stringr::str_extract(i, "(?<=autor_legend:)[a-zA-Z0-9-]*")
            replace <- autor_legend_env[[ label ]]
            if (!is.null(replace) && length(replace)) {
              line <- sub(i, replace, line)
              # only for once
              autor_legend_env[[ label ]] <- NULL
            }
          }
        }
        return(line)
      }, character(1))
    writeLines(lines, output_file)
  }
  output_formats$post_knit <- new_fun
  output_formats
}



# ==========================================================================
# xml tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

touch_docSt <- function(file = "templ_zcmu_thesis.docx",
  path = ".", to = paste0(path, "/", ".templ"))
{
  utils::unzip(paste0(path, "/", file), exdir = to)
  style <- paste0(to, "/word/styles.xml")
  attr(style, "dir") <- to
  style
}

update_docSt <- function(st, xml){
  XML::saveXML(xml, st, prefix = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n')
  file <- readLines(st)
  content <- c(file[1], "")
  xml <- paste0(gsub("^\\s*|\\s*$", "", file[-1]), collapse = "")
  pat <- unique(unlist(stringr::str_extract_all(xml, "(?<=\\s)[a-zA-Z]*=")))
  for (i in pat) {
    xml <- gsub(paste0("\\s", i), paste0(" ", "w:", i), xml)
  }
  content[2] <- xml
  writeLines(content, st)
}

pack_docSt <- function(st, file = "temple.docx"){
  if (!file.exists(file)) {
    file.create(file)
  }
  to_file <- normalizePath(file)
  file.remove(file)
  wd <- getwd()
  setwd(attr(st, "dir"))
  utils::zip(to_file, "./")
  setwd(wd)
}

readxml <- function(file){
  xml <- XML::xmlRoot(XML::xmlTreeParse(file))
  writeLines(paste0("XML length:", length(xml)))
  xml
}

ft <- function(xml, target){
  names <- XML::xmlSApply(xml, xmlName)
  unname(which(names == target))
}

sc <- function(name, value) {
  names(value) <- name
  return(value)
}

allnames <- function(xml2) {
  vapply(seq_along(xml2), function(n) XML::xmlName(xml2[[ n ]]), character(1))
}

cp_xml.style <- function(xml, target = "Table", src = "Plain Table 1"){
  src.n <- st.search(xml, src)
  target.n <- st.search(xml, target)
  xml[[ target.n ]] <- XML::removeChildren(
    xml[[ target.n ]], kids = allnames(xml[[ target.n ]])[ -1 ]
  )
  xml[[ target.n ]] <- do.call(
    XML::addChildren,
    c(list(xml[[ target.n ]]), XML::xmlChildren(xml[[ src.n ]])[ -1 ])
  )
  return(xml)
}

st.search <- function(xml, name) {
  unlist(lapply(seq_along(xml),
      function(n) {
        if (st.name(xml[[ n ]]) == name) n
      }))
}

st.name<- function(xml2){
  nametag <- xml2[[ "name" ]]
  if (!is.null(nametag)) {
    return(XML::xmlAttrs(nametag)[[ "val" ]])
  }
  return("_NULL_")
}

repl_xml.fontShd <- function(xml,
  param = "fill", as = "ededed", name = show_allTok(xml),
  target = "rPr", ttag = "shd")
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

insert_xml.tabs <- 
  function(xml,
    param = c("val", "pos"), as = c("left", 240), name = paste0("heading ", 1:9),
    target = "pPr", ttag = "tabs", namespace = "w",
    insert_ttag = XML::xmlNode("tabs",
      XML::xmlNode("tab", attrs = sc(param, as), namespace = namespace),
      namespace = namespace))
  {
    args <- as.list(environment())
    args$match.arg <- FALSE
    do.call(repl_xml.font, args)
  }

insert_xml.praNum <- 
  function(xml,
    param = "val", as = 20, name = paste0("heading ", 1:9),
    target = "pPr", ttag = "numPr", namespace = "w",
    insert_ttag = XML::xmlNode("numPr",
      XML::xmlNode("numId", attrs = sc(param, as), namespace = namespace),
      namespace = namespace))
  {
    args <- as.list(environment())
    args$match.arg <- FALSE
    do.call(repl_xml.font, args)
  }

repl_xml.fontSpace <- function(xml,
  param = "after", as = "0", name = "caption",
  main = "style", target = "pPr", ttag = "spacing", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = FALSE, insert_param = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

repl_xml.fontOutl <- function(xml,
  param = "val", as = 0, name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  main = "style", target = "pPr", ttag = "outlineLvl", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = FALSE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

repl_xml.fontAlign <- function(xml,
  param = "val", as = "center", name = "caption",
  main = "style", target = "pPr", ttag = "jc", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = FALSE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

repl_xml.fontSz <- function(xml,
  param = "val", as = "18", name = "heading 1",
  main = "style", target = "rPr", ttag = "sz", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = FALSE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

repl_xml.fontSzCs <- function(xml,
  param = "val", as = "24", name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  main = "style", target = "rPr", ttag = "szCs", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = FALSE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

dele_xml.lang <- function(xml,
  name = "Normal", target = "rPr", ttag = "lang", delete_ttag = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

dele_xml.ital <- function(xml,
  name = "caption", target = "rPr", ttag = "i", delete_ttag = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

dele_xml.szCs <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "rPr", ttag = "szCs", delete_ttag = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

dele_xml.qform <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "qFormat", delete_target = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

dele_xml.unhide <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "unhideWhenUsed", delete_target = TRUE)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  do.call(repl_xml.font, args)
}

reset_xml.fontCol <- function(xml,
  param = "val", as = "000000",
  ttag = "color", name = "TOC Heading",
  remove_ttagAttr = TRUE, ...)
{
  args <- as.list(environment())
  args$match.arg <- FALSE
  args <- c(args, list(...))
  do.call(repl_xml.font, args)
}

repl_xml.font <- function(xml, 
  param = c("eastAsia", "ascii", "hAnsi", "cs"),
  as = c("SimHei", "SimSun", "Times New Roman"),
  name = show_allName(xml),
  main = "style", target = "rPr", ttag = "rFonts", namespace = "w",
  match.arg = TRUE, delete_target = FALSE,
  insert_ttag = NULL, force_target = FALSE, delete_ttag = FALSE, insert_param = FALSE,
  remove_ttagAttr = FALSE)
{
  if (match.arg) {
    param <- match.arg(param)
    as <- match.arg(as)
  }
  xml[seq_along(xml)] <- XML::xmlSApply(xml,
    function(xml2) {
      if (XML::xmlName(xml2) == main) {
        xml2.name <- XML::xmlAttrs(xml2[[ "name" ]])[[ "val" ]]
        if (any(xml2.name == name)) {
          if (delete_target) {
            if (!is.null(xml2[[ target ]]))
              xml2 <- XML::removeChildren(xml2, target)
            return(xml2)
          }
          tag <- xml2[[ target ]][[ ttag ]]
          if (!is.null(tag)) {
            if (delete_ttag) {
              if (!is.null(xml2[[ target ]])) {
                xml2[[ target ]] <- XML::removeChildren(xml2[[ target ]], ttag)
              }
              return(xml2)
            } else if (remove_ttagAttr) {
              xml2[[ target ]][[ ttag ]] <- XML::removeAttributes(xml2[[ target ]][[ ttag ]])
            }
            tag.attr <- XML::xmlAttrs(tag)
            if (any(param == names(tag.attr)) | insert_param) {
              XML::xmlAttrs(xml2[[ target ]][[ ttag ]])[[ param ]] <- as
            }
          } else if (!is.null(insert_ttag)) {
            if (is.null(xml2[[ target ]]) & force_target) {
              node <- XML::xmlNode(target, namespace = namespace)
              xml2 <- XML::append.xmlNode(xml2, node)
            }
            if (!is.null(xml2[[ target ]])) {
              xml2[[ target ]] <- XML::append.xmlNode(xml2[[ target ]], insert_ttag)
            }
          }
        }
      }
      return(xml2)
    })
  return(xml)
}

rm_xml.lsdAttr <- function(xml2,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)), attrs = "qFormat")
{
  xml2[seq_along(xml2)] <- XML::xmlSApply(xml2,
    function(xml3) {
      if (any(XML::xmlAttrs(xml3)[[ "name" ]] == name)) {
        vals <- XML::xmlAttrs(xml3)
        vals <- vals[ names(vals) != attrs ]
        xml3 <- XML::removeAttributes(xml3)
        XML::xmlAttrs(xml3) <- vals
      }
      return(xml3)
    })
  return(xml2)
}

show_allName <- function(xml){
  obj <- XML::xmlSApply(xml,
    function(x) {
      x.name <- x[[ "name" ]]
      if (!is.null(x.name)) {
        if (!is.null(x[[ "rPr" ]]) | !is.null(x[[ "pPr" ]]))
          XML::xmlAttrs(x.name)[[ "val" ]]
      }
    })
  unlist(obj)
}

list_allName <- function(xml){
  lapply(seq_along(xml),
    function(n) {
      xml[[n]][[ "name"]]
    })
}

show_allTok <- function(xml){
  all <- show_allName(xml)
  c(all[ grep("Tok$", all) ])
}

if (requireNamespace("XML", quietly = TRUE)) {
  .xmlFont <- XML::xmlNode("rFonts",
    attrs = c(ascii = "Times New Roman",
      hAnsi = "Times New Roman",
      eastAsia = "SimSun",
      cs = "Nimbus Roman"), namespace = "w")
}

all_tableStyle <- function(xml){
  lst <- lapply(seq_along(xml),
    function(n) {
      attr <- XML::xmlAttrs(xml[[n]])
      if (!is.null(attr)) {
        if (any(names(attr) == "type")) {
          if (attr[[ "type" ]] == "table")
            XML::xmlAttrs(xml[[ n ]][[ "name" ]])[[ "val" ]]
        }
      }
    })
  unlist(lst)
}

# ==========================================================================
# others
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

mutate_get_fun <- function (x) {
  if (is.function(x)) {
    return(x)
  }
  if (grepl("::", x, fixed = TRUE)) {
    coumpounds <- strsplit(x, split = "::", x, fixed = TRUE)[[1]]
    z <- getFromNamespace(coumpounds[2], ns = coumpounds[1])
  } else {
    z <- getAnywhere(x)
    if (length(z$objs) < 1) {
      stop("could not find any function named ", shQuote(z$name), 
        " in loaded namespaces or in the search path. If the package is installed, specify name w
        ith `packagename::function_name`.")
    }
  }
  z
}

pdf_convert <- function (pdf, format = "png", pages = NULL, filenames = NULL, 
    dpi = 72, antialias = TRUE, opw = "", upw = "", verbose = TRUE) 
{
  config <- pdftools::poppler_config()
  if (!config$can_render || !length(config$supported_image_formats)) 
    stop("You version of libppoppler does not support rendering")
  format <- match.arg(format, config$supported_image_formats)
  if (is.null(pages)) 
    pages <- seq_len(pdftools::pdf_info(pdf, opw = opw, upw = upw)$pages)
  if (!is.numeric(pages) || !length(pages)) 
    stop("Argument 'pages' must be a one-indexed vector of page numbers")
  if (length(filenames) < 2) {
    input <- ifelse(is.raw(pdf), "output", sub(".pdf", "", 
        basename(pdf), fixed = TRUE))
  }
  if (is.null(filenames)) {
    filenames <- sprintf("%s_%d.%s", input, pages, format)
  }
  if (length(filenames) != length(pages)) 
    stop("Length of 'filenames' must be one or equal to 'pages'")
  antialiasing <- isTRUE(antialias) || isTRUE(antialias == 
    "draw")
  text_antialiasing <- isTRUE(antialias) || isTRUE(antialias == 
    "text")
  if (.Platform$OS.type == "unix") {
    sink("/dev/null")
  } else {
    sink("NUL")
  }
  res <- pdftools:::poppler_convert(pdftools:::loadfile(pdf), format, pages, filenames, 
    dpi, opw, upw, antialiasing, text_antialiasing, verbose)
  sink()
  res
}

# ==========================================================================
# other tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

as_df.lst <- function(lst, col.name = 'type', col.value = 'name') {
  data <- data.frame(
    x = rep(names(lst), lengths(lst), each = TRUE), 
    y = unlist(lst, use.names = FALSE)
  )
  colnames(data) <- c(col.name, col.value)
  data
}

cutWords <- function(ch, len = 15, len.post = 3, max = len + len.post + 1) {
  while( any(grepl(paste0("[^ ^\\-]{", max, "}"), ch)) ) {
    ch <- gsub(paste0("\\b([^ ^\\-]{", len, "})([^ ^\\-]{", len.post, ",})\\b"), "\\1-\\2", ch)
  }
  ch
}

formatBib <- function(bib = paste0(.expath, "/library.bib"), savename = "tmp.bib", 
  journalAbb_file = paste0(.expath, "/endlib.txt"))
{
  endlib <- tibble::as_tibble(data.table::fread(journalAbb_file, header = FALSE))
  dics <- .as_dic(endlib[[2]], endlib[[1]])
  names(dics) <- tolower(names(dics))
  lst <- read_bib(bib)
  lst <- lapply(lst,
    function(ch) {
      pos <- grepl("^\\s*journal =", ch)
      if (any(pos)) {
        journal <- stringr::str_extract(ch[pos], "(?<=\\{).*(?=\\})")
        journal <- tolower(journal)
        if (journal %in% names(dics)) {
          abb <- dics[[ journal ]]
          ch[pos] <- sub("\\{.*\\}", paste0("{", abb, "}"), ch[pos])
        } else {
          message("No journal abbrev: ", journal)
        }
      }
      pos <- grepl("^\\s*doi =", ch)
      ch[ !pos ]
    })
  writeLines(unlist(lst), savename)
  return(savename)
}
