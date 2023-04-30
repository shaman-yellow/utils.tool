# ==========================================================================
# tools for fastly write formatting thesis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inclu.fig <- function(image, land = F, saveDir = "thesis_fig", dpi = 300,
  scale = if (land) {
    list(width = 10, height = 6.2)
  } else {
    list(width = 6.2, height = 10)
  })
{
  if (!file.exists(saveDir))
    dir.create(saveDir)
  upper <- get_path(image)
  if (is.na(upper))
    upper <- "."
  file <- get_filename(image)
  ## backup for figure
  if (grepl("\\.pdf$", file)) {
    savename <- paste0(saveDir, "/", sub("\\.pdf$", ".png", file))
    if (!file.exists(savename)) {
      pdf_convert(image, filenames = savename, dpi = dpi, pages = 1)
      exists <- T
    } else exists <- F
  } else {
    savename <- paste0(saveDir, "/", file)
    if (!file.exists(savename)) {
      file.copy(image, savename)
      exists <- T
    } else exists <- F
  }
  ## trim the border
  if (!exists) {
    gc()
    pic_trim(savename)
  }
  ## record image info
  if (is.null(image_info[[ savename ]])) {
    img <- magick::image_read(savename)
    info <- magick::image_info(img)
    magick::image_destroy(img)
    image_info[[ savename ]] <- list(width = info$width, height = info$height)
    assign("image_info", image_info, envir = .env)
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
  knitr::include_graphics(savename)
}

inclu.capt <- function(img, saveDir = "thesis_fig") {
  filename <- get_filename(img)
  savename <- paste0(saveDir, "/", filename)
  if (!file.exists(savename)) {
    png(savename, 1000, 1000, res = 150)
    grid::grid.raster(png::readPNG(img))
    dev.off()
  }
  inclu.fig(savename)
}

.env <- topenv(environment())

image_info <- list()

pic_trim <- function(file){
  img <- magick::image_read(file)
  img <- magick::image_trim(img)
  magick::image_write(img, file)
  magick::image_destroy(img)
}

pretty_flex2 <- function(data,
  caption = "This is table", footer = NULL,  bold_header = F,
  weight = sc(colnames(data), rep(1, length(data))), width = 6,
  family = "Times New Roman", e.family = "SimSun", font.size = 10.5,
  caption_style = "Table Caption", form_body = F, form_header = T)
{
  args <- as.list(environment())
  do.call(pretty_flex, args)
}

pretty_flex <- function(data,
  caption = "This is table", footer = "This is footer",  bold_header = F,
  weight = sc(colnames(data), rep(1, length(data))), width = 6,
  family = "Times New Roman", e.family = "SimSun", font.size = 10.5,
  caption_style = "Table Caption", form_body = T, form_header = T)
{
  if (form_body) {
    data <- dplyr::mutate_if(
      data, function(x) ifelse(is.character(x) | is.factor(x), T, F), form
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
    }, double(1), USE.NAMES = F)
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
    }, simplify = F)
}

as_code_list <- function(names, value = rep("", length(names)), prefix = "c"){
  writeLines(paste0(prefix, "("))
  lapply(1:length(names),
    function(n) {
      end <- if (n == length(names)) "\"" else "\","
      writeLines(paste0("  ", names[[ n ]], " = \"", value[[ n ]], end))
    })
  writeLines(")")
}

as_code_list2 <- function(lst){
  lapply(1:length(lst),
    function(n) {
      writeLines(paste0("## ", names(lst[n])))
      as_code_list(rep(names(lst[n]), length(lst[[ n ]])),
        lst[[ n ]])
      writeLines("")
    })
}

# ==========================================================================
# insert text
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ftext <- officer::ftext
fp_text <- officer::fp_text
fpar <- officer::fpar
fp_par <- officer::fp_par

st.index <- fp_text(font.size = 14, font.family = "SimHei")

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
  ids.pos <- (1:length(md))[ !is.na(ids) ]
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
      for (i in 1:length(refs)) {
        md <- gsub(paste0("{@@ref:", refs[i], "}"), refs.content[i], md, fixed = T)
      }
    }
  }
  ## output
  filename <- get_filename(file)
  filepath <- get_path(file)
  writeLines(md, nfile <- paste0(filepath, "/", "_temp_", filename))
  fun_render(nfile)
}

get_chunkPos <- function(md) {
  chunk.pos <- grep("^```", md)
  unlist(lapply(seq(1, length(chunk.pos), by = 2),
      function(n) {
        chunk.pos[n]:chunk.pos[ n + 1 ]
      }))
}

fix_spell <- function(md) {
  chunk.pos <- get_chunkPos(md)
  yaml.pos <- grep("^---", md)
  yaml.pos <- 1:(yaml.pos[2])
  pos <- (1:length(md))
  pos <- pos[ !pos %in% c(chunk.pos, yaml.pos) ]
  md[ pos ] <- gsub("\"",  "\'", md[ pos ])
  md[ pos ] <- gsub("(?<=[^，^。])\\s*\'\\s*(?=[a-zA-Z])", " \'", md[ pos ], perl = T)
  md[ pos ] <- gsub("(?<=[a-zA-Z0-9])\\s*\'\\s*(?=[^，^。])", "\' ", md[ pos ], perl = T)
  md
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
        return(T)
      } else F
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
  do.call(officedown::rdocx_document, c(list('bookdown::word_document2'), args))
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
  vapply(1:length(xml2), function(n) XML::xmlName(xml2[[ n ]]), character(1))
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
  unlist(lapply(1:length(xml),
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
  args$match.arg <- F
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
    args$match.arg <- F
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
    args$match.arg <- F
    do.call(repl_xml.font, args)
  }

repl_xml.fontSpace <- function(xml,
  param = "after", as = "0", name = "caption",
  main = "style", target = "pPr", ttag = "spacing", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = F, insert_param = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

repl_xml.fontOutl <- function(xml,
  param = "val", as = 0, name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  main = "style", target = "pPr", ttag = "outlineLvl", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = F)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

repl_xml.fontAlign <- function(xml,
  param = "val", as = "center", name = "caption",
  main = "style", target = "pPr", ttag = "jc", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = F)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

repl_xml.fontSz <- function(xml,
  param = "val", as = "18", name = "heading 1",
  main = "style", target = "rPr", ttag = "sz", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = F)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

repl_xml.fontSzCs <- function(xml,
  param = "val", as = "24", name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  main = "style", target = "rPr", ttag = "szCs", namespace = "w",
  insert_ttag = XML::xmlNode(ttag, attrs = sc(param, as), namespace = namespace),
  force_target = F)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

dele_xml.lang <- function(xml,
  name = "Normal", target = "rPr", ttag = "lang", delete_ttag = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

dele_xml.ital <- function(xml,
  name = "caption", target = "rPr", ttag = "i", delete_ttag = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

dele_xml.szCs <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "rPr", ttag = "szCs", delete_ttag = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

dele_xml.qform <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "qFormat", delete_target = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

dele_xml.unhide <- function(xml,
  name = c(paste0("heading ", 1:9), paste0("toc ", 1:9)),
  target = "unhideWhenUsed", delete_target = T)
{
  args <- as.list(environment())
  args$match.arg <- F
  do.call(repl_xml.font, args)
}

reset_xml.fontCol <- function(xml,
  param = "val", as = "000000",
  ttag = "color", name = "TOC Heading",
  remove_ttagAttr = T, ...)
{
  args <- as.list(environment())
  args$match.arg <- F
  args <- c(args, list(...))
  do.call(repl_xml.font, args)
}

repl_xml.font <- function(xml, 
  param = c("eastAsia", "ascii", "hAnsi", "cs"),
  as = c("SimHei", "SimSun", "Times New Roman"),
  name = show_allName(xml),
  main = "style", target = "rPr", ttag = "rFonts", namespace = "w",
  match.arg = T, delete_target = F,
  insert_ttag = NULL, force_target = F, delete_ttag = F, insert_param = F,
  remove_ttagAttr = F)
{
  if (match.arg) {
    param <- match.arg(param)
    as <- match.arg(as)
  }
  xml[1:length(xml)] <- XML::xmlSApply(xml,
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
  xml2[1:length(xml2)] <- XML::xmlSApply(xml2,
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
  lapply(1:length(xml),
    function(n) {
      xml[[n]][[ "name"]]
    })
}

show_allTok <- function(xml){
  all <- show_allName(xml)
  c(all[ grep("Tok$", all) ])
}

.xmlFont <- XML::xmlNode("rFonts",
  attrs = c(ascii = "Times New Roman",
    hAnsi = "Times New Roman",
    eastAsia = "SimSun",
    cs = "Nimbus Roman"), namespace = "w")

all_tableStyle <- function(xml){
  lst <- lapply(1:length(xml),
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
  pdftools:::poppler_convert(pdftools:::loadfile(pdf), format, pages, filenames, 
    dpi, opw, upw, antialiasing, text_antialiasing, verbose)
}

# ==========================================================================
# other tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

as_df.lst <- function(lst, col.name = 'type', col.value = 'name') {
  data <- data.frame(
    x = rep(names(lst), lengths(lst), each = T), 
    y = unlist(lst, use.names = F)
  )
  colnames(data) <- c(col.name, col.value)
  data
}
