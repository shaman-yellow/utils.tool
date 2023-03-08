# ==========================================================================
# tools for fastly write formatting thesis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

inclu.fig <- function(image, saveDir = "thesis_fig", dpi = 300){
  if (!file.exists(saveDir))
    dir.create(saveDir)
  upper <- get_path(image)
  if (is.na(upper))
    upper <- "."
  file <- get_filename(image)
  if (grepl("\\.pdf$", file)) {
    savename <- paste0(saveDir, "/", sub("\\.pdf$", ".png", file))
    if (!file.exists(savename)) {
      pdf_convert(image, filenames = savename, dpi = dpi)
      exists <- T
    } else exists <- F
  } else {
    savename <- paste0(saveDir, "/", file)
    if (!file.exists(savename)) {
      file.copy(file, savename)
      exists <- T
    } else exists <- F
  }
  if (!exists) {
    pic_trim(savename)
  }
  knitr::include_graphics(savename)
}

pic_trim <- function(file){
  img <- magick::image_read(file)
  img <- magick::image_trim(img)
  magick::image_write(img, file)
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
# custom docx tools
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

custom_docx_document <- function(...){
  args <- list(...)
  args <- list(
    reference_docx = args$reference_docx,
    df_print = "kable",
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
