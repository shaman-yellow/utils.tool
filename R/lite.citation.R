# ==========================================================================
# manualy set figure number in artical
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pandoc.docx <- 
  function(
    file.md,
    format = ".docx",
    output = sub("\\.[a-z]*$", format, file.md),
    deal_with = T,
    script = "pandoc.sh",
    script_path = system.file("extdata", script, package = "utils.tool"),
    template = "template.tex")
  {
    ## read md
    md <- readLines(file.md)
    if (deal_with) {
      citation <- stringr::str_extract(md, "^\\[citation\\]: .*$")
      ## get all citation
      citation <- citation[!is.na(citation)]
      citation.key <- stringr::str_extract(citation, "(?<=@).{1,20}(?=:)")
      ## unique
      citation.key <- unique(citation.key)
      ## paste as pattern.set
      pattern.set <- paste0("\\{@", citation.key, ":[^{]{1,50}\\}\\{nolink=True\\}")
      record <-  lapply(pattern.set,
        function(pattern){
          cite <- stringr::str_extract(citation, pattern)
          cite <- cite[!is.na(cite)]
          ## as table
          df <- data.table::data.table(cite = cite, seq = 1:length(cite))
          ## replace
          md <- mapply_rename_col(df$cite, df$seq, md)
          ## push to parent envir
          assign("md", md, envir = parent.frame(2))
          return(df)
        })
    }
    cat(md, sep = "\n", file = ".TMP.md")
    tmp_script <- paste0(".TMP.", script)
    file.copy(script_path, tmp_script)
    system(paste0("bash ", tmp_script, " .TMP.md ", output, " ", template))
    file.remove(tmp_script)
    if (deal_with)
      return(record)
  }

read_captions <- function(file.md, sig = "Fig. \\{@.*\\|")
{
  lst <- sep_list(readLines(file.md))
  pos <- vapply(lst, function(ch) grepl(sig, ch[1]), logical(1))
  lst <- lapply(lst[pos],
    function(ch) {
      ch <- paste0(ch, collapse = " ")
      gsub("^.* \\| ", "**", ch)
    })
  lst
}

read_bib <- function(bib) {
  lst <- sep_list(readLines(bib))
  getKey.biblst(lst)
}

shunt_bib <- function(bib, keys, export = "library.bib"){
  lst <- read_bib(bib)
  lst <- lst[names(lst) %in% keys]
  lst <- unlist(lst)
  writeLines(lst, export)
}

asTex.rmd <- function(file.md, yaml = paste0(.expath, "/", "article.yml"),
  export = paste0("tex_", sub("\\.[a-z]*$", "", get_filename(file.md))),
  bib = paste0(.expath, "/library.bib"), style = paste0(.expath, "/style.csl"))
{
  dir.create(export)
  file.copy(file.md, nfile.md <- paste0(export, "/",
      sub("\\.[a-z]*$", ".Rmd", get_filename(file.md))), T)
  md <- readLines(nfile.md)
  fig.pos <- grep("^!\\[.*\\]\\(.*\\)", md)
  fig.num <- 1
  md[fig.pos] <- vapply(md[fig.pos], FUN.VALUE = character(1),
    function(ch) {
      file <- stringr::str_extract(ch, "(?<=\\]\\().*(?=\\))")
      if (!grepl("^!\\[\\]", ch)) {
        prefix <- paste0("fig", fig.num, ".")
        fig.num <<- fig.num + 1
      } else {
        prefix <- character(1)
      }
      nfile <- gsub("fig[0-9]\\.", prefix, get_filename(file))
      whether <- file.copy(file, paste0(export, "/", nfile), T)
      if (!whether) stop("The file of figure not found")
      gsub(file, nfile, ch, fixed = T)
    })
  yml <- c("---", readLines(yaml), "---\n")
  writeLines(c(yml, md), nfile.md)
  shunt_bib(bib, extract_ref(nfile.md), paste0(export, "/", "library.bib"))
  file.copy(style, paste0(export, "/", get_filename(style)), T)
  writeLines("Done")
}

sep_list <- function(lines, sep = "^\\s*$")
{
  seps <- grep(sep, lines) + 1
  group <- 1
  groups <- vapply(1:length(lines), FUN.VALUE = double(1),
    function(n) {
      if (any(n == seps)) 
        group <<- group + 1
      group
    })
  split(lines, groups)
}

getKey.biblst <- function(lst){
  names(lst) <- vapply(lst, FUN.VALUE = character(1),
    function(lines) {
      paste0("@", stringr::str_extract(lines[1], "(?<=\\{).*(?=,$)"))
    })
  lst
}

getFilePath.biblst <- function(lst){
  lapply(lst,
    function(lines){
      n.file <- grep("^\\s*file", lines)
      if (length(n.file) == 1) {
        stringr::str_extract(lines[n.file], "(?<=\\{).*(?=\\})")
      } else {
        NULL
      }
    })
}

fix_ch2en <- function(md){
  md <- gsub("（", " (", md)
  md <- gsub("）", ") ", md)
  md <- gsub("：", ": ", md)
  md <- gsub("[\u2002]", " ", md)
  md
}

cite_extract <- function(data) {
  data <- dplyr::mutate(
    data, cite = stringr::str_extract(cite, "(?<=:)[a-zA-Z0-9_.]{1,}(?=\\})")
  )
  .as_dic(data$seq, data$cite)
}

fts_copy <- function(ft, seq.lst, path, prefix = "fig", sub_target = "fts") {
  lapply(ft, function(x){
    name <- gsub("\\..*$", "", x[2])
    file.copy(paste0(path, "/", x[1]),
      paste0(path, "/", sub_target, "/", prefix, seq.lst[[ name ]], ".", x[2]), T)
  })
}

insert_tocg <- function(md, fig) {
  pos <- grep("^## Abstract", md)
  blanks <- grep("^\\s*$", md)
  pos <- blanks[blanks > pos + 1][1]
  md[pos] <- paste0(md[pos], "\n", "\\includegraphics", "{", fig, "}", "\n")
  md
}

insert_figs <- function(md, ids, figs, captions,
  patterns = paste0("Fig\\. ", 1:length(ids)),
  at_ref = paste0("Fig. {@fig:", ids, "}"),
  figs_command = paste0("![", gsub("\n", "", captions), "]",
    "(", figs, ")", "{#fig:", ids, "}"))
{
  fb <- function(n, blanks) {
    blanks[blanks > n][1]
  }
  insert <- function(i, ch) {
    paste0(ch, "\n", figs_command[i], "\n")
  }
  rp <- function(ch, i){
    gsub(patterns[i], at_ref[i], ch)
  }
  blanks <- grep("^\\s*$", md)
  for (i in 1:length(ids)) {
    ns <- grep(patterns[i], md)
    md[ns] <- unlist(vapply(md[ns], rp, character(1), i = i))
    insert.pos <- fb(ns[1], blanks)
    md[insert.pos] <- insert(i, md[insert.pos])
  }
  return(md)
}

# ==========================================================================
# post modification (reference version dependent)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## Science bulletin

extract_ref_pan2md.sci_bull <- function(md){
  ref <- unlist(stringr::str_extract_all(md, "\\\\\\[\\^\\[.*?\\)\\\\\\]"))
  generate_ref_key(ref)
}

## Vancouver (superscript)

revise_pan2md.vanco <- function(file, output = "tmp.md"){
  md <- readLines(file)
  md <- revise_symbol_pan2md(md)
  ## ref
  md <- gsub("\\[\\^", "^[", md)
  md <- gsub("\\^\\]\\(\\\\l\\)", "](\\\\l)^", md)
  md <- gsub("\\(\\\\l\\)(?!,\\[|\\-\\-|\\^)", "(\\\\l)^", md, perl = T)
  md <- gsub("(?<!\\^|l\\),|\\-\\-)\\[(?=[0-9])", "^[", md, perl = T)
  return(md)
}

extract_ref_pan2md.vanco <- function(md) {
  ref <- unlist(stringr::str_extract_all(md, "\\^\\[.*?\\)\\^"))
  generate_ref_key(ref)
}

# ==========================================================================
# utilites for revise
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

revise_pan2md <- function(file, output = "tmp.md") {
  md <- readLines(file)
  md <- revise_symbol_pan2md(md)
}

revise_symbol_pan2md <- function(md) {
  md <- gsub("\\\\<", "&lt;", md)
  md <- gsub("\\\\>", "&gt;", md)
  md
}

## the sep is ",|--"
generate_ref_key <- function(ref, pattern = "(?<=\\[)[0-9]{1,}(?=\\])") {
  ref_num <- stringr::str_extract_all(ref, pattern)
  ref_sep <- stringr::str_extract_all(ref, ",|--")
  key <- lapply(1:length(ref_num),
    function(n) {
      ch <- rep("", length(ref_num[[n]]) * 2 - 1)
      for (i in 1:length(ch)) {
        if (i %% 2 == 1)
          ch[i] <- ref_num[[n]][(i + 1) / 2]
        else
          ch[i] <- ref_sep[[n]][i / 2]
      }
      set <- integer(0)
      i <- 1
      if (length(ch) > 1) {
        while(i < length(ch)) {
          if (ch[ i + 1 ] == ",") {
            set <- c(set, as.integer(ch[i]))
          } else if (ch[ i + 1] == "--") {
            set <- c(set, as.integer(ch[i]):as.integer(ch[i + 2]))
          }
          i <- i + 2
        }
      } else {
        set <- as.integer(ch)
      }
      set
    })
  names(key) <- ref
  key
}

extract_ref <- function(file) {
  md <- readLines(file)
  key <- stringr::str_extract_all(md, "(?<=\\[|; |^)@[0-9|a-z|A-Z|_\\-]*(?=;|\\])")
  key <- unique(unlist(key))
  key[!grepl("@fig|@tab|@.*fig$|@.*tab$", key)]
}

comb_ref <- function(n, ref) {
  com <- paste0(ref[n], collapse = "; ")
  paste0(" [", com, "]")
}

get_path <- function(path_str){
  if (!grepl("/", path_str))
    return(".")
  stringr::str_extract(path_str, ".*(?=/)")
}

get_filename <- function(path_str){
  if (!grepl("/", path_str))
    return(path_str)
  stringr::str_extract(path_str, "(?<=/)[^/]*$")
}


# ==========================================================================
# modify citation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
