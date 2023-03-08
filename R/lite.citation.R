# ==========================================================================
# manualy set figure number in artical
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pandoc.docx <- 
  function(
    file.md,
    path = "/mnt/data/wizard/Documents/article/",
    deal_with = T)
  {
    local.path <- getwd()
    setwd(path)
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
    system(paste0("bash pandoc.sh .TMP.md"))
    file.rename(".TMP.md.docx", paste0(file.md, ".docx"))
    setwd(local.path)
    if (deal_with)
      return(record)
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

# ==========================================================================
# modify citation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
