## R
pandoc.docx <- 
  function(
           file.md,
           path = "~/Documents/article/"
           ){
    local.path <- getwd()
    setwd(path)
    ## read md
    md <- readLines(file.md)
    citation <- stringr::str_extract(md, "^\\[citation\\]: .*$")
    ## get all citation
    citation <- citation[!is.na(citation)]
    ## ------------------------------------- 
    citation.key <- stringr::str_extract(citation, "(?<=@).{1,20}(?=:)")
    ## unique
    citation.key <- unique(citation.key)
    ## paste as pattern.set
    pattern.set <- paste0("\\{@", citation.key, ":[^{]{1,50}\\}\\{nolink=True\\}")
    ## ------------------------------------- 
    record <-  lapply(pattern.set, function(pattern){
                        cite <- stringr::str_extract(citation, pattern)
                        cite <- cite[!is.na(cite)]
                        ## as table
                        df <- data.table::data.table(cite = cite, seq = 1:length(cite))
                        ## replace
                        md <- extra::mapply_rename_col(df$cite, df$seq, md)
                        ## push to parent envir
                        assign("md", md, envir = parent.frame(2))
                        return(df)
           })
    cat(md, sep = "\n", file = ".TMP.md")
    ## ------------------------------------- 
    system(paste0("bash pandoc.sh .TMP.md"))
    file.rename(".TMP.md.docx", paste0(file.md, ".docx"))
    ## ------------------------------------- 
    setwd(local.path)
    return(record)
  }
