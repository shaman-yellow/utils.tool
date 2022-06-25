generate_slidy <- 
  function(
           name,
           path = "/mnt/data/wizard/Documents/zcmu_reports"
           ){
    ## ------------------------------------- 
    file <- paste0(path, "/", name, ".Rmd")
    meta <- generate_slidy_meta()
    cat(meta, file = file)
  }
generate_slidy_meta <- 
  function(
           ...
           ){
    main <- 
"---
title: 'title'
author: 'Lichuang Huang'
institute: 'ZCMU'
date: \"`r format(Sys.Date())`\"
output:
  beamer_presentation:
    toc: true
    df_print: tibble
    theme: 'Dresden'
    colortheme: 'beaver'
    fonttheme: 'structurebold'
fontsize: 10pt
---"
return(main)
  }
make_slidy <- 
  function(
           name,
           path = "/mnt/data/wizard/Documents/zcmu_reports"
           ){
    file <- paste0(path, "/", name, ".Rmd")
    rmarkdown::render(file)
  }
