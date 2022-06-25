build_docu_from_script <- 
  function(
           path = "~/outline",
           pattern = "\\.R$",
           ref.docu = "mcnebula.workflow.Rmd",
           ref.path = system.file("extdata", ref.docu, package = "utils.tool")
           ){
    script <- find_latest_script(path, pattern)$file[1]
    cat("[INFO]: Latest script is:", script, "\n")
    ## ---------------------------------------------------------------------- 
    cluster <- msource(path = path, pattern = pattern,
                       script = script, source = F)
    return(cluster)
  }
