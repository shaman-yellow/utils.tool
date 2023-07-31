# build_docu_from_script <- function(
  #   path = "~/outline",
  #   pattern = "\\.R$",
  #   output = ifelse(file.exists("mcnebula_results"), "mcnebula_results", "."),
  #   ref.docu = "mcnebula.workflow.Rmd",
  #   ref.path = system.file("extdata", ref.docu, package = "utils.tool")
  #   )
  # {
  #   ## ---------------------------------------------------------------------- 
  #   ref.cluster <- msource(script = ref.path, block = "^@.*", source = F)
  #   ref <- readLines(ref.path)
  #   ## metadata of inset
  #   df <- data.table::data.table(tag = stringr::str_extract(ref, "^@.*")) %>% 
  #     dplyr::mutate(nrow = 1:length(tag)) %>% 
  #     dplyr::filter(!is.na(tag)) %>% 
  #     dplyr::mutate(tag.name = stringr::str_extract(tag, "(?<=@).{1,}"))
  #   ## ------------------------------------- 
  #   ## build empty list
  #   merge <- list()
  #   length(merge) <- length(ref.cluster) * 2 - 1
  #   for(i in 1:length(merge)){
  #     if(i %% 2 == 1){
  #       merge[[i]] <- ref.cluster[[(i + 1) / 2]]
  #     }else{
  #       names(merge)[i] <- df$tag.name[i / 2]
  #     }
  #   }
  #   merge <- lapply(merge, function(vec){
  #     if(class(vec) == "numeric"){
  #       ref[vec]
  #     }
  #   })
  #   ## ---------------------------------------------------------------------- 
  #   script <- find_latest_script(path, pattern)$file[1]
  #   cat("[INFO]: Latest script is:", script, "\n")
  #   ## catch tag in script
  #   code <- readLines(script)
  #   tag.script <- data.table::data.table(
  #     start = grep("## \\^[a-z]{1,}_{1,}", code),
  #     end = grep("## \\$[a-z]{1,}_{1,}", code)
  #     ) %>% 
  #   dplyr::mutate(tag = stringr::str_extract(code[start], "(?<=## \\^)[a-z]{1,}(?=_)"),
  #     start = start + 1,
  #     end = end -1)
  #   ## ------------------------------------- 
  #   insert.line <-
  #     apply(tag.script, 1, simplify = F,
  #       function(vec){
  #         start <- as.numeric(vec[["start"]])
  #         end <- as.numeric(vec[["end"]])
  #         code <- code[start:end]
  #         cluster <- msource(CODE = code, source = F) %>% 
  #           lapply(
  #             function(num){
  #               code <- code[num] %>% 
  #                 .[!grepl("^## -{1,}", .)]
  #               code.intr <- code[grepl("^## \\[introduction\\] ", code)] %>% 
  #                 gsub("^## \\[introduction\\] ", "", .)
  #               ## ---------------------------------------------------------------------- 
  #               code.run <- code[!grepl("^## \\[[a-z]{1,}\\] ", code)]
  #               ## show figure of table
  #               if(T %in% grepl("^## @[a-z]{1,}", code)){
  #                 code.run <- gsub("^\\s{0,}# |^## @[a-z]{1,}", "", code.run)
  #                 return(c(code.intr, code.run))
  #               }
  #               if(length(code.run) != 0){
  #                 code.run <- c(
  #                   "",
  #                   "```{r, eval = F}",
  #                   code.run,
  #                   "```",
  #                   ""
  #                 )
  #               }
  #               code.eval <- code[grepl("^## \\[echo\\]", code)]
  #               if(length(code.eval) != 0){
  #                 code.eval <- c(
  #                   "",
  #                   "```{r, echo = T}",
  #                   code.eval,
  #                   "```",
  #                   ""
  #                 )
  #               }
  #               gather <- c(code.intr, code.run, code.eval)
  #             })
  #       })
  #   names(insert.line) <- tag.script$tag
  #   for(i in names(insert.line)){
  #     merge[[i]] <- insert.line[[i]]
  #   }
  #   md <- unlist(merge, use.names = F, recursive = T)
  #   ## ---------------------------------------------------------------------- 
  #   savename <- paste0(output, "/reports.Rmd")
  #   writeLines(md, savename)
  #   ## output pdf
  #   rmarkdown::render(savename)
  #   cat("Done\n")
  # }

