pretty_table <- 
  function(
           df,
           title = "compounds summary",
           subtitle = "LC-MS",
           footnote = "Compounds summary",
           font = "Times New Roman",
           filename = "tmp.html",
           path = ".",
           return_gt = T,
           shorter_name = T,
           spanner = F,
           group = T,
           default = F
           ){
    ## ---------------------------------------------------------------------- 
    title = paste0("**", Hmisc::capitalize(title), "**")
    subtitle = paste0("**", Hmisc::capitalize(subtitle), "**")
    colnames(df) <- Hmisc::capitalize(colnames(df))
    ## ---------------------------------------------------------------------- 
    if(!default){
      t <- base_gt_solid_line(df, title = title, subtitle = subtitle,
                              footnote = footnote, font = font, group = group)
    }
    ## ------------------------------------- 
    if(default){
      if(group){
        t <- gt(df, groupname_col = "Info")
      }else{
        t <- gt(df)
      }
      t <- opt_table_font(t, font=list(font))
      t <- tab_header(t, title = md(title), subtitle = md(subtitle))
      t <- opt_align_table_header(t, align = "left")
      t <- tab_footnote(t, footnote = footnote,
                        locations = cells_title(groups = c("title")))
    }
    ## ---------------------------------------------------------------------- 
    if(shorter_name){
      t <- cols_width(t, Name ~ px(300))
    }
    ## ---------------------------------------------------------------------- 
    if(spanner){
      columns <- colnames(df) %>% 
        .[grepl("#", .)]
      t <- tab_spanner_delim(t, columns = columns,
                             delim = "#")
    }
    ## ---------------------------------------------------------------------- 
    gtsave(t, filename, path)
    if(return_gt == T)
      return(t)
  }
## ---------------------------------------------------------------------- 
mapply_rename_col <- 
  function(
           mutate_set,
           replace_set,
           names,
           fixed = F
           ){
    envir <- environment()
    mapply(base_mapply_rename_col, mutate_set, replace_set,
           MoreArgs = list(envir = envir, fixed = fixed))
    return(names)
  }
base_mapply_rename_col <- 
  function(
           mutate,
           replace,
           envir,
           fixed = F,
           names = get("names", envir = envir)
           ){
    names <- gsub(mutate, replace, names, perl = ifelse(fixed, F, T), fixed = fixed)
    assign("names", names, envir = envir)
  }
## ---------------------------------------------------------------------- 
# add_spanner <- 
#   function(
#            t,
#            names,
#            group,
#            col_spanner
#            ){
#     envir <- environment()
#     mapply(base_add_spanner, group, col_spanner,
#            MoreArgs = list(envir = envir))
#     return(t)
#   }
# base_add_spanner <- 
#   function(
#            the_group,
#            spanner,
#            envir,
#            names = get("names", envir = envir),
#            t = get("t", envir = envir)
#            ){
#     columns <- names %>%
#       .[grepl(the_group, .)]
#     t <- t %>%
#       tab_spanner(label = spanner,
#                   columns = columns)
#     assign("t", t, envir = envir)
#   }
## ---------------------------------------------------------------------- 
base_gt_solid_line <- 
  function(
           df,
           title,
           subtitle,
           footnote,
           font,
           group
           ){
    ## ------------------------------------- 
    if(group){
      ## final group, the end row names
      tmp <- dplyr::mutate(df, row = 1:nrow(df))
      end_row <- dplyr::filter(tmp, Info == tail(unique(tmp$Info), n = 1))$row %>% 
        tail(n = 1)
      t <- gt(df, groupname_col = "Info")
    }else{
      t <- gt(df)
      end_row <- nrow(df)
    }
    ## ------------------------------------- 
    ## set font
    t <- opt_table_font(t, font=list(font))
    ## ------------------------------------- 
    ## set title or footnote etc., annotation
    t <- tab_header(t, title = md(title),
                    subtitle = md(subtitle))
    t <- tab_footnote(t, footnote = footnote,
                      locations = cells_title(groups = c("title")))
    ## ------------------------------------- 
    ## set alignment
    t <- opt_align_table_header(t, align = "left")
    t <- cols_align(t, align = "left",
                    columns = everything()) 
    t <- tab_style(t, style = cell_text(v_align = "top"),
                   locations = cells_column_labels(columns = everything())) 
    # t <- tab_style(t, style = cell_text(v_align = "top"),
                   # locations = cells_body(columns = everything())) 
    ## ------------------------------------- 
    ## set lines
    t <- opt_table_lines(t, extent = c("none"))
    ## head lines
    t <- tab_style(t, style = cell_borders(sides = c("top", "bottom"),
                                           color = "black",
                                           weight = px(1.5),
                                           style = "solid"),
                   locations = cells_column_labels()) 
    ## end lines
    t <- tab_style(t, style = cell_borders(sides = c("bottom"),
                                           color = "black",
                                           weight = px(1.5),
                                           style = "solid"),
                   locations = cells_body(columns = everything(),
                                          rows = eval(parse(text = end_row))))
    ## ------------------------------------- 
    ## set group rows
    t <- tab_style(t, style = cell_text(align = "center",
                                        weight = "bold"),
                   locations = cells_row_groups(groups = everything()))
    ## set lines
    t <- tab_style(t, style = cell_borders(sides = c("top", "bottom"),
                                           color = "grey",
                                           weight = px(1),
                                           style = "solid"),
                   locations = cells_row_groups(groups = everything()))
    return(t)
  }
