# ==========================================================================
# format table as .tex .html ...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @import gt
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
    title = paste0("**", Hmisc::capitalize(title), "**")
    subtitle = paste0("**", Hmisc::capitalize(subtitle), "**")
    colnames(df) <- Hmisc::capitalize(colnames(df))
    if(!default){
      t <- gt_solid_line(df, title = title, subtitle = subtitle,
                              footnote = footnote, font = font, group = group)
    }
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
    if(shorter_name){
      t <- cols_width(t, Name ~ px(300))
    }
    if(spanner){
      columns <- colnames(df) %>% 
        .[grepl("#", .)]
      t <- tab_spanner_delim(t, columns = columns,
                             delim = "#")
    }
    gtsave(t, filename, path)
    if(return_gt == T)
      return(t)
  }

gt_solid_line <- 
  function(
           df,
           title,
           subtitle,
           footnote,
           font,
           group
           ){
    if(group){
      tmp <- dplyr::mutate(df, row = 1:nrow(df))
      end_row <- dplyr::filter(tmp, Info == tail(unique(tmp$Info), n = 1))$row %>% 
        tail(n = 1)
      t <- gt(df, groupname_col = "Info")
    }else{
      t <- gt(df)
      end_row <- nrow(df)
    }
    t <- opt_table_font(t, font=list(font))
    t <- tab_header(t, title = md(title),
                    subtitle = md(subtitle))
    t <- tab_footnote(t, footnote = footnote,
                      locations = cells_title(groups = c("title")))
    t <- opt_align_table_header(t, align = "left")
    t <- cols_align(t, align = "left",
                    columns = everything()) 
    t <- tab_style(t, style = cell_text(v_align = "top"),
                   locations = cells_column_labels(columns = everything())) 
    # t <- tab_style(t, style = cell_text(v_align = "top"),
                   # locations = cells_body(columns = everything())) 
    t <- opt_table_lines(t, extent = c("none"))
    t <- tab_style(t, style = cell_borders(sides = c("top", "bottom"),
                                           color = "black",
                                           weight = px(1.5),
                                           style = "solid"),
                   locations = cells_column_labels()) 
    t <- tab_style(t, style = cell_borders(sides = c("bottom"),
                                           color = "black",
                                           weight = px(1.5),
                                           style = "solid"),
                   locations = cells_body(columns = everything(),
                                          rows = eval(parse(text = end_row))))
    t <- tab_style(t, style = cell_text(align = "center",
                                        weight = "bold"),
                   locations = cells_row_groups(groups = everything()))
    t <- tab_style(t, style = cell_borders(sides = c("top", "bottom"),
                                           color = "grey",
                                           weight = px(1),
                                           style = "solid"),
                   locations = cells_row_groups(groups = everything()))
    return(t)
  }
