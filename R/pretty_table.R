# ==========================================================================
# format table as .tex .html ...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @import gt
pretty_table <- 
  function(
    df, title = "compounds summary", subtitle = "LC-MS",
    footnote = "...", spanner = F, default = F,
    filename = "tmp.html", path = tempdir(),
    font = "Times New Roman")
  {
    title = paste0("**", Hmisc::capitalize(title), "**")
    subtitle = paste0("**", Hmisc::capitalize(subtitle), "**")
    colnames(df) <- Hmisc::capitalize(colnames(df))
    if(!default){
      t <- gt_solid_line(df, title = title, subtitle = subtitle,
        footnote = footnote, font = font)
    }
    if(default){
      t <- opt_table_font(gt(df), font=list(font))
      t <- tab_header(t, title = md(title), subtitle = md(subtitle))
      t <- opt_align_table_header(t, align = "left")
      t <- tab_footnote(t, footnote = footnote,
        locations = cells_title(groups = c("title")))
    }
    if(spanner){
      columns <- colnames(df) %>% 
        .[grepl("#", .)]
      t <- tab_spanner_delim(t, columns = columns,
        delim = "#")
    }
    gtsave(t, filename, path)
    return(t)
  }

footnote <- function(gt, text, columns){
  tab_footnote(gt, footnote = text,
    locations = cells_column_labels(columns = !!columns))
}

gt_solid_line <- 
  function(df, title = "Table", subtitle = "Table", footnote = "...",
    font = "Times New Roman")
  {
    t <- opt_table_font(gt(df), font = list(font))
    t <- tab_header(t, title = md(title),
      subtitle = md(subtitle))
    t <- tab_footnote(t, footnote = footnote,
      locations = cells_title(groups = c("title")))
    t <- opt_align_table_header(t, align = "left")
    t <- cols_align(t, align = "left",
      columns = everything()) 
    t <- tab_style(t, style = cell_text(v_align = "top"),
      locations = cells_column_labels(columns = everything())) 
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
        rows = eval(parse(text = nrow(df)))))
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

# ==========================================================================
# for kable usage
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

kroup_row <- function(data, group = "Group", order = unique(data[[ group ]]))
{
  lst <- split(data.frame(data), data[[ group ]])
  lst <- lapply(lst,
    function(data) {
      data[[ group ]] <- c(data[[ group ]][1], rep("", nrow(data) - 1))
      data
    })
  lst <- lapply(order, function(name) lst[[ name ]])
  tibble::as_tibble(data.table::rbindlist(lst))
}


