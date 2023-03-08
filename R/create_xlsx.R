# ==========================================================================
# create .xlsx via openxlsx2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

xl_dim <- function(rows, cols) {
  if (max(cols) > 26) {
    stop( "max(cols) > 26" )
  }
  paste0(LETTERS[min(cols)], min(rows),
    ":", LETTERS[max(cols)], max(rows))
}

xl_table <- function(
  data, title = "Table",
  group_by = "Group",
  font = "Times New Roman",
  wb = openxlsx2::wb_workbook()) 
{
  wb$add_worksheet()
  ## title
  wb$merge_cells(rows = 1, cols = seq_along(data))
  wb$add_data(x = title)
  title_dim <- xl_dim(1, seq_along(data))
  wb$add_font(, title_dim, font, bold = "double")
  ## data
  wb$add_data_table(x = data, startRow = 2, withFilter = F, na.strings = "")
  data_dim <- xl_dim(2:(nrow(data) + 2), 1:ncol(data))
  wb$add_font(, data_dim, font)
  wb$add_cell_style(, data_dim, horizontal = "left", vertical = "top")
  ## group
  group_col <- grep(group_by, colnames(data))
  if (length(group_col) != 0) {
    group <- split(1:nrow(data) + 1, data[[ group_by ]])
    for (i in group) {
      wb$merge_cells(rows = i + 1, cols = group_col)
    }
  }
  ## width
  nchar <- rbind(nchar(colnames(data)), apply(data, 2, nchar))
  nchar.max <- apply(nchar, 2, function(x) max(x, na.rm = T))
  nchar.max <- vapply(nchar.max, function(x) if (x > 30) 30 else x, numeric(1))
  for (i in 1:ncol(data)) {
    wb$set_col_widths(, cols = i, width = nchar.max[i] * 1 + 3)
  }
  ## border
  header_dim <- xl_dim(2, seq_along(data))
  wb$add_border( , header_dim,
    left_border = NULL, right_border = NULL,
    top_border = "double", bottom_border = "double"
  )
  end_dim <- xl_dim(nrow(data) + 2, seq_along(data))
  wb$add_border(, end_dim,
    left_border = NULL, right_border = NULL,
    top_border = NULL, bottom_border = "double"
  )
  return(wb)
}

xl_save <- function(wb, path) {
  openxlsx2::wb_save(wb, path)
}
