
summary_week.bosai <- function(
  time = Sys.Date(),
  orders = get_orders(),
  month = lubridate::month(time),
  year = lubridate::year(time),
  week = lubridate::week(time),
  path = .prefix(),
  templ_dir = .prefix("summary/"),
  rm = F)
{
  dir <- file.path(path, paste0("summary_", week))
  targets <- c(ass = "summary.xlsx")
  if (rm) {
    unlink(list.files(dir, full.names = T, all.files = T, recursive = T), T, T)
    dir.create(dir)
    file.copy(file.path(templ_dir, targets), dir)
  }
  if (!dir.exists(dir)) {
    dir.create(dir)
    file.copy(file.path(templ_dir, targets), dir)
  }
  targets[] <- file.path(dir, targets)
  wb <- openxlsx2::wb_load(targets[[ "ass" ]])
  ## prepare orders (this week or next week)
  orders <- dplyr::filter(orders,
    is.na(lubridate::week(finish)) | (lubridate::week(finish) == !!week)
  )
  ## modify the assess table
  data_ass <- dplyr::mutate(orders, expect = finish)
  data_ass <- dplyr::select(data_ass,
    id, client, type, title, start, end, expect, finish, note
  )
  fun <- function(wb, data) {
    openxlsx2::wb_add_data(wb, 1, data, col_names = F,
      dims = do.call(openxlsx2::wb_dims, pos.data_ass), na.strings = "")
  }
  pos.data_ass <- list(3, 1)
  wb <- fun(wb, dplyr::filter(data_ass, !is.na(finish)))
  pos.data_ass <- list(15, 1)
  wb <- fun(wb, dplyr::filter(data_ass, is.na(finish)))
  openxlsx2::wb_save(wb, targets[[ "ass" ]])
  ## check
  browseURL(normalizePath(targets[[ "ass" ]]))
}

td <- function(character_date) {
  as.Date(character_date, tryFormats = "%Y%m%d")
}
