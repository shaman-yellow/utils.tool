# ==========================================================================
# render as report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write_articlePdf("index.Rmd", "output.Rmd", "")

id <- "IN2023072803-3+销售：周燕青+客户：戴心怡+斑痕增生+生信分析"
file.copy("./output.pdf", report <- paste0(id, ".pdf"), T)

package_results(head = NULL, masterZip = NULL, report = report)
file.rename("./client.zip", paste0(id, ".zip"))

info <- list(
  type = "固定业务",
  title = "...",
  status = "完成",
  coef = .25,
  date = Sys.Date(),
  receive_date = od_get_date(),
  info = od_get_info(),
  id = od_get_id(),
  score = od_get_score(),
  member = "黄礼闯"
)

