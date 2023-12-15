# ==========================================================================
# render as report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

browseURL(write_articlePdf("index.Rmd", "output.Rmd", ""))

idname <- gidn()
file.copy("./output.pdf", report <- paste0(idname, ".pdf"), T)

package_results(head = NULL, masterZip = NULL, report = report)
file.rename("./client.zip", paste0(idname, ".zip"))

info <- items(
  belong = odate(12),
  coef = NA
)
