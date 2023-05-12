# ==========================================================================
# plot or ...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

bibnetwork <- function(db, save, topic = "keywords", n = 30,
  labelsize = 1.2, seed = 20, cluster = "none", width = 12, height = 7,
  title = "", ...)
{
  NetMatrix <- biblioNetwork(
    db, analysis = "co-occurrences", network = topic, sep = ";"
  )
  pdf(save, width, height)
  set.seed(seed)
  net <- networkPlot(
    NetMatrix, normalize = "association",
    weighted = T, n = n, Title = title,
    type = "fruchterman", size = T, edgesize = 5, labelsize = labelsize,
    cluster = cluster, ...
  )
  dev.off()
  return(net)
}
