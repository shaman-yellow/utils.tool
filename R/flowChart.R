# ==========================================================================
# help with drawing flow chart
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sortSugi <- function(values, dec = T) {
  num <- values
  i <- 1L
  order <- order(values, decreasing = dec)
  num[order] <- vapply(num[order],
    FUN.VALUE = integer(1),
    function(value) {
      x <- i
      i <<- i + 1L
      return(x)
    })
  num
}
