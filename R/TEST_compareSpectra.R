# ## from package ProtGenerics
# setGeneric("compareSpectra", function(x, y, ...)
#     standardGeneric("compareSpectra"))
# ## from MSnbase
# setClass("Spectrum",
#          representation = representation(
#              msLevel="integer",
#              peaksCount="integer",
#              rt="numeric",
#              acquisitionNum="integer",
#              scanIndex = "integer",
#              tic = "numeric",
#              mz = "numeric",
#              intensity = "numeric",
#              fromFile = "integer",
#              centroided = "logical",
#              smoothed = "logical",
#              polarity="integer",
#              "VIRTUAL"),
#          contains=c("Versioned"),
#          prototype = prototype(
#              rt = numeric(),
#              polarity = NA_integer_,
#              acquisitionNum = NA_integer_,
#              msLevel = NA_integer_,
#              centroided = NA,
#              smoothed = NA,
#              peaksCount = 0L,
#              tic = 0,
#              scanIndex = integer(),
#              mz = numeric(),
#              intensity = numeric()),
#          validity = function(object)
#              validSpectrum(object)
#          )

# setClass("Spectrum2",
#          representation = representation(
#              merged="numeric",
#              precScanNum="integer",
#              precursorMz="numeric",
#              precursorIntensity = "numeric",
#              precursorCharge = "integer",
#              collisionEnergy = "numeric"),
#          contains=c("Spectrum"),
#          prototype = prototype(
#              merged = 1,
#              acquisitionNum = integer(),
#              precScanNum = integer(),
#              precursorMz = numeric(),
#              precursorIntensity = numeric(),
#              msLevel = as.integer(2),
#              precursorCharge = integer(),
#              collisionEnergy = numeric()),
#          validity = function(object) {
#              msg <- validMsg(NULL, NULL)
#              msl <- object@msLevel
#              if (msl < as.integer(2))
#                  msg <- validMsg(msg,
#                                  paste0("Object of class ",
#                                         class(object),
#                                         " but msLevel is ", msl,
#                                         " (should be > 1)"))
#              if (is.null(msg)) TRUE
#              else msg
#          })


# setMethod("compareSpectra", c("Spectrum", "Spectrum"),
#           function(x, y, fun=c("common", "cor", "dotproduct"),
#                    ...) {
#               compare_Spectra(x, y, fun = fun, ...)
#           })

# #' calculate similarity between spectra (between their intensity profile)
# #' @param x spectrum1 (MSnbase::Spectrum)
# #' @param y spectrum2 (MSnbase::Spectrum)
# #' @param fun similarity function (must take two spectra and ... as arguments)
# #' @param ... further arguments passed to "fun"
# #' @return double, similarity score
# #' @noRd
# compare_Spectra <- function(x, y,
#                             fun=c("common", "cor", "dotproduct"),
#                             ...) {
#   if (is.character(fun)) {
#     fun <- match.arg(fun)
#     if (fun == "cor" || fun == "dotproduct") {
#       binnedSpectra <- bin_Spectra(x, y, ...)
#       inten <- lapply(binnedSpectra, intensity)
#       return(do.call(fun, inten))
#     } else if (fun == "common") {
#       return(numberOfCommonPeaks(x, y, ...))
#     }
#   } else if (is.function(fun)) {
#     return(fun(x, y, ...))
#   }
#   return(NA)
# }


# setMethod("mz", "Spectrum", function(object) object@mz)
# setMethod("intensity", "Spectrum", function(object) object@intensity)

# bin_Spectra <- function(object1, object2, binSize = 1L,
#                         breaks = seq(floor(min(c(mz(object1), mz(object2)))),
#                                      ceiling(max(c(mz(object1), mz(object2)))),
#                                      by = binSize)) {
#     breaks <- .fix_breaks(breaks, range(mz(object1), mz(object2)))
#     list(bin_Spectrum(object1, breaks = breaks),
#          bin_Spectrum(object2, breaks = breaks))
# }

# bin_Spectrum <- function(object, binSize = 1L,
#                          breaks = seq(floor(min(mz(object))),
#                                       ceiling(max(mz(object))),
#                                       by = binSize),
#                          fun = sum,
#                          msLevel.) {
#     ## If msLevel. not missing, perform the trimming only if the msLevel
#     ## of the spectrum matches (any of) the specified msLevels.
#     if (!missing(msLevel.)) {
#         if (!(msLevel(object) %in% msLevel.))
#             return(object)
#     }
#     bins <- .bin_values(object@intensity, object@mz, binSize = binSize,
#                         breaks = breaks, fun = fun)
#     object@mz <- bins$mids
#     object@intensity <- bins$x
#     object@tic <- sum(object@intensity)
#     object@peaksCount <- length(object@mz)
#     if (validObject(object))
#         return(object)
# }

# #' The function aggregates `x` for `toBin` falling into bins defined
# #' by `breaks` using the `fun` function.
# #'
# #' @details
# #'
# #' This is a combination of the code from the former bin_Spectrum.
# #'
# #' @param x `numeric` with the values that should be binned.
# #'
# #' @param toBin `numeric`, same length than `x`, with values to be used for the
# #'     binning.
# #'
# #' @param binSize `numeric(1)` with the size of the bins.
# #'
# #' @param breaks `numeric` defining the breaks/bins.
# #'
# #' @param fun `function` to be used to aggregate values of `x` falling into the
# #'     bins defined by `breaks`.
# #'
# #' @return `list` with elements `x` and `mids` being the aggregated values
# #'     of `x` for values in `toBin` falling within each bin and the bin mid
# #'     points.
# #'
# #' @author Johannes Rainer, Sebastian Gibb
# #'
# #' @noRd
# .bin_values <- function(x, toBin, binSize = 1, breaks = seq(floor(min(toBin)),
#                                                             ceiling(max(toBin)),
#                                                             by = binSize),
#                         fun = max) {
#     if (length(x) != length(toBin))
#         stop("lengths of 'x' and 'toBin' have to match.")
#     fun <- match.fun(fun)
#     breaks <- .fix_breaks(breaks, range(toBin))
#     nbrks <- length(breaks)
#     idx <- findInterval(toBin, breaks)
#     ## Ensure that indices are within breaks.
#     idx[which(idx < 1L)] <- 1L
#     idx[which(idx >= nbrks)] <- nbrks - 1L

#     ints <- double(nbrks - 1L)
#     ints[unique(idx)] <- unlist(lapply(base::split(x, idx), fun),
#                                 use.names = FALSE)
#     list(x = ints, mids = (breaks[-nbrks] + breaks[-1L]) / 2L)
# }



# #' calculate the dot product between two vectors
# #'
# #' Stein, S. E., and Scott, D. R. (1994).
# #' Optimization and testing of mass spectral library search algorithms for
# #' compound identification.
# #' Journal of the American Society for Mass Spectrometry, 5(9), 859-866.
# #' doi: https://doi.org/10.1016/1044-0305(94)87009-8
# #'
# #' Lam, H., Deutsch, E. W., Eddes, J. S., Eng, J. K., King, N., Stein, S. E.
# #' and Aebersold, R. (2007)
# #' Development and validation of a spectral library searching method for peptide
# #' identification from MS/MS.
# #' Proteomics, 7: 655-667.
# #' doi: https://doi.org/10.1002/pmic.200600625
# #'
# #' @param x double
# #' @param y double
# #' @return double, length == 1
# #' @noRd
# dotproduct <- function(x, y) {
#   as.vector(x %*% y) / (sqrt(sum(x*x)) * sqrt(sum(y*y)))
# }



# #' Simple function to ensure that breaks (for binning) are span al leat the
# #' expected range.
# #'
# #' @param brks `numeric` with *breaks* such as calculated by `seq`.
# #'
# #' @param rng `numeric(2)` with the range of original numeric values on which
# #'     the breaks were calculated.
# #'
# #' @noRd
# .fix_breaks <- function(brks, rng) {
#     ## Assuming breaks being sorted.
#     if (brks[length(brks)] <= rng[2])
#         brks <- c(brks, max((rng[2] + 1e-6),
#                             brks[length(brks)] + mean(diff(brks))))
#     brks
# }

# ## ---------------------------------------------------------------------- 

# numberOfCommonPeaks <- function(x, y, tolerance=25e-6, relative=TRUE) {
#   sum(commonPeaks(x, y, tolerance=tolerance, relative=relative))
# }
# commonPeaks <- function(x, y, method=c("highest", "closest"),
#                          tolerance=25e-6, relative=TRUE) {
#   m <- matchPeaks(x, y, method=match.arg(method), tolerance=tolerance,
#                   relative=relative)

#   m[which(is.na(m))] <- 0L

#   as.logical(m)
# }

# matchPeaks <- function(x, y, method=c("highest", "closest", "all"),
#                        tolerance=25e-6, relative=TRUE) {
#   method <- match.arg(method)

#   if (inherits(y, "Spectrum")) {
#     y <- mz(y)
#   }

#   if (peaksCount(x) == 0 || length(y) == 0) {
#     return(integer(peaksCount(x)))
#   }

#   m <- relaxedMatch(mz(x), y, nomatch=NA, tolerance=tolerance,
#                     relative=relative)

#   if (anyDuplicated(m)) {
#     if (method == "highest") {
#       o <- order(intensity(x), decreasing=TRUE)
#     } else if (method == "closest") {
#       o <- order(abs(mz(x)-y[m]))
#     } else {
#       o <- 1:length(x)
#     }
#     sortedMatches <- m[o]

#     if (method != "all") {
#       sortedMatches[which(duplicated(sortedMatches))] <- NA
#     }
#     m[o] <- sortedMatches
#   }

#   as.integer(m)
# }

# relaxedMatch <- function(x, table, nomatch=NA_integer_, tolerance=25e-6,
#                          relative=TRUE) {

#   if (relative) {
#     if (tolerance > 1L) {
#       stop(sQuote("tolerance"),
#            " must be smaller than 1 for relative deviations.")
#     }
#     tolerance <- table*tolerance
#   }

#   match.closest(x, table, tolerance=tolerance, nomatch=nomatch)
# }

