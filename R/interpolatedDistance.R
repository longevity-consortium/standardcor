#' Convert a (standardized) correlation value r to a model-based distance
#'
#' Interpolates an adjacency value in 0..1 from a correlation value -1 <= r <= 1
#' in a tabulated soft thresholding function, and averages 1-adjacency with squared Euclidean distance.
#' This is an alternative to the "softer" method of adding the probability of the correlation under
#' the null model to squared Euclidean distance; this version uses the FN/(FN+FP) logic used in
#' the WCNA adjacency function.
#'
#' @param r A numeric object containing (standardized) correlation coefficients
#' @param adjTable A tabulated soft-threshold adjacency function, as produced by nullModelAdjacencyTable()
#' @param signed When TRUE, use a signed interpretation of correlations; default FALSE (unsigned).
#' @return A numeric object containing adjacencies for a correlation network model
#'
#' @export
interpolatedDistance <- function(r, adjTable, signed=FALSE) {
  adj <- interpolatedAdjacency(r, adjTable)
  return ((squaredEuclidean(r, signed) + (1-adj))/2)
}

# end of interpolatedDistance.R
