#' Convert a (standardized) correlation value r to an adjacency
#'
#' Interpolates an adjacency value in 0..1 from a correlation value -1 <= r <= 1
#' in a tabulated soft thresholding function.
#'
#' @param r A numeric object containing (standardized) correlation coefficients
#' @param adjTable A tabulated soft-threshold adjacency function, as produced by nullModelAdjacencyTable()
#' @return A numeric object containing the adjacencies for a correlation network
#'
#' @export
interpolatedAdjacency <- function(r, adjTable) {
  n     <- dim(adjTable)[1]
  x.min <- adjTable$x[1]
  lo <- min(max(floor((r-x.min) / (adjTable$x[n]-x.min)), 1), n-1)
  a  <- adjTable$x[lo]
  b  <- adjTable$x[lo+1]
  return (ifelse(a < x,
          ifelse(x < b,
                 (adjTable$y[lo]*(b-x) + adjTable$y[lo+1]*(x-a)) / (b-a),
                  adjTable$y[lo+1]),
                  adjTable$y[lo]))
}

# end of interpolatedAdjacency.R
