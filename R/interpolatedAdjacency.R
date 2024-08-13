#' Convert a (standardized) correlation value r to an adjacency
#'
#' Interpolates an adjacency value in 0..1 from a correlation value -1 <= r <= 1
#' in a tabulated soft thresholding function.
#'
#' @param r A numeric object containing (standardized) correlation coefficients
#' @param adjTable A tabulated soft-threshold adjacency function, as produced by nullModelAdjacencyTable()
#' @return A numeric object containing adjacencies for a correlation network model
#'
#' @export
interpolatedAdjacency <- function(r, adjTable) {
  n     <- dim(adjTable)[1]
  x.min <- adjTable$x[1]
  lo <- sapply(r,
    function(x, n, x.min, x.max) {
      k <- floor(n*(x-x.min) / (x.max-x.min))
      return(ifelse(k < 1, 1, ifelse(k < n-1, k, n-1)))
    }, n-1, x.min, adjTable$x[n])
  a  <- adjTable$x[lo]
  b  <- adjTable$x[lo+1]
  return (ifelse(a < r,
          ifelse(r < b,
                 (adjTable$y[lo]*(b-r) + adjTable$y[lo+1]*(r-a)) / (b-a),
                  adjTable$y[lo+1]),
                  adjTable$y[lo]))
}


# end of interpolatedAdjacency.R
