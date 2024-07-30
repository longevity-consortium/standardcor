#' Calculate power distance
#'
#' Calculates distances for WCNA from correlation coefficient(s) using a power function.
#'
#' Correlation and distance are transformed as follows: for type = "unsigned"= 1 - abs(r)^k; for type = "signed" = ((1 - r)/2)^k
#'
#' @param r Correlation coefficient(s), -1 <= r <= 1.
#' @param k The power parameter. 
#' @param unsigned If TRUE, this function computes distances based on unsigned associations (|r|), otherwise it computes distances based on signed correlations (r).
#' @return A matrix of distances
#' @export
powerDistance <- function(r, k, unsigned=TRUE) {
  if (unsigned) {
    d = 1 - abs(r)^k
  } else {
    d = ((1 - r)/2)^k
  }
  return(d)
}

# end of powerDistance.R
