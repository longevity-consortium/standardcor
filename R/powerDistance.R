#' Calculate power distance
#'
#' Calculates distances for WCNA from correlation coefficient(s) using a power function.
#'
#' Correlation and distance are transformed as follows: for type = "unsigned"= 1 - abs(r)^k; for type = "signed" = ((1 - r)/2)^k
#'
#' @param r Correlation coefficient(s), -1 <= r <= 1.
#' @param k The power parameter. 
#' @param signed If TRUE, computes distance from correlation (r) instead of association (|r|). Defaults to FALSE.
#' @return A matrix of distances
#' @export
powerDistance <- function(r, k, signed=FALSE) {
  if (signed) {
    d = ((1 - r)/2)^k
  } else {
    d = 1 - abs(r)^k
  }
  return(d)
}

# end of powerDistance.R
