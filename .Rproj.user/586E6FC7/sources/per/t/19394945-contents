#' Calculate power distance
#'
#' Calculates distance network from correlation coefficient matrix using a power function.
#'
#' Correlation and distance are transformed as follows: for type = "unsigned"= 1 - abs(r)^k; for type = "signed" = 1 - r)/2)^k
#'
#' @param mat A matrix of correlation coefficients
#' @param k The power to raise the distance to
#' @param type The type of distance to calculate. Either "unsigned" or "signed"
#' @return A matrix of distances
#' @export
#'
#' @export
powerDistance <- function(mat, k, type="unsigned") {
  if (type == "unsigned") {
    return(1 - abs(r)^k)
  } else {
    return(((1 - r)/2)^k)
  }
}
