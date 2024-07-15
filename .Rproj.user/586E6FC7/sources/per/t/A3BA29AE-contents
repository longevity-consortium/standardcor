#' Calculate Euclidean distance
#'
#' Calculates distance network from correlation coefficient matrix using a Euclidean function.
#'
#'
#'
#' @param mat A matrix of correlation coefficients
#' @param v The first parameter of the beta distribution
#' @param w The second parameter of the beta distribution
#' @param type The type of distance to calculate
#' @return A matrix of distances
#'
#' @export
euclidDistance <- function(mat, v=1, w=1, type="unsigned") {
  if (type == "unsigned") {
    return((ue2d(r) + dbeta((1+abs(r))/2, v, v)) / (1 + dbeta(1/2, v, v)))
  } else {
    return((e2d(r) + pbeta((1+r)/2, v, w, lower.tail = FALSE))/2)
  }
}

e2d <- function(r) { return ((1 - r)/2) }
ue2d <- function(r) { return(1-abs(r)) }
