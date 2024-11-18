#' Calculate Euclidean distance
#'
#' Calculates the Euclidean distances equivalent to a matrix of correlation coefficients. The Pearson correlation coefficient between two vectors is the cosine of the angle formed between the two vectors when placed with their tales at the origin. The corresponding Euclidean distance between the two vectors is the length of the vector difference between two unit vectors in the same direction as the input vectors, a value that is 0 when the two vectors point in the same direction, 1 when they are perpendicular, and 2 when they point in opposite directions.
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
