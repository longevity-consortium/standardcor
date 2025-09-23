#' Squared Euclidean Distances
#'
#' Returns the signed or unsigned squared Euclidean distance equivalent to (a matrix of) correlation coefficients,
#' scaled to the range 0 to 1.
#' 
#' The Pearson correlation coefficient of a paired sample { (x_i, y_i) } is defined geometrically
#' in the 2D plane defined in n dimensions by the origin, the point x = <x_i>, and the point y = <y_i>.
#' Using 1 = <1,0> to represent the unit vector in the x-direction in the plane, let the y-direction
#' in the plane be represented by unit vector u. The correlation coefficient r = r(x,y) = r({(x_i, y_i)})
#' is the cosine of the angle from 1 to u.
#' 
#' The corresponding Euclidean distance is the length of u - 1, but it turns out that the square of the
#' Euclidean distance is simply 1 - r. Note that while a correlation coefficient may be positive or
#' negative, both positive and negative correlation indicate that the two values in the pair are
#' associated, in the sense that knowing one can help predict the other. We can therefore use the
#' (signed) correlation coefficient or its absolute value, a measure of (unsigned) association, to
#' convert correlation to distance. The range of (squared) signed correlation distances defined above
#' is [0, 2], but this function rescales the signed version to have the same [0,1] range as the
#' (squared) unsigned version.
#' 
#' @param r A matrix of correlation coefficients
#' @param type The type of distance to calculate: 'unsigned' (Default) or 'signed'. The unsigned squared distance of independent samples is 1, and perfectly correlated and anticorrelated samples reach distance 0. For signed squared distances, perfectly anticorrelated samples have distance 1, uncorrelated samples have distance 1/2, and only perfectly correlated samples reach distance 0.
#' @return A numeric (matrix, array, or scalar, depending on r) of distances
#'
#' @export
squaredEuclidean <- function(r, signed-FALSE) {
  if (signed) {
    d <- (1 - r) / 2 # naturally 0 to 2, scaled down to 0 to 1
  } else {
    d <- 1 - abs(r)  # naturally 0 to 1
  }
  return(d)
}

# end of squaredEuclidean.R
