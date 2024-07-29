#' Calculate Dimensionally-adjusted correlation distances
#'
#' The standard (Euclidean) transformation of correlation coefficients to distances
#' is based on the 2-dimensional geometry of the plane formed by the correlated
#' vectors and the origin. However, simulations of random vectors demonstrate that
#' the probability of correlations between random vectors reaching a given correlation
#' coefficient depends on the number of dimensions spanned by the random vectors.
#' This function allows for transforming correlation coefficients to distances
#' which reflect the null distribution, which corresponds to a Beta distribution
#' that can be estimated from the bulk of the correlations under the assumption
#' that the majority of the correlations are from the null model.
#'
#' @param r A numeric object containing correlation coefficients
#' @param v The first shape parameter of a 2-parameter Beta distribution
#' @param w The second shape parameter of a 2-parameter Beta distribution
#' @param mix The relative weight of squared Euclidean distance to the Beta distance penalty. When mix = 0, only the Beta distance is used and any correlation significantly beyond the null distribution is essentially at zero distance. When mix = 1, the only correlation corresponding to zero distance is +/-1, and the Beta distance acts as a penalty on top of the squared Euclidean distance. While intermediate values are allowed, they are not recommended.
#' @param unsigned If TRUE (the default), computes distance from association (|r|) instead of correlation (r).
#' @return A numeric object containing the dimensionally-adjusted distances
#'
#' @export
betaDistance <- function(r, v=1, w=1, mix=1, unsigned=TRUE) {
  if (unsigned) {
    d = dbeta((1+abs(r))/2, v, v) / dbeta(1/2, v, v)
    d = d + mix*(1-abs(r))
  } else {
    d = pbeta((1+r)/2, v, w, lower.tail = FALSE)
    d = d + mix*(1-r)/2
  }
  return(d / (1 + mix))
}

# end of betaDistance.R
