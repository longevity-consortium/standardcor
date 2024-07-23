#' Center correlations
#'
#'
#' The null distribution model is r ~ 2 Beta(v,w) - 1
#' The centering (target) distribution is r_centered ~ 2 Beta(nu, nu) - 1
#'
#' @param r The correlation to center
#' @param v The first parameter of the Beta distribution
#' @param w The second parameter of the Beta distribution
#' @param nu The parameter of the target Beta distribution
#' @return The centered correlation
#' @export
centerBeta <- function(r, v, w, nu) {
  return(2*qbeta(pbeta((1 + r)/2, v, w), nu, nu) - 1)
}
