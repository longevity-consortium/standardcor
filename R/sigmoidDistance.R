#' Calculate sigmoid distance
#'
#' Calculates distances from correlation coefficients using a sigmoid function. WCNA converts correlations to distances and distances to adjacencies to construct an adjacency matrix for a network model of analyte relationships from pairwise correlations; this function would be used in the first step. The distances can be thought of as providing a soft thresholding of the correlations. This sigmoid function provides a threshold parameter (tau0), which specifies the correlation value corresponding to adjacency = 0.5; and a rate parameter (alpha > 0), which controls the "hardness" of the threshold. 
#'
#' @param r Correlation coefficient(s)
#' @param alpha The alpha parameter of the sigmoid function (the log of the slope at tau0) 
#' @param tau0 The tau0 parameter of the sigmoid function (the correlation value corresponding to distance = 0.5)
#' @param signed If TRUE, computes distance from correlation (r) instead of association (|r|). Defaults to FALSE.
#' @param stretch If TRUE, the distances are rescaled to the range 0..1. (default: stretch=FALSE)
#' @return The distance(s) corresponding to the correlation coefficient(s)
#'
#' @export
sigmoidDistance <- function(r, alpha, tau0, signed=FALSE, stretch=FALSE) {
  sigmoid <- function(r,alpha,tau0) { return (1 - (1 / (1 + exp(-alpha*(r - tau0))))) }
  if (signed) {
    s <- sigmoid( r, alpha, tau0)
  } else {
    s <- sigmoid(abs(r), alpha, tau0)
  }
  if (stretch) {
    if (signed) {
      s.min <- sigmoid( 1, alpha, tau0)
      s.max <- sigmoid(-1, alpha, tau0)
    } else {
      s.min <- sigmoid( 1 , alpha, tau0)
      s.max <- sigmoid( 0 , alpha, tau0)
    }
    if (abs(s.max - s.min) > 1.0e-7) {
      s <- (s - s.min) / (s.max - s.min)
    }
  }
  return (s)
}

# end of sigmoidDistance.R
