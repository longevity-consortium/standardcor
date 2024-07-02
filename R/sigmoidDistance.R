#' Calculate sigmoid distance
#'
#' Calculates distance network from correlation coefficient matrix using a sigmoid function.
#'
#'
#' @param mat A matrix of correlation coefficients
#' @param alpha The alpha parameter of the sigmoid function
#' @param tau0 The tau0 parameter of the sigmoid function
#' @param type The type of distance to calculate
#' @return A matrix of distances
#'
#' @export
sigmoidDistance <- function(mat, alpha, tau0, type="unsigned") {
  if (type == "unsigned") {
    s     <- 1 / (1 + exp(-alpha*(abs(r) - abs(tau0))))
    s.min <- 1 / (1 + exp(-alpha*(1 - abs(tau0))))
    s.max <- 1 / (1 + exp(-alpha*(0 - abs(tau0))))
    if (abs(s.max - s.min) > 1.0e-7) {  # Stretched to span [0..1]
      s <- (s - s.min) / (s.max - s.min)
      }
    return (s)
  } else {
    s     <- 1 / (1 + exp(-alpha*( r - tau0)))
    s.min <- 1 / (1 + exp(-alpha*( 1 - tau0)))
    s.max <- 1 / (1 + exp(-alpha*(-1 - tau0)))
    if (abs(s.max - s.min) > 1.0e-7) {  # Stretched to span [0..1]
      s <- (s - s.min) / (s.max - s.min)
      }
    return (s)
  }
}
