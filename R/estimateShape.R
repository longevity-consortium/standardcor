#' Estimate Beta parameters
#'
#' This function estimates the Beta parameters for a distribution given the correlations.
#'
#' @param corSet Table of correlation coefficients.
#' @param left Left parameter for the Beta distribution.
#' @param right Right parameter for the Beta distribution.
#' @param plot If TRUE, a plot of the distribution is shown.
#' @param fine Number of bins for the histogram.
#' @param trim Proportion of the data to be trimmed.
#' @return A vector with the two Beta parameters.
#'
#' @export
estimateBetas <- function(corSet, left = 0, right = 0, plot=FALSE, fine = NULL, trim=0.01, ...) {
  cors <- corSet[! is.na(corSet)]
  n    <- length(cors)
  good <- abs(cors) < n/(1+n) # exclude extreme values near +/-1
  if (length(which(good)) < 5) { good <- rep(TRUE,n) }
  n.good <- length(which(good))

  # Method of moments
  x    <- (1+cors)/2 # compress range from (-1,1) to [0,1]
  mu <- max(1.0e-7, min(mean(x[good]), 1-1.0e-7))
  s2 <- mad(x[good])^2
  z1 <- mu*(1-mu)/s2 - 1

  # Untrimmed method of moments estimate
  v1 <- mu*z1
  w1 <- (1-mu)*z1

  # Estimating the mode
  obs   <- sort(unique(round(x[good],3)))
  n.obs <- length(obs)
  width <- round(min(max(5,min(sqrt(n.obs),100))))
  # println("evw: ", width, n.obs, length(x[good]))
  obs.density <- unlist(sapply(c(1:(1+n.obs-width)),
                               density.in,width,n.obs,obs,x[good]))
  # println("evw: ",length(obs.density), length(which(is.na(obs.density))))
  k <- which.max(obs.density)
  mode.density <- obs.density[k]
  mode.region <- c(k:(k+width-1))
  mode <- ( sum(obs.density[mode.region]*obs[mode.region])
            / sum(obs.density[mode.region]))

  # Estimating z = v+w
  z2 <- max(2*mode.density^2,3) # Heuristic starting point
  for (tries in c(1:10)) {
    v <- 1 + max(mode*(z2-2),0)
    zest.dens <- dbeta(mode, v, max(z2 - v, 1))
    f <- zest.dens / mode.density
    if (abs(f-1) < 0.01) break
    z2 <- max(z2/f,2)
    if (10 == tries) { # optimization failed
      println("evw: warning suboptimal",round(100*f,2),
              "%, v+w",round(z1,2),",",round(z2,2))
    }
  }
  # Mode-based estimate (unchanged by trimming)
  v2 <- 1 + mode * (z2 - 2)
  w2 <- z2 - v

  # Determine trimming limits
  v  <- (v1 + v2) / 2
  w  <- (w1 + w2) / 2
  dig  <- 1
  qq   <- qbeta(c(trim/2,1-trim/2), v, w)
  good <- (qq[1] < x) & (x < qq[2])
  if (length(which(good)) < 5) { good <- rep(TRUE,length(x)) }

  # Method of moments on trimmed distribution
  mu <- max(1.0e-7, min(mean(x[good]), 1-1.0e-7))
  s2 <- mad(x[good])^2
  z3 <- mu*(1-mu)/s2 - 1
  v3   <-   mu  * z3 + right
  w3   <- (1-mu)* z3 + left

  # Final estimate
  v <- (v2+v3)/2
  w <- (w2+w3)/2

  if (plot) {
    if (is.null(fine)) {
      fine <- min(max(5,round(sqrt(length(corSet)))),100)
    }
    Bs   <- (c(-fine:(1+fine))-0.5)/fine
    r    <- c(-1,Bs[abs(Bs) < 1],1)
    H    <- hist(corSet, breaks=Bs, prob=TRUE,
                 xlab="Correlation", plot=plot, ...)
    box()
    abline(h=0)
    abline(v=0)
    abline(v=2*qq-1,lty=2,col='gray')
    lines(r, dbeta((1+r)/2,v,w)/2, lwd=3,col='MediumBlue')
  }
  return(c(v,w))
}
# Internal function
#' @export
density.in <- function(k,width,n,obs,data) {
  kk     <- k+width-1
  left  <- ifelse(k  > 1, (obs[k-1] + obs[k])   /2, 0)
  right <- ifelse(kk < n, (obs[kk]  + obs[kk+1])/2, 1)
  count  <- length(which(left <= data & data <= right))
  denom  <- length(data)*(right - left)
  stopifnot((! is.na(denom)) | (count > 0) | (! is.na(denom)) | (denom > 0))
  return(count/denom)
}
