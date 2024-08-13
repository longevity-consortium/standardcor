#' Compute soft threshold adjacency table for unsigned correlations from a symmetric Beta distribution null model
#'
#' Tabulates a soft thresholding function for converting correlations directly to adjacencies suitable for a network model
#' of a set of correlations among a set of analytes, for example a table of transcriptomic, proteomic, metabolomic, or
#' other 'Omics data. The tabulated function uses a background model for random correlations expressed with the parameter
#' nu: Pr{r | null model} = Pr{x <= Beta(nu, nu) | x = (1+r)/2}, to estimate the number of correlations that are not
#' from the background model, and then estimate for any value r the probability that, given that the value is misclassified,
#' then it should be classified as a false negative. From an estimated number of false negatives (FN) and false positives (FP),
#' this probability is FN / (FN + FP). This serves as an appropriate soft-threshold function, with a probability of 1/2 when
#' we estimate that the number of false negatives equals the number of false positives.
#'
#' @param uniqueCor A numeric object containing the full set of non-self, unique correlation values to consider. If the analyte-analyte correlation matrix is Z, the non-self, unique correlations can be found using Z[row(Z) < col(Z)].
#' @param v The parameter for the null model of random correlations, Pr{r | null model} = Pr{x <= Beta(v, v) | x = (1+r)/2}.
#' @param scale Center of the soft threshold relative to the estimated number of correlations beyond those that fit the null model. Defaults to 2; higher values increase the connectivity of the network model, lower values decrease it.
#' @return A data frame with 2 columns, x and y, tabulating y as a function of x with -1 <= x <= 1 and 0 <= y <= 1.
#'
#' @export
nullModelAdjacencyTable <- function(uniqueCor, v, scale = 2, bins=100) {

  quickSearch <- function(x,L,n) {
    lo <- 1
    hi <- n
    while (lo + 1 < hi) {
      k  <- min(max(1+lo,round((lo+hi)/2)),hi-1)
      if (x <= L[k]) {
        hi <- ifelse(k < hi, k, hi-1)
      } else {
        lo <- ifelse(lo < k, k, lo+1)
      }
    }
    r    <- ifelse(x == L[lo], lo, hi)
    r.lo <- r
    while (r.lo > 1) {
      if (L[r.lo-1] < L[r]) break
      r.lo <- r.lo - 1
    }
    r.hi <- r
    while (r.hi < n) {
      if (L[r] < L[r.hi+1]) break
      r.hi <- r.hi + 1
    }
    return(c(r.lo,r.hi))
  }

  A <- sort(abs(uniqueCor))
  n <- length(A)
  m <- round((1+sqrt(1 + 8*n))/2)
  IQR <- 2*qbeta(p=c(1,3)/4, v, v) - 1
  outside <- length(which(uniqueCor < IQR[1] | IQR[2] < uniqueCor))
  n.bg <- (n - outside)
  est.outliers <- scale * m * (1 - n.bg / n)
  Bs   <- (c(1:bins)-1/2) / bins
  adj  <- rep(1,bins)
  for (i in c(1:bins)) {
    k  <- mean(quickSearch(Bs[i], A, n))
    P  <- 1 + n - k
    k  <- round(k)
    FP <- n.bg * pbeta((1 + A[k])/2, v, v, lower.tail=FALSE)
    FN <- max(est.outliers - max(0, n - (k + FP)), 1)
    adj[i] <- FN / (FN + FP)
  }
  return(as.data.frame(cbind(x = c(-rev(Bs),Bs),y = c(rev(adj),adj))))
}

# end of nullModelAdjacencyTable.R
