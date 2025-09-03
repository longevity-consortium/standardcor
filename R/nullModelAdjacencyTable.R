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
#' @param scale Center of the soft threshold relative to the estimated number of correlations beyond those that fit the null model. Defaults to 1; higher values increase the connectivity of the network model, lower values decrease it.
#' @return A data frame with 2 columns, x and y, tabulating y as a function of x with -1 <= x <= 1 and 0 <= y <= 1.
#'
#' @export
nullModelAdjacencyTable <- function(uniqueCor, v, scale = 1, signed=FALSE, bins=100) {

  # return the position of x in the sorted array L
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

  ###
  # Here is the logic for estimating the number of outliers (n.outliers) versus
  # samples from the null model (n.bg):
  #    n.obs   = n.bg     + n.outliers
  #    outside = n.bg / 2 + n.outliers
  # Therefore,
  #    n.outliers = 2 (n.bg / 2 + n.outliers) - (n.bg + n.outliers)
  #               = 2*outside - n.obs
  ###
  A       <- sort(abs(uniqueCor))
  n.obs   <- length(A)
  IQR     <- 2*qbeta(p=c(1,3)/4, v, v) - 1
  outside <- length(which(uniqueCor < IQR[1] | IQR[2] < uniqueCor))
  est.outliers <- scale * (2*outside - n.obs)
  n.bg <- n.obs - est.outliers

  if (signed) { # Asymmetric, -1 < Bs[i] < 1
    Bs   <- (c((1-bins):bins)-1/2) / bins
    adj  <- rep(1,length(Bs))
    for (i in c(1:length(Bs))) {
      k <- round(mean(quickSearch(Bs[i], A, n.obs)))
      FP <- n.bg * pbeta((1 + A[k])/2, v, v, lower.tail=FALSE)
      outside.high <- length(which(IQR[2] < uniqueCor))
      est.high <- scale * (outside.high - n.bg / 4)
      FN <- max(est.high - max(0, (n.obs-k) - FP), 1) # observations above Bs[i] - FP
      adj[i] <- FN / (FN + FP)
    }
    AdjTable <- as.data.frame(cbind(x = Bs,y = adj)) # Asymmetric

  } else { # Symmetric, mirroring results for 0 < Bs[i] < 1
    Bs   <- (c(1:bins)-1/2) / bins
    adj  <- rep(1,length(Bs))
    for (i in c(1:length(Bs))) {
      k <- round(mean(quickSearch(Bs[i], A, n.obs)))
      FP <- n.bg * pbeta((1 + A[k])/2, v, v, lower.tail=FALSE)
      FN <- max(est.outliers - max(0, (n.obs-k) - FP), 1) # observations above Bs[i] - FP
      adj[i] <- FN / (FN + FP)
    }
    AdjTable <- as.data.frame(cbind(x = c(-rev(Bs),Bs),y = c(rev(adj),adj))) # Symmetric
  }
  return(AdjTable)
}

# end of nullModelAdjacencyTable.R
