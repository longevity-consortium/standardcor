#' Sparse implementation of Spearman Correlation
#'
#' Computes Spearman correlations either among the columns of one sparse matrix (X) or between the columns of two sparse matrices (X, Y). This implementation is specific to Spearman correlation and highly efficient; the efficiency results from avoiding unnecessary computations involving the ranks of zero entries, as well as the equivalence of rank correlations when the ranks are shifted by a constant.
#' Copied without code changes, and some textual change to the comments, from the public blog repository
#'     https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb
#' owned by Saket Choudhary of Mumbai, India
#' 23 Mar 2024
#' @param  X A sparse data matrix
#' @param  Y A second sparse data matrix. If only X is provided, the columns of X will be correlated with each other; when Y is also specified, the columns of X are cross-correlated with the columns of Y.
#' @param cov An optional covariance matrix to use in computing the correlations. Passed to qlcMatrix::corSparse().
#' @return A matrix M of Spearman correlations. When is.null(Y), M[i,j] = cor(X[,i],X[,j],method='s'), otherwise M[i,j] = cor(X[,i],Y[j],method='s'). This implementation is considerably faster than others when X and Y are sparse.
#' @export
SparseSpearmanCor2 <- function(X, Y = NULL, cov = FALSE) {

  # Get sparsified ranks
  rankX <- SparsifiedRanks2(X)
  if (is.null(Y)){
    # Calculate pearson correlation on rank matrices
    return (corSparse(X=rankX, cov=cov))
  }
  rankY <- SparsifiedRanks2(Y)
  return(corSparse( X = rankX, Y = rankY, cov = cov))
}

# end of SparseSpearmanCor2.R
