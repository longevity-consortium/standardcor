#' Replace non-zero entries in a sparse matrix with non-zero ranks
#'
#' Creates a rank matrix for a sparse matrix X using the following approach:
#' 1. Use non-zero enries in a column to calculate the ranks
#' 2. Add (z-1)/2 to the ranks (only non-zero entries are changed). z is the number of zeros
#' in the column
#' Since all the entries are shifted by the same constant (the zeros
#' are already shifted), the covariance matrix of this shifted matrix is
#' the same as the rank matrix of the entire matrix (where the zeros would
#' all also have a rank = (z+1)/2) where z is the number of zeros
#' This rank matrix can then be used to calculate Spearman's correlation coefficient
#' (via pearson correlation)
#' Copied without code changes and only minor changes to the comments, from the public blog repository
#'     https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb
#' owned by Saket Choudhary of Mumbai, India
#' 23 Mar 2024
#' @param  X A sparse data matrix
#' @param  Y A second sparse data matrix. If only X is provided, the columns of X will be correlated with each other; when Y is also specified, the columns of X are cross-correlated with the columns of Y.
#' @return A sparse matrix with ranks of non-zero entries
#' @export
SparsifiedRanks2 <- function(X) {
  X <- as(X, "sparseMatrix")
  if (class(X)[1] != "dgCMatrix") {
    X <- as(object = X, Class = "dgCMatrix")
  }
  non_zeros_per_col <- diff(x = X@p)
  n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
  offsets <- (n_zeros_per_col - 1) / 2
  x <- X@x
  ## split entries to columns
  col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
  ## calculate sparsified ranks and do shifting
  sparsified_ranks <- unlist(x = lapply(X = seq_along(col_lst),
                                        FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]))
  ## Create template rank matrix
  X.ranks <- X
  X.ranks@x <- sparsified_ranks
  return(X.ranks)
}

# end of SparsifiedRanks2.R
