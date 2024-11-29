#' Build a standardized correlation matrix using null models for each pair of component datasets
#'
#' Uses the provided list of null models (modelL) to standardize each type of
#' correlations to the target common null model specified by the v.std parameter,
#' (r + 1)/2 ~ Beta(v.std, v.std), where r is the
#' standardized Spearman correlation coefficient between any two analytes
#' in the modeled list of 'Omics datasets. These analtyes are listed by
#' dataset in the components of the list parameter analtyeL, and must
#' be unique and contain all analytes provided in each raw correlation matrix
#' in modelL. The resulting standardized correlation matrix is built blockwise and is full symmetric.
#'
#' @param modelL A list-of-lists structure providing two inputs for each pair of 'Omics datasets A and B: a matrix modelL[[A]][[B]][['cor']] of raw Spearman correlation values, and shape parameters modelL[[A]][[B]][['shape']] == c(v,w) specifying a null model (1+r_raw)/2 ~ Beta(v,w) for the raw correlations.
#' @param analtyeL For each 'Omics dataset A, analyteL[[A]] lists the analytes provided by dataset A. These analyte identifiers are required to be unique across all datasets.
#' @param v.std 
#' @return A symmetric matrix containing standardized Spearman correlation coefficients for every pair of analytes across all datasets.
#'
#' @export
standardizeFromModel <- function(modelL, analyteL, v.std = 32) {
  cacheL <- list()
  Analytes <- unique(unlist(analyteL,use.names = FALSE))
  N <- length(Analytes)
  Z <- matrix(0, nrow=N, ncol=N)
  rownames(Z) <- Analytes
  colnames(Z) <- Analytes
  for (ds.i in names(modelL)) {
    Analytes.i <- analyteL[[ds.i]]
    for (ds.j in names(modelL[[ds.i]])) {
      Analytes.j <- analyteL[[ds.j]]
      dataSpec <- modelL[[ ds.i ]][[ ds.j ]][[ 'cor' ]]
      ###
      #  dataSpec can be:
      #  -- a matrix of correlation coefficients
      #  -- the name of an RDS file containing the correlation coefficients
      #  -- the analytes do not need to be all of the entries in the matrix
      ###
      if ('character' %in% class(dataSpec)) {
        print(paste("standardizeFromModel: file",dataSpec))
        if (! dataSpec %in% names(cacheL)) {
          cacheL[[dataSpec]] <- readRDS(dataSpec)
        }
        Zij <- cacheL[[dataSpec]]
        print(class(Zij))
        print(paste("standardizeFromModel: file",dataSpec,"size",paste(dim(Zij),collapse=" x ")))
      } else {
        Zij <- dataSpec
        print(paste("standardizeFromModel: size",paste(dim(Zij),collapse=" x ")))
      }
      shape <- modelL[[ds.i]][[ds.j]][['shape']]
      if (1 == length(shape)) {
        if ('none' == shape) { # Skip standardization
          Zc <- Zij[Analytes.i,Analytes.j]
        } else {
          Zc <- centerBeta(Zij[Analytes.i,Analytes.j], shape[1], shape[1], v.std)
        }
      } else {
        Zc <- centerBeta(Zij[Analytes.i,Analytes.j], shape[1], shape[2], v.std)
      }
      if (ds.i == ds.j) {
        Z[Analytes.i, Analytes.i] <- Zc
      } else {
        Z[Analytes.i, Analytes.j] <- Zc
        Z[Analytes.j, Analytes.i] <- t(Zc)
      }
    }
  }
  rm(cacheL)
  return(Z)
}

# end of standardizeFromModel.R
