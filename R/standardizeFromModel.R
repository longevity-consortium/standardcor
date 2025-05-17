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
  Analytes <- unique(unlist(analyteL, use.names = FALSE))
  N <- length(Analytes)
  Z <- matrix(0, nrow=N, ncol=N)
  rownames(Z) <- Analytes
  colnames(Z) <- Analytes
  # print(paste("Standardized matrix will be",N,"x",N))
  for (ds.i in names(modelL)) {
    # print(ds.i)
    Analytes.i <- analyteL[[ds.i]]
    stopifnot(length(intersect(Analytes,Analytes.i)) == length(Analytes.i))

    for (ds.j in names(modelL[[ds.i]])) {
      # print(paste(ds.i, ds.j))
      Analytes.j <- analyteL[[ds.j]]
      stopifnot(length(intersect(Analytes,Analytes.j)) == length(Analytes.j))

      ###
      #  dataSpec can be:
      #  -- a matrix of correlation coefficients
      #  -- the name of an RDS file containing the correlation coefficients
      #  -- the analytes do not need to be all of the entries in the matrix
      ###
      dataSpec <- modelL[[ ds.i ]][[ ds.j ]][[ 'cor' ]]
      # print(paste("cor is a",class(dataSpec),collapse=","))
      if (is.matrix(dataSpec)) {
        # print(paste("matrix",paste(dim(dataSpec),collapse=" x ")))
        Zij <- dataSpec
      } else {
        # print(paste('cor',dataSpec))
        if (! dataSpec %in% names(cacheL)) {
          cacheL[[dataSpec]] <- as.matrix(readRDS(dataSpec))
        }
        Zij <- cacheL[[dataSpec]]
      }
      stopifnot(is.matrix(Zij))
      stopifnot(length(intersect(Analytes.i, rownames(Zij))) == length(Analytes.i))
      stopifnot(length(intersect(Analytes.j, colnames(Zij))) == length(Analytes.j))
      # print(paste("From ",paste(dim(Zij),collapse=" x "),sep=" "))
      # print(paste("Using",length(Analytes.i),"x",length(Analytes.j)))

      shape <- modelL[[ds.i]][[ds.j]][['shape']]
      stopifnot(length(shape) < 3)
      # print(paste("shape",shape,collapse=" "))
      if (1 == length(shape)) {
        if ('none' == shape) { # Skip standardization
          Zc <- Zij[Analytes.i,Analytes.j]
        } else {
          Zc <- centerBeta(Zij[Analytes.i,Analytes.j], shape[1], shape[1], v.std)
        }
      } else {
        Zc <- centerBeta(Zij[Analytes.i,Analytes.j], shape[1], shape[2], v.std)
      }
      # if (! is.matrix(Zc)) {
      #   print(class(Zc))
      #   print(paste(Analytes.j,collapse=", "))
      # }
      stopifnot(is.matrix(Zc))
      rownames(Zc) <- Analytes.i
      colnames(Zc) <- Analytes.j
      Z[Analytes.i, Analytes.j] <- Zc

      if (ds.i != ds.j) {
        stopifnot(length(intersect(rownames(Z),Analytes.j)) == length(Analytes.j))
        stopifnot(length(intersect(colnames(Z),Analytes.i)) == length(Analytes.i))
        Z[Analytes.j, Analytes.i] <- t(Zc)
      }
    }
  }
  nas <- length(which(is.na(Z)))
  if (0 < nas) print(paste("WARNING: standardized matrix contains",nas,"NA entries"))

  return(Z)
}

# end of standardizeFromModel.R
