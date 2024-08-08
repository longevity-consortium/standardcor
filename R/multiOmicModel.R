#' Compute Spearman correlations and normalization distributions within and between multiple datasets
#'
#' Manages computation of Spearman correlations within and between multiple, potentially large ('Omics) datasets,
#' returning a matrix of correlations between pairs of analytes. Each listed dataset is expected to be a
#' numeric matrix containing only numeric data, however the input list can contain either the matrix itself
#' (referred to here as internal datasets) or for large datasets, the name of an RDS file containing the
#' dataset (external datasets). This approach cannot handle an unlimited amount of data; at the very least,
#' all internal datasets and every pair of external datasets must fit into memory simultaneously.
#' This function computes the null models, and has been separated from the block-wise construction of the
#' standardized correlation matrix between all analytes to provide the opportunity for the null model
#' parameters computed here to be adjusted as needed prior to use for standardization.
#'
#' @param OmicsL A named list of datasets. Each dataset may be a matrix (internal dataset) or the name of an RDS file containing a matrix (external dataset). The matrices are interpreted as rows representing samples and columns representing analytes; the analyte names must be unique across all datasets.
#' @param common When TRUE (the default), only samples common to all datasets will be used to compute correlations. Otherwise, correlations will be computed across all samples shared by each individual pair of datasets.
#' @param min.samples A lower limit on the number of samples from which analyte-analyte correlations may be computed. This parameter is used to ensure that datasets overlap sufficiently on samples to provide at least a minimal ability to compute correlations between analytes. When a common set of samples is required across all datasets, this parameter is the minimum length of the common set of samples; otherwise the min.samples is the minimum number of samples shared by each pair of 'Omics datasets independently.
#' @param common When TRUE (the default), only samples common to all datasets will be used to compute correlations. Otherwise, correlations will be computed across all samples shared by each individual pair of datasets. While requiring all correlations to depend on a common set of samples, setting common to FALSE permits all available data to contribute to each cross-correlation; it is up to the user to decide when this is appropriate.
#' @return A list of two parts: 'modelL', a list of Beta distributon null models and raw correlation values for every pair of input ('Omics) datasets, and 'analyteL', a list of the analytes contributed by each 'Omics dataset. The standardized correlation matrix can be computed from these two lists using standardizeFromModel(), either directly or after adjustment of the individual null model parameters as needed. The rows and columns of the final standardized correlation matrix will be labeled by the full set of analytes listed in analyteL, which must be unique; no effort is made to render the analyte names unique in either of these functions.
#'
#' @export
multiOmicModel <- function(OmicsL, common = TRUE, min.samples = 5, annotate=FALSE, ...) {
  sampleNames <- function(M) {
    external <-
    if (! is.matrix(M)) {
      M.in <- readRDS(M)
      names <- rownames(M.in)
      rm(M.in)
    } else {
      names <- rownames(M)
    }
    return(names)
  }
  model.L <- list()
  nSets <- length(OmicsL)
  if (common) { # Determine the (sub)set of samples common to all datasets
    common.samples <- sampleNames(OmicsL[[1]])
    for (i in c(1:nSets)) {
      common.samples <- intersect(common.samples, sampleNames(OmicsL[[i]]))
    }
    stopifnot(min.samples <= length(common.samples))
    gc("no")
  }

  analyte.L <- list()
  for (i in c(1:nSets)) {
    # Get the ith dataset
    ds.i <- names(OmicsL)[i]
    model.L[[ds.i]] <- list()
    external.i <- ! is.matrix(OmicsL[[ds.i]])
    if (external.i) {
      data.i <- readRDS(OmicsL[[ds.i]]) # need to try...catch this
    } else { # dataset is in-memory
      data.i <- OmicsL[[ds.i]]
    }
    if (common) {
      samples.i <- common.samples
    } else {
      samples.i <- rownames(data.i)
    }
    stopifnot(min.samples <= length(samples.i))

    for (j in c(i:nSets)) {
      ds.j <- names(OmicsL)[j]
      model.L[[ds.i]][[ds.j]] <- list()
      if (j == i) {  # ---------- Self-correlation
        if (1 == i) analyte.L[[ds.i]] <- colnames(data.i)
        Z <- SparseSpearmanCor2(data.i[samples.i, ])
        shape <- estimateShape(Z[row(Z) < col(Z)], main = paste(ds.i, ds.j, sep=" x "), ...)
        if (annotate) title(paste(round(shape,1),collapse=", "),line=0.5)
      } else {       # ---------- Cross-correlation
        # Get the jth dataset
        external.j <- ! is.matrix(OmicsL[[ds.j]])
        if (external.j) {
          data.j <- readRDS(OmicsL[[ds.j]]) # need to try...catch this
        } else { # internal dataset
          data.j <- OmicsL[[ds.j]]
        }
        if (1 == i) analyte.L[[ds.j]] <- colnames(data.j)
        if (common) {
          samples.ij <- common.samples
        } else {
          samples.ij <- intersect(samples.i, rownames(data.j))
        }
        stopifnot(min.samples <= length(samples.ij))
        Z <- SparseSpearmanCor2(data.i[samples.ij, ], data.j[samples.ij,])
        shape <- estimateShape(Z, main = paste(ds.i, ds.j, sep=" x "), ...)
        if (annotate) title(paste(round(shape,1),collapse=", "), line=0.5)
        if (external.j) { # external dataset
          rm(data.j)
          gc("no")
        }
      }
      model.L[[ds.i]][[ds.j]][['cor']]   <- Z
      model.L[[ds.i]][[ds.j]][['shape']] <- shape # left = 0, right = 0
    }
    if (external.i) { # external dataset
      rm(data.i)
      gc("no")
    }
  }
  return(list(modelL = model.L, analyteL = analyte.L))
}

# end of multiOmicModel.R
