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
  cacheL <- list()
  model.L <- list()
  analyte.L <- list()
  nSets <- length(OmicsL)
  unset <- rnorm(1) # A value that common.samples will never be as the result of intersections

  # Build the list of analytes
  common.samples <- unset
  for (ds in names(OmicsL)) {
    if (is.matrix(OmicsL[[ds]])) {
      cacheL[[ds]] <- OmicsL[[ds]]
    } else {
      cacheL[[ds]] <- as.matrix(readRDS(OmicsL[[ds]]))
    }
    M <- cacheL[[ds]]
    samples <- unique(rownames(M))
    analytes <- unique(colnames(M))
    stopifnot(! is.null(samples))
    if (common) {
      if (unset == common.samples) {
        common.samples <- samples
      } else {
        common.samples <- intersect(common.samples,samples)
      }
    }
    stopifnot(! is.null(analytes))
    analyte.L[[ds]] <- analytes
  }
  # Analytes <- unique(unlist(analyteL))
  stopifnot((! common) | (! is.numeric(common.samples)))
  stopifnot((! common) | (min.samples <= length(common.samples)))

  for (i in c(1:nSets)) {
    # Get the ith dataset
    ds.i <- names(OmicsL)[i]
    data.i <- cacheL[[ds.i]]
    if (common) {
      samples.i <- common.samples
    } else {
      samples.i <- rownames(data.i)
    }
    analytes.i <- analyteL[[ds.i]]
    model.L[[ds.i]] <- list()

    for (j in c(i:nSets)) {
      ds.j <- names(OmicsL)[j]
      analytes.j <- analyteL[[ds.j]]

      if (j == i) {  # ---------- Self-correlation
        Zij <- SparseSpearmanCor2(data.i[samples.i, ])
        rownames(Zij) <- analytes.i
        colnames(Zij) <- analytes.i
        est.vw <- estimateShape(Zij[row(Zij) < col(Zij)], main = paste(ds.i, ds.j, sep=" x "), ...)
        if (annotate) title(paste(round(est.vw,2),collapse=", "),line=0.5)

      } else {       # ---------- Cross-correlation
        data.j <- cache[[ds.j]]
        if (common) {
          samples.j <- common.samples
        } else {
          samples.j <- rownames(data.j)
        }
        samples.ij <- intersection(samples.i,samples.j)
        stopifnot(min.samples <= length(samples.ij))
        Zij <- SparseSpearmanCor2(data.i[samples.ij, ], data.j[samples.ij,])
        rownames(Zij) <- analytes.i
        colnames(Zij) <- analytes.j
        est.vw <- estimateShape(Zij, main = paste(ds.i, ds.j, sep=" x "), ...)
        if (annotate) title(paste(round(est.vw,2),collapse=", "), line=0.5)
      }
      model.L[[ds.i]][[ds.j]] <- list(cor = Zij, shape = est.vw)
    }
  }
  return(list(modelL = model.L, analyteL = analyte.L))
}

# end of multiOmicModel.R
