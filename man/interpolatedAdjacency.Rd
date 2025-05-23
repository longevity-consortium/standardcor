% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolatedAdjacency.R
\name{interpolatedAdjacency}
\alias{interpolatedAdjacency}
\title{Convert a (standardized) correlation value r to an adjacency}
\usage{
interpolatedAdjacency(r, adjTable)
}
\arguments{
\item{r}{A numeric object containing (standardized) correlation coefficients}

\item{adjTable}{A tabulated soft-threshold adjacency function, as produced by nullModelAdjacencyTable()}
}
\value{
A numeric object containing adjacencies for a correlation network model
}
\description{
Interpolates an adjacency value in 0..1 from a correlation value -1 <= r <= 1
in a tabulated soft thresholding function.
}
