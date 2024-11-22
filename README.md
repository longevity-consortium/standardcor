# StandardCor
This package provides utilities for interpreting correlation coefficients as distances or adjacency weights for clustering.
Pairwise correlations have a geometric interpretation as the cosine of an angle between two vectors of length N, and it is
important to note that regardless how large N is, an angle is an inherently 2-dimensional concept.  The standard approaches
to convert pairwise correlations to distances retain this 2-dimensional behavior.  But when all pairs of a set of vectors is
considered, in general the vectors will be most similar in planes with different orientations, and the distribution across
all pairs will reflect the structure of N-dimensional space. This package therefore uses a Beta distribution null model
of the bulk, spurious correlations to extend the standard approaches, producing N-dimensional distances. The null model
also provides an alternative, statistical approach to defining adjacency weights as a first step in constructing a network model
of correlations for Weighted Correlation Network Analysis.

The most important contribution of this package, though, is for co-clustering datasets collected with different high-throughput ('omics)
technologies, for which the distribution of spurious correlations in general have different statistical characteristics. By fitting
each distribution separately, this package provides a means to standardize correlation values across the different technologies,
as well as the cross-correlations between data measured on different technologies. This standardization addresses differences in
the variance of spurious correlations (noise levels) across 'omics platforms, reducing the observed tendency for correlations on
noisier platforms to drive clustering more than correlations on less noisy platforms, and more than correlations across platforms.

This package therefore supports correlation-based clustering by providing a conceptual framework, methods for fitting the
distribution of spurious correlations, methods for standardizing correlation values across 'omics platforms, and new methods for converting
correlations to distances or network adjacencies.

## Installation:
`if (!require("remotes", quietly = TRUE))
   install.packages("remotes")
 remotes::install_github("longevity-consortium/standardcor")`
