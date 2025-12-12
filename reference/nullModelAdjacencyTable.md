# Compute soft threshold adjacency table for unsigned correlations from a symmetric Beta distribution null model

Tabulates a soft thresholding function for converting correlations
directly to adjacencies suitable for a network model of a set of
correlations among a set of analytes, for example a table of
transcriptomic, proteomic, metabolomic, or other 'Omics data. The
tabulated function uses a background model for random correlations
expressed with the parameter nu: Prr \| null model = Prx \<= Beta(nu,
nu) \| x = (1+r)/2, to estimate the number of correlations that are not
from the background model, and then estimate for any value r the
probability that, given that the value is misclassified, then it should
be classified as a false negative. From an estimated number of false
negatives (FN) and false positives (FP), this probability is FN / (FN +
FP). This serves as an appropriate soft-threshold function, with a
probability of 1/2 when we estimate that the number of false negatives
equals the number of false positives.

## Usage

``` r
nullModelAdjacencyTable(uniqueCor, v, scale = 2, bins = 100)
```

## Arguments

- uniqueCor:

  A numeric object containing the full set of non-self, unique
  correlation values to consider. If the analyte-analyte correlation
  matrix is Z, the non-self, unique correlations can be found using
  Z\[row(Z) \< col(Z)\].

- v:

  The parameter for the null model of random correlations, Prr \| null
  model = Prx \<= Beta(v, v) \| x = (1+r)/2.

- scale:

  Center of the soft threshold relative to the estimated number of
  correlations beyond those that fit the null model. Defaults to 2;
  higher values increase the connectivity of the network model, lower
  values decrease it.

## Value

A data frame with 2 columns, x and y, tabulating y as a function of x
with -1 \<= x \<= 1 and 0 \<= y \<= 1.
