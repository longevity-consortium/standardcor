# Calculate sigmoid distance

Calculates distances from correlation coefficients using a sigmoid
function. WCNA converts correlations to distances and distances to
adjacencies to construct an adjacency matrix for a network model of
analyte relationships from pairwise correlations; this function would be
used in the first step. The distances can be thought of as providing a
soft thresholding of the correlations. This sigmoid function provides a
threshold parameter (tau0), which specifies the correlation value
corresponding to adjacency = 0.5; and a rate parameter (alpha \> 0),
which controls the "hardness" of the threshold.

## Usage

``` r
sigmoidDistance(r, alpha, tau0, unsigned = TRUE, stretch = FALSE)
```

## Arguments

- r:

  Correlation coefficient(s)

- alpha:

  The alpha parameter of the sigmoid function (the log of the slope at
  tau0)

- tau0:

  The tau0 parameter of the sigmoid function (the correlation value
  corresponding to distance = 0.5)

- stretch:

  If TRUE, the distances are rescaled to the range 0..1. (default:
  stretch=FALSE)

- type:

  The interpretation of the correlation coefficient signed (correlation)
  or unsigned (association, default). Any value other than "unsigned" is
  interpreted as signed.

## Value

The distance(s) corresponding to the correlation coefficient(s)
