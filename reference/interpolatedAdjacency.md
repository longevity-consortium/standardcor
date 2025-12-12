# Convert a (standardized) correlation value r to an adjacency

Interpolates an adjacency value in 0..1 from a correlation value -1 \<=
r \<= 1 in a tabulated soft thresholding function.

## Usage

``` r
interpolatedAdjacency(r, adjTable)
```

## Arguments

- r:

  A numeric object containing (standardized) correlation coefficients

- adjTable:

  A tabulated soft-threshold adjacency function, as produced by
  nullModelAdjacencyTable()

## Value

A numeric object containing adjacencies for a correlation network model
