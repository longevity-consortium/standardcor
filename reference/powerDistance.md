# Calculate power distance

Calculates distances for WCNA from correlation coefficient(s) using a
power function.

## Usage

``` r
powerDistance(r, k, unsigned = TRUE)
```

## Arguments

- r:

  Correlation coefficient(s), -1 \<= r \<= 1.

- k:

  The power parameter.

- unsigned:

  If TRUE, this function computes distances based on unsigned
  associations (\|r\|), otherwise it computes distances based on signed
  correlations (r).

## Value

A matrix of distances

## Details

Correlation and distance are transformed as follows: for type =
"unsigned"= 1 - abs(r)^k; for type = "signed" = ((1 - r)/2)^k
