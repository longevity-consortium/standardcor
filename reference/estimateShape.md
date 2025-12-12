# Estimate parameters of a Beta distribution null model

This function estimates the shape parameters v and w for a Beta(v,w)
distribution to fit the bulk (background) distribution of a given set of
correlations. The objective of the estimating process is to match the
density of the distribution at the mode and the breadth of the
distribution around the mode, with the goal of estimating a model of the
bulk (null) distribution under the assumption that the correlation
coefficents provided as corSet contains a modest number of values not
drawn from the null distribution. Note that it is not possible to
provide an infallible definition of the desired null model; this method
is therefore necessarily heuristic and will not provide a satisfactory
result in all cases. When this occurs, the user is encouraged to use the
"left" and "right" parameters as needed to provide a satisfactory null
model.

## Usage

``` r
estimateShape(
  corSet,
  left = 0,
  right = 0,
  plot = FALSE,
  fine = NULL,
  trim = 0.01,
  ...
)
```

## Arguments

- corSet:

  A vector of non-unique, non-self correlation coefficients.

- left:

  An adjustment to the estimated value of w. This parameter is provided
  so the effect of changing the estimated v parameter can be assessed in
  the plot. Note that the mean of Beta(v,w) is v/(v+w), adding to w
  moves the mean to the left.

- right:

  An adjustment to the estimated value of v. Adding to the value of v
  moves the mean of the Beta distribution to the right; Adding to both v
  and w decreases the variance.

- plot:

  If TRUE, a plot of the distribution is shown, along with the fitted
  distribution Beta(v,w). If ev and ew are the estimated parameters, v =
  ev + right, w = ew + left.

- fine:

  Half the number of bins for the histogram. When fine = NULL (the
  default), an appropriate value between 5 and 100 is used.

- trim:

  The mean of the correlations is estimated robustly by making an
  initial Beta model, then trimming values below quantile trim/2 and
  above quantile 1-trim/2 of the initial model (0 \<= trim \< 1). If
  trimming reduces the number of correlations that remain to less than
  five, trimming is disabled and the full set of correlations is used.

## Value

A vector c(v, w) containing the two Beta parameters (see plot=TRUE
above).
