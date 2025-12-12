# Center correlations

The null distribution model is r ~ 2 Beta(v,w) - 1 The centering
(target) distribution is r_centered ~ 2 Beta(nu, nu) - 1

## Usage

``` r
centerBeta(r, v, w, nu)
```

## Arguments

- r:

  The correlation to center

- v:

  The first parameter of the Beta distribution

- w:

  The second parameter of the Beta distribution

- nu:

  The parameter of the target Beta distribution

## Value

The centered correlation
