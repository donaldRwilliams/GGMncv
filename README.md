
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GGMncv

[![Build
Status](https://travis-ci.org/donaldRwilliams/GGMncv.svg?branch=master)](https://travis-ci.org/donaldRwilliams/GGMncv)

The goal of GGMncv is to provide variance non-convex penalties for
estimating Gaussian graphical models. These are known to overcome the
various limitation of \(L\)

## Installation

You can install the released version of GGMncv from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GGMncv")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMncv")
```

## Penalties

The following are implemented in `GGMncv`

1.  Arctan (X). This is currently the default (`penalty = "atan"`).

2.  Seamless L0 (`penalty = "selo"`)

3.  Exponential (`penalty = "exp"`)

4.  Smoothly clipped absolute deviation (`penalty = scad`)

5.  Minimax concave penality (`penalty = mcp`)

## Example
