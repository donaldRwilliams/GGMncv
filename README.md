
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GGMncv

[![Build
Status](https://travis-ci.org/donaldRwilliams/GGMncv.svg?branch=master)](https://travis-ci.org/donaldRwilliams/GGMncv)

The goal of GGMncv is to provide non-convex penalties for estimating
Gaussian graphical models. These are known to overcome the various
limitations of lasso, including (but not limited to) consistent model
selection, nearly unbiased estimates, and a lower false positive rate.

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

2.  Seamless L0 (`penalty = "selo"`; Dicker, Huang, and Lin 2013)

3.  Exponential (`penalty = "exp"`)

4.  Smoothly clipped absolute deviation (`penalty = scad`)

5.  Minimax concave penalty (`penalty = mcp`)

## Example

<div id="refs" class="references">

<div id="ref-dicker2013variable">

Dicker, Lee, Baosheng Huang, and Xihong Lin. 2013. “Variable Selection
and Estimation with the Seamless-L 0 Penalty.” *Statistica Sinica*,
929–62.

</div>

</div>
