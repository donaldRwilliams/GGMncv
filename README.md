
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

1.  Atan (`penalty = "atan"`; Wang and Zhu, n.d.). This is currently the
    default.

2.  Seamless L0 (`penalty = "selo"`; Dicker, Huang, and Lin 2013)

3.  Exponential (`penalty = "exp"`; Wang, Fan, and Zhu 2018)

4.  Smoothly clipped absolute deviation (`penalty = scad`; Fan and Li
    2001)

5.  Minimax concave penalty `penalty = mcp`; Zhang and others (2010)\]

Note that options 1-3 are continuous approximations to the L0 penalty,
that is, best subsets model selection. However, the solution is
computationally efficient and solved with the local linear approximation
described in Fan, Feng, and Wu (2009) or the one-step approach described
in Zou and Li (2008).

## Tuning Parameter

`GGMncv` is currently tuning free. This is accomplished by setting the
tuning parameter to `sqrt(log(p)/n)` (see for example Zhang, Ren, and
Chen 2018; Li et al. 2015; Jankova, Van De Geer, and others 2015).

## Example

<div id="refs" class="references">

<div id="ref-dicker2013variable">

Dicker, Lee, Baosheng Huang, and Xihong Lin. 2013. “Variable Selection
and Estimation with the Seamless-L 0 Penalty.” *Statistica Sinica*,
929–62.

</div>

<div id="ref-fan2009network">

Fan, Jianqing, Yang Feng, and Yichao Wu. 2009. “Network Exploration via
the Adaptive Lasso and Scad Penalties.” *The Annals of Applied
Statistics* 3 (2): 521.

</div>

<div id="ref-fan2001variable">

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave
Penalized Likelihood and Its Oracle Properties.” *Journal of the
American Statistical Association* 96 (456): 1348–60.

</div>

<div id="ref-jankova2015confidence">

Jankova, Jana, Sara Van De Geer, and others. 2015. “Confidence Intervals
for High-Dimensional Inverse Covariance Estimation.” *Electronic Journal
of Statistics* 9 (1): 1205–29.

</div>

<div id="ref-li2015flare">

Li, Xingguo, Tuo Zhao, Xiaoming Yuan, and Han Liu. 2015. “The Flare
Package for High Dimensional Linear Regression and Precision Matrix
Estimation in R.” *Journal of Machine Learning Research: JMLR* 16: 553.

</div>

<div id="ref-wang2018variable">

Wang, Yanxin, Qibin Fan, and Li Zhu. 2018. “Variable Selection and
Estimation Using a Continuous Approximation to the \[L\_0\] Penalty.”
*Annals of the Institute of Statistical Mathematics* 70 (1): 191–214.

</div>

<div id="ref-wang2016variable">

Wang, Yanxin, and Li Zhu. n.d. “Variable Selection and Parameter
Estimation with the Atan Regularization Method.” *Journal of Probability
and Statistics* 2016.

</div>

<div id="ref-zhang2010nearly">

Zhang, Cun-Hui, and others. 2010. “Nearly Unbiased Variable Selection
Under Minimax Concave Penalty.” *The Annals of Statistics* 38 (2):
894–942.

</div>

<div id="ref-zhang2018silggm">

Zhang, Rong, Zhao Ren, and Wei Chen. 2018. “SILGGM: An Extensive R
Package for Efficient Statistical Inference in Large-Scale Gene
Networks.” *PLoS Computational Biology* 14 (8): e1006369.

</div>

<div id="ref-zou2008one">

Zou, Hui, and Runze Li. 2008. “One-Step Sparse Estimates in Nonconcave
Penalized Likelihood Models.” *Annals of Statistics* 36 (4): 1509.

</div>

</div>
