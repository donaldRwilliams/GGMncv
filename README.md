
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/hex.png" width = 250 />

# GGMncv: Gaussian Graphical Models with Non-Convex Penalties

[![Build
Status](https://travis-ci.org/donaldRwilliams/GGMncv.svg?branch=master)](https://travis-ci.org/donaldRwilliams/GGMncv)

The goal of GGMncv is to provide non-convex penalties for estimating
Gaussian graphical models. These are known to overcome the various
limitations of lasso, including (but not limited to) inconsistent model
selection, biased\[1\] estimates, and a high false positive rate.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

<!-- released version of GGMncv from -->

<!-- [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("GGMncv") -->

<!-- ``` -->

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMncv")
```

## Penalties

The following are implemented in `GGMncv`:

1.  Atan (`penalty = "atan"`; Wang and Zhu 2016). This is currently the
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

The methods in **GGMncv** are currently tuning free. This is
accomplished by setting the tuning parameter to `sqrt(log(p)/n)` (see
for example Zhang, Ren, and Chen 2018; Li et al. 2015; Jankova and Van
De Geer 2015).

## Example

A GGM can be fitted as follows

``` r
library(GGMncv)

# data
Y <- GGMncv::ptsd[,1:10]

# polychoric
S <- psych::polychoric(Y)$rho

# fit model
fit <- GGMncv(S, n = nrow(Y), 
              penalty = "atan", 
              LLA = TRUE)

# print
fit

#>       1     2     3     4     5     6     7     8     9    10
#> 1  0.000 0.255 0.000 0.309 0.101 0.000 0.000 0.000 0.073 0.000
#> 2  0.255 0.000 0.485 0.000 0.000 0.000 0.122 0.000 0.000 0.000
#> 3  0.000 0.485 0.000 0.185 0.232 0.000 0.000 0.000 0.000 0.000
#> 4  0.309 0.000 0.185 0.000 0.300 0.000 0.097 0.000 0.000 0.243
#> 5  0.101 0.000 0.232 0.300 0.000 0.211 0.166 0.000 0.000 0.000
#> 6  0.000 0.000 0.000 0.000 0.211 0.000 0.234 0.079 0.000 0.000
#> 7  0.000 0.122 0.000 0.097 0.166 0.234 0.000 0.000 0.000 0.000
#> 8  0.000 0.000 0.000 0.000 0.000 0.079 0.000 0.000 0.000 0.114
#> 9  0.073 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.261
#> 10 0.000 0.000 0.000 0.243 0.000 0.000 0.000 0.114 0.261 0.000
```

## Bootstrapping

It might be tempting to perform a bootstrap and then attempt to
construct confidence intervals for the edges. However, in general, these
“confidence” intervals do not have the correct properties to be
considered confidence intervals (see
[Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)). This
sentiment is echoed in Section 3.1, “Why standard bootstrapping and
subsampling do not work As,” of Bühlmann, Kalisch, and Meier (2014):

> The (limiting) distribution of such a sparse estimator is non-Gaussian
> with point mass at zero, and this is the reason why standard bootstrap
> or subsampling techniques do not provide valid confidence regions or
> p-values (pp. 7-8).

For this reason, it is common to not provide the standard errors for
penalized models. **GGMncv** follows the idea of behind the
**penalized** `R` package:

> It is a very natural question to ask for standard errors of regression
> coefficients or other estimated quantities. In principle such standard
> errors can easily be calculated, e.g. using the bootstrap. Still, this
> package deliberately does not provide them. The reason for this is
> that standard errors are not very meaningful for strongly biased
> estimates such as arise from penalized estimation methods (p.18,
> Goeman, Meijer, and Chaturvedi 2018)

Thus, at this time, confidence intervals are not provided for the
partial correlations. However, **GGMncv** does include the so-called
variable inclusion “probability” for each relation (see p. 1523 in Bunea
et al. 2011; and Figure 6.7 in Hastie, Tibshirani, and Wainwright 2015).
These are computed using a non-parametric bootstrap strategy.

### Variable Inclusion “Probability”

``` r
# data
Y <- GGMncv::ptsd[,1:10]

# polychoric
S <- psych::polychoric(Y)$rho

# fit model
fit <- GGMncv(S, n = nrow(Y), 
              penalty = "atan", 
              vip = TRUE)

# plot
plot(fit, size = 4)
```

![](man/figures/vip.png)

## References

<div id="refs" class="references">

<div id="ref-bunea2011penalized">

Bunea, Florentina, Yiyuan She, Hernando Ombao, Assawin Gongvatana, Kate
Devlin, and Ronald Cohen. 2011. “Penalized Least Squares Regression
Methods and Applications to Neuroimaging.” *NeuroImage* 55 (4): 1519–27.

</div>

<div id="ref-Buhlmann2014">

Bühlmann, Peter, Markus Kalisch, and Lukas Meier. 2014.
“High-Dimensional Statistics with a View Toward Applications in
Biology.” *Annual Review of Statistics and Its Application* 1 (1):
255–78. <https://doi.org/10.1146/annurev-statistics-022513-115545>.

</div>

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

<div id="ref-goeman2018l1">

Goeman, Jelle, Rosa Meijer, and Nimisha Chaturvedi. 2018. “L1 and L2
Penalized Regression Models.”

</div>

<div id="ref-hastie2015statistical">

Hastie, Trevor, Robert Tibshirani, and Martin Wainwright. 2015.
*Statistical Learning with Sparsity: The Lasso and Generalizations*. CRC
press.

</div>

<div id="ref-jankova2015confidence">

Jankova, Jana, and Sara Van De Geer. 2015. “Confidence Intervals for
High-Dimensional Inverse Covariance Estimation.” *Electronic Journal of
Statistics* 9 (1): 1205–29.

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

Wang, Yanxin, and Li Zhu. 2016. “Variable Selection and Parameter
Estimation with the Atan Regularization Method.” *Journal of Probability
and Statistics*.

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

1.  Note that the penalties in **GGMncv** should be *nearly* unbiased.
