
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/hex.png" width = 200 />

# GGMncv: Gaussian Graphical Models with Non-Convex Penalties

[![Build
Status](https://travis-ci.org/donaldRwilliams/GGMncv.svg?branch=master)](https://travis-ci.org/donaldRwilliams/GGMncv)

The goal of GGMncv is to provide non-convex penalties for estimating
Gaussian graphical models. These are known to overcome the various
limitations of lasso, including (but not limited to) inconsistent model
selection (Zhao and Yu 2006), biased <span id="a1">[\[1\]](#f1)</span>
estimates (Zhang 2010), and a high false positive rate (see for example
Williams and Rast 2020; Williams et al. 2019).

Note that these limitations of lasso are well-known. In the case of
false positives, for example, it has been noted that

> The lasso is doing variable screening and, hence, I suggest that we
> interpret the second ‘s’ in lasso as ‘screening’ rather than
> ‘selection’. Once we have the screening property, the task is to
> remove the false positive selections (p. 278, Tibshirani 2011).

Hence, contrary to the popular view of lasso (least absolute shrinkage
“screening” operator), using it for model *selection* is not
necessarily ideal. There are various ways to remove the false positives,
including thresholding after model selection (i.e., removing small
relations, Loh and Wainwright 2012) and two-stage procedures (Zou 2006).
The approach in **GGMncv**, on the other hand, selects the graph with
non-convex penalization (with *L*<sub>1</sub> as a special case).

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

4.  Smooth integration of counting and absolute deviation (`penalty =
    "sica"`; Lv, Fan, and others 2009)

5.  *L*<sub>q</sub> (`penalty = "lq"`; e.g., Knight and Fu 2000)

6.  Log (`penalty = "log"`; Mazumder, Friedman, and Hastie 2011)

7.  Smoothly clipped absolute deviation (`penalty = "scad"`; Fan and Li
    2001)

8.  Minimax concave penalty (`penalty = "mcp"`; Zhang 2010)

Options 1-4 are continuous approximations to the *L*<sub>0</sub>
penalty, that is, best subsets model selection. However, the solution is
computationally efficient and solved with the local linear approximation
described in Fan, Feng, and Wu (2009) or the one-step approach described
in Zou and Li (2008).

Note that computing the non-convex solution is a challenging task.
However, section 3.3 in Zou and Li (2008) indicates that the one-step
approach is a viable approximation for a variety of non-convex
penalties, assuming the initial estimates are “good enough”.

## Tuning Parameter

### Tuning Free

The default approach in **GGMncv** is tuning free. This is accomplished
by setting the tuning parameter to `sqrt(log(p)/n)` (see for example
Zhang, Ren, and Chen 2018; Li et al. 2015; Jankova and Van De Geer
2015).

### Selection

It is also possible to select the tuning parameter with BIC. This is
accomplished by setting `select = TRUE`.

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

**GGMncv** does not provide confidence intervals. This is because, in
general, “confidence” intervals from penalized approaches do not have
the correct properties to be considered confidence intervals [see
Wikipedia](https://en.wikipedia.org/wiki/Confidence_interval)). This
sentiment is echoed in Section 3.1, “Why standard bootstrapping and
subsampling do not work,” of Bühlmann, Kalisch, and Meier (2014):

> The (limiting) distribution of such a sparse estimator is non-Gaussian
> with point mass at zero, and this is the reason why standard bootstrap
> or subsampling techniques do not provide valid confidence regions or
> p-values (pp. 7-8).

For this reason, it is common to not provide standard errors (and thus
confidence intervals) for penalized models
<span id="a2">[\[2\]](#f2)</span>. **GGMncv** follows the idea of behind
the **penalized** `R` package:

> It is a very natural question to ask for standard errors of regression
> coefficients or other estimated quantities. In principle such standard
> errors can easily be calculated, e.g. using the bootstrap. Still, this
> package deliberately does not provide them. The reason for this is
> that standard errors are not very meaningful for strongly biased
> estimates such as arise from penalized estimation methods (p.18,
> Goeman, Meijer, and Chaturvedi 2018)

However, **GGMncv** does include the so-called variable inclusion
“probability” for each relation (see p. 1523 in Bunea et al. 2011; and
Figure 6.7 in Hastie, Tibshirani, and Wainwright 2015). These are
computed using a non-parametric bootstrap strategy.

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

## Citing **GGMncv**

It is important to note that **GGMncv** merely provides a software
implementation of other researchers work. There are no methological
innovations. Hence, in addition to citing the package
`citation("GGMncv")`, it is important to give credit to the primary
sources.

  - Atan

> (<span class="citeproc-not-found" data-reference-id="article">**???**</span>){wang2016variable,
> title={Variable selection and parameter estimation with the Atan
> regularization method}, author={Wang, Yanxin and Zhu, Li},
> journal={Journal of Probability and Statistics}, year={2016},
> publisher={Hindawi} }\`

## Footnotes

1.  <span id="f1"></span> Note that the penalties in **GGMncv** should
    provide *nearly* unbiased estimates. [(return)](#a1)

2.  <span id="f2"></span> It is possible to compute confidence intervals
    for lasso with the methods included in the **SILGGM** `R` package.
    These do not use the bootstrap (Zhang, Ren, and Chen 2018)  
    [(return)](#a2)

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
Penalized Regression Models.” *Vignette R Package Penalized.*

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

<div id="ref-knight2000asymptotics">

Knight, Keith, and Wenjiang Fu. 2000. “Asymptotics for Lasso-Type
Estimators.” *Annals of Statistics*, 1356–78.

</div>

<div id="ref-li2015flare">

Li, Xingguo, Tuo Zhao, Xiaoming Yuan, and Han Liu. 2015. “The Flare
Package for High Dimensional Linear Regression and Precision Matrix
Estimation in R.” *Journal of Machine Learning Research: JMLR* 16: 553.

</div>

<div id="ref-loh2012structure">

Loh, Po-Ling, and Martin J Wainwright. 2012. “Structure Estimation for
Discrete Graphical Models: Generalized Covariance Matrices and Their
Inverses.” In *Advances in Neural Information Processing Systems*,
2087–95.

</div>

<div id="ref-lv2009unified">

Lv, Jinchi, Yingying Fan, and others. 2009. “A Unified Approach to Model
Selection and Sparse Recovery Using Regularized Least Squares.” *The
Annals of Statistics* 37 (6A): 3498–3528.

</div>

<div id="ref-mazumder2011sparsenet">

Mazumder, Rahul, Jerome H Friedman, and Trevor Hastie. 2011. “Sparsenet:
Coordinate Descent with Nonconvex Penalties.” *Journal of the American
Statistical Association* 106 (495): 1125–38.

</div>

<div id="ref-tibshirani2011regression">

Tibshirani, Robert. 2011. “Regression Shrinkage and Selection via the
Lasso: A Retrospective.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 73 (3): 273–82.

</div>

<div id="ref-wang2018variable">

Wang, Yanxin, Qibin Fan, and Li Zhu. 2018. “Variable Selection and
Estimation Using a Continuous Approximation to the L0 Penalty.” *Annals
of the Institute of Statistical Mathematics* 70 (1): 191–214.

</div>

<div id="ref-wang2016variable">

Wang, Yanxin, and Li Zhu. 2016. “Variable Selection and Parameter
Estimation with the Atan Regularization Method.” *Journal of Probability
and Statistics*.

</div>

<div id="ref-williams2020back">

Williams, Donald R, and Philippe Rast. 2020. “Back to the Basics:
Rethinking Partial Correlation Network Methodology.” *British Journal of
Mathematical and Statistical Psychology* 73 (2): 187–212.

</div>

<div id="ref-williams2019nonregularized">

Williams, Donald R, Mijke Rhemtulla, Anna C Wysocki, and Philippe Rast.
2019. “On Nonregularized Estimation of Psychological Networks.”
*Multivariate Behavioral Research* 54 (5): 719–50.

</div>

<div id="ref-zhang2010nearly">

Zhang, Cun-Hui. 2010. “Nearly Unbiased Variable Selection Under Minimax
Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942.

</div>

<div id="ref-zhang2018silggm">

Zhang, Rong, Zhao Ren, and Wei Chen. 2018. “SILGGM: An Extensive R
Package for Efficient Statistical Inference in Large-Scale Gene
Networks.” *PLoS Computational Biology* 14 (8): e1006369.

</div>

<div id="ref-zhao2006model">

Zhao, Peng, and Bin Yu. 2006. “On Model Selection Consistency of Lasso.”
*Journal of Machine Learning Research* 7 (Nov): 2541–63.

</div>

<div id="ref-zou2006adaptive">

Zou, Hui. 2006. “The Adaptive Lasso and Its Oracle Properties.” *Journal
of the American Statistical Association* 101 (476): 1418–29.

</div>

<div id="ref-zou2008one">

Zou, Hui, and Runze Li. 2008. “One-Step Sparse Estimates in Nonconcave
Penalized Likelihood Models.” *Annals of Statistics* 36 (4): 1509.

</div>

</div>
