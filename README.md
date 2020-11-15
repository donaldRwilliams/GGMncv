
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/hex.png" width = 200 />

# GGMncv: Gaussian Graphical Models with Non-Convex Penalties

[![CRAN
Version](http://www.r-pkg.org/badges/version/GGMncv)](https://cran.r-project.org/package=GGMncv)
[![Downloads](https://cranlogs.r-pkg.org/badges/GGMncv)](https://cran.r-project.org/package=GGMncv)
[![Build
Status](https://travis-ci.org/donaldRwilliams/GGMncv.svg?branch=master)](https://travis-ci.org/donaldRwilliams/GGMncv)

The primary goal of GGMncv is to provide non-convex penalties for
estimating Gaussian graphical models. These are known to overcome the
various limitations of lasso (least absolute shrinkage “screening”
operator), including inconsistent model selection (Zhao and Yu 2006),
biased estimates (Zhang 2010)<span id="a1">[\[1\]](#f1)</span>, and a
high false positive rate (see for example Williams and Rast 2020;
Williams et al. 2019).

Note that these limitations of lasso are well-known. In the case of
false positives, for example, it has been noted that

> The lasso is doing variable screening and, hence, I suggest that we
> interpret the second ‘s’ in lasso as ‘screening’ rather than
> ‘selection’. Once we have the screening property, the task is to
> remove the false positive selections (p. 278, Tibshirani 2011).

An additional goal of **GGMncv** is to provide methods for making
statistical inference in **regularized** Gaussian graphical models. This
is accomplished with the de-sparsified graphical lasso estimator
introduced in Jankova and Van De Geer (2015). This is described in the
section [De-Sparsified Estimator](#de-sparsified-estimator). The next
[section](#comparing-ggms) shows how the de-sparsified estimator can be
used to compare GGMs.

## Installation

You can install the released version (`1.1.0`) of GGMncv from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GGMncv")
```

You can install development version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("donaldRwilliams/GGMncv")
```

## Penalties

The following are implemented in `GGMncv`:

1.  Atan (`penalty = "atan"`; Wang and Zhu 2016). This is currently the
    default.

2.  Seamless *L*<sub>0</sub> (`penalty = "selo"`; Dicker, Huang, and Lin
    2013)

3.  Exponential (`penalty = "exp"`; Wang, Fan, and Zhu 2018)

4.  Smooth integration of counting and absolute deviation (`penalty =
    "sica"`; Lv and Fan 2009)

5.  Log (`penalty = "log"`; Mazumder, Friedman, and Hastie 2011)

6.  *L*<sub>q</sub> (`penalty = "lq"`, *0* \< *q* \< *1*; e.g., Knight
    and Fu 2000)

7.  Smoothly clipped absolute deviation (`penalty = "scad"`; Fan and Li
    2001)

8.  Minimax concave penalty (`penalty = "mcp"`; Zhang 2010)

9.  Adaptive lasso (`penalty = "adapt"`; Zou 2006)

10. Lasso (`penalty = "lasso"`; Tibshirani 1996)

Options 1-5 are continuous approximations to the *L*<sub>0</sub>
penalty, that is, best subsets model selection. However, the solution is
computationally efficient and solved with the local linear approximation
described in Fan, Feng, and Wu (2009) or the one-step approach described
in Zou and Li (2008).

### Penalty Function

The basic idea of these penalties is to provide “tapering,” in which
regularization is less severe for large effects. The following is an
example for the Atan penalty
(![](https://latex.codecogs.com/svg.latex?%5Cgamma) is the
hyperparameter)

![](man/figures/pen_func.png)

Note that (1) the penalty provides a “smooth” function that ranges from
*L*<sub>0</sub> (best subsets) and *L*<sub>1</sub> (lasso)
regularization; and (2) the penalty “tapers” off for large effects.

### Computation

Computing the non-convex solution is a challenging task. However,
section 3.3 in Zou and Li (2008) indicates that the one-step approach is
a viable **approximation** for a variety of non-convex penalties,
assuming the initial estimates are “good
enough”<span id="a2">[\[2\]](#f2)</span>. To this end, the initial
values can either be the sample based inverse covariance matrix or a
custom matrix specified with `initial`.

## Tuning Parameter

### Selection

The tuning parameter can be selected with several information criteria
(IC), including `aic`, `bic` (currently the default),`ebic`, `ric`, in
addition to any of the *generalized* information criteria provided in
section 5 of Kim, Kwon, and Choi (2012).

Information criterion can be understood as penalizing the likelihood,
with the difference being in the severity of the penalty. -2 times the
log-likelihood is defined as

![](https://latex.codecogs.com/svg.latex?-2l_n%28%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D%29%20%3D%20-2%20%5CBig%5B%5Cfrac%7Bn%7D%7B2%7D%20%5Ctext%7Blogdet%7D%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D%20-%20%5Ctext%7Btr%7D%28%5Chat%7B%5Ctext%7B%5Cbf%7BS%7D%7D%7D%20%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D%29%20%5CBig%5D)

where
![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D)
is the estimated precision matrix and
![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ctext%7B%5Cbf%7BS%7D%7D%7D)
is the sample based covariance matrix. The included criterion then add
the following penalties:

  - GIC<sub>1</sub> (BIC):
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%20%5Ctext%7Blog%7D%28n%29)

Note that
![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C)
refers to the cardinality of the edge set, that is, the number of edges.

  - GIC<sub>2</sub>:
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%20p%5E%7B1/3%7D)

*p* denotes the number of nodes or columns in the data matrix.

  - GIC<sub>3</sub> (RIC):
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%202%20%5Ctext%7Blog%7D%28p%29)

  - GIC<sub>4</sub>:
    ![](https://latex.codecogs.com/svg.latex?2%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%5Ctext%7Blog%7D%28p%29%20+%20%5Ctext%7Bloglog%7D%28p%29)

  - GIC<sub>5</sub> (BIC with divergent dimensions):
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%20%5Ctext%7Bloglog%7D%28n%29%5Ctext%7Blog%7D%28p%29)

  - GIC<sub>6</sub>:
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%20%5Ctext%7Blog%7D%28n%29%5Ctext%7Blog%7D%28p%29)

  - AIC:
    ![](https://latex.codecogs.com/svg.latex?2%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C)

Although cross-validation is not implemented for selecting the tuning
parameter, AIC can be used to approximate leave-one-out
cross-validation.

  - EBIC:
    ![](https://latex.codecogs.com/svg.latex?%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%5Ctext%7Blog%7D%28n%29%20+%204%20%7C%5Ctext%7B%5Cbf%7BE%7D%7D%7C%20%5Cgamma%20%5Ctext%7Blog%7D%28p%29%2C%5C%3B%20%5C0%20%5Cleq%20%5Cgamma%20%5Cleq%201)

The tuning parameter is selected by setting `select = TRUE` and then the
desired IC with, for example, `ic = "gic_3"`.

### Tuning Free

A tuning free option is also available. This is accomplished by setting
the tuning parameter to
![](https://latex.codecogs.com/svg.latex?%5Csmall%20%5Csqrt%7Blog%28p%29/n%7D)
(see for example Zhang, Ren, and Chen 2018; Li et al. 2015; Jankova and
Van De Geer 2015) and then selecting
![](https://latex.codecogs.com/svg.latex?%5Cgamma)

## Example: Structure Learning

A GGM can be fitted as follows

``` r
library(GGMncv)

# data
Y <- GGMncv::ptsd[,1:10]

# polychoric
R <- cor(Y, method = "spearman")

# fit model
fit <- ggmncv(R = R, n = nrow(Y), 
              penalty = "atan")

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

Note that the object `fit` can be plotted with the `R` package
[**qgraph**](https://CRAN.R-project.org/package=qgraph).

## Example: Out-of-Sample Prediction

The **GGMncv** package can also be used for prediction, given the
correspondence between the inverse covariance matrix and multiple
regression (Kwan 2014).

``` r
Y <- scale(Sachs)

# test data
Ytest <- Y[1:100,]

# training data
Ytrain <- Y[101:nrow(Y),]

# default: atan and tuning free
fit <- ggmncv(cor(Ytrain), 
              n = nrow(Ytrain))

# predict
pred <- predict(fit, train_data = Ytrain, newdata = Ytest)

# print mse
round(apply((pred - Ytest)^2, 2, mean), 2)

#>  Raf  Erk Plcg  PKC  PKA PIP2 PIP3  Mek  P38  Jnk  Akt 
#> 0.18 0.27 0.59 0.42 0.39 0.47 0.69 0.16 0.15 0.69 0.26 
```

## Solution Path

When `select = 'lambda'`, the solution path for the partial correlations
can be plotted.

### Atan Penalty

Here is the current default penalty

``` r
# data
Y <- ptsd

# fit model
fit <- GGMncv(cor(Y), n = nrow(Y), 
              select = TRUE, 
              store = TRUE)
              
# plot path
plot(fit, 
     alpha = 0.75)
```

![](man/figures/atan_path.png)

The dotted line is denotes the selected lambda. Notice how the larger
partial correlations “escape” regularization, at least to some degree,
compared to the smaller partial correlations.

### Lasso Penalty

Next *L*<sub>1</sub> regularization is implemented by setting `penalty =
"lasso"`.

``` r
# data
Y <- ptsd

# fit model
fit <- GGMncv(cor(Y), n = nrow(Y), 
              penalty = "lasso"
              store = TRUE)

# plot path
plot(fit, 
     alpha = 0.75)
```

![](man/figures/lasso_path.png)

This solution is much different than above. For example, it is clear
that the large partial correlations are heavily penalized, whereas this
was not so for the atan penalty. The reason this is not ideal is that,
if the partial correlations are large, it makes sense that they should
not be penalized that much. This property of non-convex regularization
should provide *nearly* unbiased estimates, which can improve, say,
predictive accuracy.

Also notice that the atan penalty provides a sparser solution.

## Bootstrapping

**GGMncv** does not provide confidence intervals based on bootstrapping.
This is because, in general, these kinds “confidence” intervals from
penalized approaches do not have the correct properties to be considered
confidence intervals [(see
Wikipedia)](https://en.wikipedia.org/wiki/Confidence_interval). This
sentiment is echoed in Section 3.1, “Why standard bootstrapping and
subsampling do not work,” of Bühlmann, Kalisch, and Meier (2014):

> The (limiting) distribution of such a sparse estimator is non-Gaussian
> with point mass at zero, and this is the reason why standard bootstrap
> or subsampling techniques do not provide valid confidence regions or
> p-values (pp. 7-8).

For this reason, it is common to **not** provide standard errors (and
thus confidence intervals) for penalized models
<span id="a3">[\[3\]](#f3)</span>. For example, this is from the
**penalized** `R` package:

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

Additionally, more recent work does allow for obtaining confidence
intervals and *p*-values with the de-sparsified method. For the
graphical lasso, the former are not available for the partial
correlations so currently only *p*-values are provided [(Statistical
Inferece)](#statistical-inference).

### Variable Inclusion “Probability”

``` r
# data
Y <- GGMncv::ptsd[,1:5]

# edge inclusion
eips <- boot_eip(Y)

# plot
plot(eips, size = 4)
```

![](man/figures/vip.png)

## Statistical Inference

It might be tempting to think these approaches lead to rich *inference*.
This would be a mistake–they suffer from all of the problems inherent to
automated procedures for model selection (e.g., Berk et al. 2013; Lee et
al. 2016).

And note that:

1.  Simply **not** detecting an effect does not provide evidence for the
    null hypothesis.

2.  There is **not** necessarily a difference between an effect that was
    and an effect that was not detected.

Supporting these claims would require a valid confidence interval that
has been corrected for model selection and/or regularization. With these
caveats in mind, data driven model selection in **GGMncv** can be used
for explicit data mining or prediction.

### De-Sparsified Estimator

To make inference, **GGMncv** computes the de-sparsified estimator,
![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D),
introduced in Jankova and Van De Geer (2015), that is

![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D%20%3D%202%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D%20-%20%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D%20%5Chat%7B%5Ctext%7B%5Cbf%7BR%7D%7D%7D%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D)

where
![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Cboldsymbol%7B%5CTheta%7D%7D)
is the estimated precision matrix and
![](https://latex.codecogs.com/svg.latex?%5Chat%7B%5Ctext%7B%5Cbf%7BR%7D%7D%7D)
is the sample based correlation matrix. As the name implies, this
removes the zeros and corrects the bias from regularization. The
asymptotic variance is then given as

![](https://latex.codecogs.com/svg.latex?%5Ctext%7BVar%7D%5B%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D%5D%20%3D%20%7B%5Ctext%7Bdiag%7D%28%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D%29%20%5Ctext%7Bdiag%7D%28%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D%29%5E%5Cprime%20+%20%5Chat%7B%5Ctext%7B%5Cbf%7BT%7D%7D%7D%5E2%7D)

which readily allows for computing *p*-values for each off-diagonal
element of the de-sparsified estimator.

This is implemented with

``` r
# data
Y <- ptsd

# fit model
fit <- ggmncv(cor(Y), n = nrow(Y))

# make inference
fdr_ggm <- inference(fit, method = "fdr")

# print
fdr_ggm

#> Statistical Inference
#> fdr: 0.05
#> ---

#>   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
#> 1  0 1 0 1 0 0 0 0 0  0  1  0  0  0  0  0  0  0  0  0
#> 2  1 0 1 0 0 0 0 0 0  0  0  0  0  0  0  1  0  0  0  0
#> 3  0 1 0 1 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0
#> 4  1 0 1 0 1 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0
#> 5  0 0 0 1 0 1 0 0 0  0  0  0  0  0  1  0  0  1  0  0
#> 6  0 0 0 0 1 0 1 0 0  0  0  0  0  0  0  0  0  0  0  0
#> 7  0 0 0 0 0 1 0 0 0  0  0  1  0  0  0  0  0  0  0  0
#> 8  0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  0  0
#> 9  0 0 0 0 0 0 0 0 0  0  1  1  0  0  0  0  0  0  0  0
#> 10 0 0 0 0 0 0 0 0 0  0  1  0  0  0  0  0  0  0  0  0
#> 11 1 0 0 0 0 0 0 0 1  1  0  0  0  0  1  0  0  0  0  0
#> 12 0 0 0 0 0 0 1 0 1  0  0  0  1  0  0  0  0  0  1  0
#> 13 0 0 0 0 0 0 0 0 0  0  0  1  0  1  0  0  0  0  1  0
#> 14 0 0 0 0 0 0 0 0 0  0  0  0  1  0  0  0  0  0  0  0
#> 15 0 0 0 0 1 0 0 0 0  0  1  0  0  0  0  1  0  0  0  0
#> 16 0 1 0 0 0 0 0 0 0  0  0  0  0  0  1  0  0  0  0  0
#> 17 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  1  0  0
#> 18 0 0 0 0 1 0 0 0 0  0  0  0  0  0  0  0  1  0  0  0
#> 19 0 0 0 0 0 0 0 0 0  0  0  1  1  0  0  0  0  0  0  1
#> 20 0 0 0 0 0 0 0 0 0  0  0  0  0  0  0  0  0  0  1  0
```

Note that the object `fdr_ggm` includes the de-sparsified precision
matrix, the partial correlation matrix, and *p*-values for each
relation. Furthermore, there is a function called `desparsify()` that
can be used to obtain the de-sparsified estimator without computing the
*p*-values.

## Comparing GGMs

Because the de-sparsified estimator provides the variance for each
relation, this readily allows for comparing GGMs. This is accomplished
by computing the difference and then the variance of that difference.
Assuming there is two groups, `Y_g1` and `Y_g2`, this is implemented
with

``` r
fit1 <- ggmncv(Y_g1, n = nrow(Y_g1))
fit2 <- ggmncv(Y_g2, n = nrow(Y_g2))

ggm_diff <- ggm_compare(fit1, fit2)
```

The object `ggm_diff` includes the partial correlation differences,
*p*-values, and the adjacency matrix.

## Citing **GGMncv**

It is important to note that **GGMncv** merely provides a software
implementation of other researchers work. There are no methological
innovations, although this is the most comprehensive `R` package for
estimating GGMs with non-convex penalties. Hence, in addition to citing
the package `citation("GGMncv")`, it is important to give credit to the
primary sources. The references can be found in
[(Penalties)](#penalties).

Additionally, please cite Williams (2020) which is survey of these
approaches that is meant to accompany **GGMncv**.

## Footnotes

1.  <span id="f1"></span> Note that the penalties in **GGMncv** should
    provide *nearly* unbiased estimates [(return)](#a1).

2.  <span id="f2"></span> In low-dimensional settings, assuming that *n*
    is sufficiently larger than *p*, the sample covariance matrix
    provides adequate initial estimates. In high-dimensional settings
    (*n* \< *p*), the initial estimates are obtained from lasso
    [(return)](#a2).

3.  <span id="f3"></span> It is possible to compute confidence intervals
    for lasso with the methods included in the **SILGGM** `R` package
    (Zhang, Ren, and Chen 2018). These do not use the bootstrap
    [(return)](#a3).

## References

<div id="refs" class="references">

<div id="ref-berk2013valid">

Berk, Richard, Lawrence Brown, Andreas Buja, Kai Zhang, Linda Zhao, and
others. 2013. “Valid Post-Selection Inference.” *The Annals of
Statistics* 41 (2): 802–37.

</div>

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

<div id="ref-kim2012consistent">

Kim, Yongdai, Sunghoon Kwon, and Hosik Choi. 2012. “Consistent Model
Selection Criteria on High Dimensions.” *The Journal of Machine Learning
Research* 13: 1037–57.

</div>

<div id="ref-knight2000asymptotics">

Knight, Keith, and Wenjiang Fu. 2000. “Asymptotics for Lasso-Type
Estimators.” *Annals of Statistics*, 1356–78.

</div>

<div id="ref-kwan2014regression">

Kwan, Clarence CY. 2014. “A Regression-Based Interpretation of the
Inverse of the Sample Covariance Matrix.” *Spreadsheets in Education* 7
(1): 4613.

</div>

<div id="ref-lee2016exact">

Lee, Jason D, Dennis L Sun, Yuekai Sun, Jonathan E Taylor, and others.
2016. “Exact Post-Selection Inference, with Application to the Lasso.”
*The Annals of Statistics* 44 (3): 907–27.

</div>

<div id="ref-li2015flare">

Li, Xingguo, Tuo Zhao, Xiaoming Yuan, and Han Liu. 2015. “The Flare
Package for High Dimensional Linear Regression and Precision Matrix
Estimation in R.” *Journal of Machine Learning Research: JMLR* 16: 553.

</div>

<div id="ref-lv2009unified">

Lv, Jinchi, and Yingying Fan. 2009. “A Unified Approach to Model
Selection and Sparse Recovery Using Regularized Least Squares.” *The
Annals of Statistics* 37 (6A): 3498–3528.

</div>

<div id="ref-mazumder2011sparsenet">

Mazumder, Rahul, Jerome H Friedman, and Trevor Hastie. 2011. “Sparsenet:
Coordinate Descent with Nonconvex Penalties.” *Journal of the American
Statistical Association* 106 (495): 1125–38.

</div>

<div id="ref-tibshirani1996regression">

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the
Lasso.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 58 (1): 267–88.

</div>

<div id="ref-tibshirani2011regression">

———. 2011. “Regression Shrinkage and Selection via the Lasso: A
Retrospective.” *Journal of the Royal Statistical Society: Series B
(Statistical Methodology)* 73 (3): 273–82.

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

<div id="ref-williams2020beyond">

Williams, Donald R. 2020. “Beyond Lasso: A Survey of Nonconvex
Regularization in Gaussian Graphical Models.” *PsyArXiv*.

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
