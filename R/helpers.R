#' @importFrom numDeriv grad
#' @importFrom Rcpp sourceCpp
#' @importFrom MASS mvrnorm
#' @importFrom stats deriv
#' @useDynLib GGMncv, .registration=TRUE

# Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized likelihood and
# its oracle properties. Journal of the American statistical Association, 96(456), 1348-1360.
scad_deriv <- function(Theta, lambda, gamma = 3.7){
  Theta <- abs(Theta)
  lambda_mat <- lambda * (Theta <= lambda) + (pmax(gamma * lambda - Theta, 0) / (gamma - 1)) * (Theta >lambda)
  return(lambda_mat)
}

scad_pen <- function(Theta, lambda, gamma = 3.7){
  Theta <- abs(Theta)
  pen_mat <- (Theta <= lambda) * lambda * Theta  +
    (lambda < Theta & Theta <= gamma * lambda) * -(Theta^2 - (2*gamma*lambda*Theta) + lambda^2) / ((gamma-1)*2) +
    (Theta > gamma * lambda) * ((gamma + 1) * lambda ^2)/2
  return(pen_mat)
}

# Wang, Y., & Zhu, L. (2016). Variable selection and parameter estimation
# with the Atan regularization method. Journal of Probability and Statistics, 2016.
atan_deriv <- function(Theta, lambda, gamma = 0.01){
  Theta <- abs(Theta)
  lambda_mat <- lambda * ((gamma * (gamma + 2/pi)) / (gamma^2 + Theta^2))
  return(lambda_mat)
}

atan_pen <- function(Theta, lambda, gamma = 0.01){
  Theta <- abs(Theta)
  pen_mat <- lambda * (gamma + 2/pi) * atan(Theta/gamma)
  return(pen_mat)
  }

# Dicker, L., Huang, B., & Lin, X. (2013). Variable selection and estimation
# with the seamless-L0 penalty. Statistica Sinica, 929-962.
selo_pen <- function(Theta, lambda, gamma = 0.01) {
  Theta <- abs(Theta)
  pen_mat <- (lambda / log(2)) * log((Theta / (Theta + gamma)) + 1)
  return(pen_mat)
}

selo_deriv <- function(Theta, lambda, gamma = 0.01){
  p <- ncol(Theta)
  Theta <- abs(Theta)
  Theta <- ifelse(Theta == 0, 1e-5, Theta)
  lambda_mat <- matrix(numDeriv::grad(selo_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}

# Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty.
# The Annals of statistics, 38(2), 894-942.
mcp_deriv <- function(Theta, lambda, gamma = 3){
  Theta <- abs(Theta)
  lambda_mat <- (Theta <= gamma * lambda) * (lambda - Theta/gamma)
  return(lambda_mat)
  }

mcp_pen <- function(Theta, lambda, gamma = 3){
  Theta <- abs(Theta)
  pen_mat <- (Theta <= lambda * gamma) *  (lambda * Theta - (Theta^2/(2*gamma))) +
    (Theta > lambda * gamma) * (0.5 * gamma * lambda^2)
  return(pen_mat)
}

# Wang, Y., Fan, Q., & Zhu, L. (2018). Variable selection and estimation using a
# continuous approximation to the $$ L_0 $$ penalty. Annals of the
# Institute of Statistical Mathematics, 70(1), 191-214.
exp_deriv <- function(Theta, lambda, gamma = 0.01){
  Theta <- abs(Theta)
  lambda_mat <- (lambda/gamma) * exp(-(Theta/gamma))
  return(lambda_mat)
  }

exp_pen <- function(Theta, lambda, gamma = 0.01){
  Theta <- abs(Theta)
  pen_mat <- lambda * (1 - exp(-(Theta/gamma)))
  return(pen_mat)
}

# Mazumder, R., Friedman, J. H., & Hastie, T. (2011). Sparsenet: Coordinate descent
# with nonconvex penalties. Journal of the American Statistical Association, 106(495), 1125-1138.
log_pen <- function(Theta, lambda, gamma = 0.01){
  gamma <- 1/gamma
  Theta <- abs(Theta + 0.0001)
  pen_mat <- ((lambda / log(gamma+ 1)) * log(gamma * Theta + 1))
  return(pen_mat)
}

log_deriv <- function(Theta, lambda, gamma = 0.01){
  p <- ncol(Theta)
  Theta <- abs(Theta)
  lambda_mat <- matrix(numDeriv::grad(log_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}

lasso_deriv <- function(Theta, lambda, gamma = 0){
  p <- ncol(Theta)
  lambda_mat <- matrix(lambda, p, p)
}

lq_pen <- function(Theta, lambda, gamma = 0.5){
  Theta <- abs(Theta)
  epsilon <- 0.0001
  pen_mat <- lambda * ((Theta + epsilon)^gamma)
  return(pen_mat)
}

lq_deriv <- function(Theta, lambda, gamma = 0.5){
  Theta <- abs(Theta)
  epsilon <- 0.0001
  lambda_mat <- lambda * gamma * (Theta + epsilon)^(gamma - 1)
  return(lambda_mat)
}

# Lv, J., & Fan, Y. (2009). A unified approach to model selection and sparse recovery using
# regularized least squares. The Annals of Statistics, 37(6A), 3498-3528.
sica_pen <- function(Theta, lambda, gamma){
  Theta <- abs(Theta + 0.0001)
  pen_mat <- lambda * (((gamma + 1) * Theta) /(Theta+gamma))
  return(pen_mat)
}

sica_deriv <- function(Theta, lambda, gamma = 0.01){
  p <- ncol(Theta)
  Theta <- abs(Theta)
  lambda_mat <- matrix(numDeriv::grad(sica_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}

# Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American
# statistical association, 101(476), 1418-1429.
adapt_deriv <- function(Theta, lambda, gamma = 0.5){
  Theta <- abs(Theta + 0.0001)
  # note: 1-gamma for consistency (-> 0 large weight)
  lambda_mat <- lambda * Theta^(-(1-gamma))
  return(lambda_mat)
}


# Kim, Y., Kwon, S., & Choi, H. (2012). Consistent model selection criteria on high
# dimensions. The Journal of Machine Learning Research, 13, 1037-1057.
gic_helper <- function(Theta, R, edges, n, p, type = "bic", ...) {
  log.like <- (n / 2) * (log(det(Theta)) - sum(diag(R %*% Theta)))

  neg_ll <- -2 * log.like

  if (type == "bic" | type == "gic_1") {
    ic <- neg_ll + edges * log(n)
  } else if (type == "aic") {
    ic <- neg_ll + edges * 2
  } else if (type == "gic_2") {
    ic <- neg_ll + edges * p ^ (1 / 3)
  } else if (type == "ric" | type == "gic_3") {
    ic <- neg_ll + edges * 2 * log(p)
  } else if (type == "gic_4") {
    ic <- neg_ll + edges * 2 * (log(p) + log(log(p)))
  } else if (type == "gic_5") {
    ic <- neg_ll + edges * log(log(n)) * log(p)
  } else if (type == "gic_6") {
    ic <- neg_ll + edges * log(n) * log(p)
  } else if (type == "ebic") {
    dots <- list(...)
    if (is.null(dots$ebic_gamma)) {
      gamma <- 0.5
    } else {
      gamma <- dots$ebic_gamma
    }
    ic <- neg_ll + edges * log(n) + 4 * edges * gamma * log(p)
  } else {
    stop("ic not found. see documentation")
  }

  return(ic)

}


# note this is implemented in c++ for speed.
# This is for testing purposes
htf <- function(Sigma, adj, tol = 1e-10) {
  S <- Sigma
  p <- ncol(S)
  diag(adj) <- 0
  W <- W_previous <- S
  n_iter <- 0
  repeat {
    for (i in 1:p) {
      beta <- rep(0, p - 1)
      pad_index <- which(adj[i, -i] == 1)
      if (length(pad_index) == 0) {
        w_12 <- beta
      }
      else {
        W_11 <- W[-i, -i]
        s_12 <- S[i, -i]
        W_11_star <- W_11[pad_index, pad_index]
        s_12_star <- s_12[pad_index]
        beta[pad_index] <- solve(W_11_star) %*% s_12_star
        w_12 <- W_11 %*% beta
      }
      W[-i, i] <- W[i, -i] <- w_12
    }
    max_diff <- max(W_previous[upper.tri(W)] - W[upper.tri(W)])
    if (max_diff < tol) {
      break
    }
    else {
      W_previous <- W
      n_iter <- n_iter + 1
    }
  }
  returned_object <- list(Theta = round(solve(W), 4), Sigma = round(W, 4))
  returned_object
}

coef_helper <- function(Theta){
  p <- ncol(Theta)
  betas <- round(t(sapply(1:p, function(x) Theta[x,-x] / Theta[x,x])), 3)  * -1
  return(betas)
}

check_gamma <- function(penalty, gamma){

  if (penalty == "scad") {
    if (any(gamma <= 3)) {
      stop("gamma must be greater than 3")
    }
  } else if (penalty == "mcp") {
    if (any(gamma <= 1)) {
      stop("gamma must be greater than 1")
    }

  } else {
    if (any(gamma < 0)) {
      stop("gamma must be positive")
    }
  }
}

# taken from
# Kuismin, M., & Sillanpää, M. J. (2016). Use of Wishart prior and simple extensions for
# sparse precision matrix estimation. PloS one, 11(2), e0148171.
lw_helper <- function(x, n){

  p <- ncol(x)

  if (isSymmetric(as.matrix(x))) {

    Y <- MASS::mvrnorm(
      n = n,
      mu = rep(0, p),
      Sigma = x,
      empirical = TRUE
    )

  } else {

    Y <- x
  }
  # Y (n x p)    : n iid observations on p random variables
  # Sigma (p x p): invertible covariance matrix estimator
  #
  # Shrinks towards one-parameter matrix:
  #    all variances are the same
  #    all covariances are zero

  # Modified from the MATLAB-code downloaded from the website of Michael Wolf in the Department of Economics of the University of Zurich.
  # Based on the presentation in the article of Ledoit & Wolf (2004): "Honey, I Shrunk The Sample Covariance Matrix". The Journal of Portfolio Management Summer, Vol. 30, No. 4, 110-119.

  # This version: 2/2015

  ############################################################################

  # This file is released under the BSD 2-clause license.

  # Copyright (c) 2014, Olivier Ledoit and Michael Wolf
  # All rights reserved.

  # Redistribution and use in source and binary forms, with or without
  # modification, are permitted provided that the following conditions are
  # met:

  # 1. Redistributions of source code must retain the above copyright notice,
  # this list of conditions and the following disclaimer.

  # 2. Redistributions in binary form must reproduce the above copyright
  # notice, this list of conditions and the following disclaimer in the
  # documentation and/or other materials provided with the distribution.

  # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
  # IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
  # THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  # PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
  # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  # LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  # NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  ##############################################################################
  # de-mean returns

  Y <- scale(Y, scale = FALSE)

  # compute S covariance matrix
  S <- crossprod(Y)/n

  # compute prior

  meanvar <- mean(diag(S))
  prior <- meanvar*diag(1,p)

  # what we call b
  X <- Y^2
  phiMat <- (crossprod(X)/n) - S^2
  phi <- sum(phiMat)

  # what we call c
  gamma = sum(abs(S-prior)^2)

  # compute shrinkage constant
  kappa <- phi/gamma
  shrinkage <- max(0,min(1,kappa/n))

  Sigma <- shrinkage*prior+(1-shrinkage)*S
  Theta <- solve(cov2cor(Sigma))
  return(Theta)

}

globalVariables(c("VIP",
                  "new1",
                  "Y",
                  "cs",
                  "value",
                  "X1",
                  "X2",
                  "ic",
                  "lambda",
                  "coef",
                  "EIP",
                  "dat",
                  "thetas",
                  "pen"))
