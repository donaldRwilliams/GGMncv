#' @importFrom numDeriv grad
#' @importFrom Rcpp sourceCpp
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
atan_deriv <- function(Theta, lambda, gamma = 0.1){
  Theta <- abs(Theta)
  lambda_mat <- lambda * ((gamma * (gamma + 2/pi)) / (gamma^2 + Theta^2))
  return(lambda_mat)
}

atan_pen <- function(Theta, lambda, gamma = 0.1){
  Theta <- abs(Theta)
  pen_mat <- lambda * (gamma + 2/pi) * atan(Theta/gamma)
  return(pen_mat)
  }

# Dicker, L., Huang, B., & Lin, X. (2013). Variable selection and estimation
# with the seamless-L0 penalty. Statistica Sinica, 929-962.
selo_pen <- function(Theta, lambda, gamma = 0.1) {
  Theta <- abs(Theta)
  pen_mat <- (lambda / log(2)) * log((Theta / (Theta + gamma)) + 1)
  return(pen_mat)
}

selo_deriv <- function(Theta, lambda, gamma = 0.1){
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
exp_deriv <- function(Theta, lambda, gamma = 0.1){
  Theta <- abs(Theta)
  lambda_mat <- (lambda/gamma) * exp(-(Theta/gamma))
  return(lambda_mat)
  }

# Mazumder, R., Friedman, J. H., & Hastie, T. (2011). Sparsenet: Coordinate descent
# with nonconvex penalties. Journal of the American Statistical Association, 106(495), 1125-1138.
log_pen <- function(x, lambda, gamma){
  gamma <- 1/gamma
  inv <- abs(x)
  pen_mat <- ((lambda / log(gamma+ 1)) * log(gamma * inv + 1))
  return(pen_mat)
}

log_deriv <- function(Theta, lambda, gamma = 0.1){
  p <- ncol(Theta)
  Theta <- abs(Theta)
  lambda_mat <- matrix(numDeriv::grad(log_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}


lasso_deriv <- function(Theta, lambda, gamma = 0){
  p <- ncol(Theta)
  lambda_mat <- matrix(lambda, p, p)
}

lq_deriv <- function(Theta, lambda, gamma = 0.5){
  Theta <- abs(Theta)
  epsilon <- 0.0001
  lambda_mat <- lambda * gamma * (Theta + epsilon)^(gamma - 1)
  lambda_mat
}

# Lv, J., & Fan, Y. (2009). A unified approach to model selection and sparse recovery using
# regularized least squares. The Annals of Statistics, 37(6A), 3498-3528.
sica_pen <- function(x, lambda, gamma){
  Theta <- x
  Theta <- abs(Theta)
  pen_mat <- lambda * (((gamma + 1) * Theta) /(Theta+gamma))
  return(pen_mat)
}

sica_deriv <- function(Theta, lambda, gamma = 0.1){
  p <- ncol(Theta)
  Theta <- abs(Theta)
  lambda_mat <- matrix(numDeriv::grad(sica_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}

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

globalVariables(c("VIP", "new1", "Y", "cs", "value", "X1", "X2", "ic", "lambda"))
