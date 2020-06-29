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
  lambda_mat <- matrix(numDeriv::grad(selo_pen, x = Theta, lambda = lambda, gamma = gamma), p, p)
  return(lambda_mat)
}

# Zhang, C. H. (2010). Nearly unbiased variable selection under minimax concave penalty.
# The Annals of statistics, 38(2), 894-942.
mcp_deriv <- function(Theta, lambda, gamma = 2){
  Theta <- abs(Theta)
  lambda_mat <- (Theta <= gamma * lambda) * (lambda - Theta/gamma)
  return(lambda_mat)
  }

mcp_pen <- function(Theta, lambda, gamma = 2){
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
