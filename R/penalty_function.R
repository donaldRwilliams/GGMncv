scad_penalty <- function(inv, lambda, gamma = 3.7){
  lam_mat <-ifelse(abs(inv) <= lambda, lambda,
                   ifelse(lambda < abs(inv) & abs(inv) < lambda * gamma,
                           (lambda * gamma - abs(inv)) / (gamma - 1),
                           0) )
  lam_mat
}

mcp_penalty <- function(inv, lambda, gamma = 2){
  abs_inv <- abs(inv)
  lam_mat <- ifelse(abs_inv <= lambda * gamma,
                    lambda - (abs_inv/gamma), 0)
  lam_mat <- ifelse(lam_mat < 0, 0, lam_mat)
  lam_mat
}


atan_penalty <- function(inv, lambda, gamma = 0.005){

  lam_mat <- lambda * ( (gamma * (gamma + 2/pi)) / (gamma^2 + inv^2))

  lam_mat

  }

# no good so far
truncated_penalty <- function(inv, lambda, gamma){

  abs_inv <- abs(inv)

  lam_mat <- lambda * (abs_inv < gamma)
  lam_mat

  }

mlog_penalty <- function(inv, lambda, gamma){
  abs_inv <- abs(inv)
  lam_mat <- gamma * (abs_inv < gamma) + (lambda * gamma / abs_inv) * (abs_inv >= gamma)
  lam_mat
}

