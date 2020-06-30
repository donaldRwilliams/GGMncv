#' GGMncv
#'
#' @param Sigma Correlation matrix
#'
#' @param n Numeric. Sample size.
#'
#' @param penalty Character string. Which penalty should be used (defaults to \code{selo}).
#'
#' @param nonreg Logical. Should the precision matrix be reestimated, given the adjacency matrix
#'               (defaults to TRUE)? This provides the \strong{nonregularized}, maximum likelhood estimate,
#'               with constraints.
#'
#' @param threshold Logical. Should a threshold be applied (\code{1/sqrt(n-p-5)} ) ? Note that this
#'                  may be useful for small sample sizes to improve specificity (i.e., lower false positive rate),
#'                  but not needed otherwise (see details).
#'
#' @param LLA Logical. Should the local linear approximation be used for maxmizing the penalized likelihood ?
#'            The default is \code{FALSE} which results in a one-step approach (see details).
#'
#' @param steps Numeric. How many steps if \code{LLA = TRUE} ? Set to \code{Inf} for iterating until
#'              the algorithm converges.
#'
#' @return An object of class \code{ggmncv}, including:
#' @export
GGMncv <- function(Sigma, n,
                   penalty = "selo",
                   nonreg = TRUE,
                   LLA = FALSE,
                   steps = NULL,
                   threshold = FALSE){
  p <- ncol(Sigma)
  lambda <- sqrt(log(p)/n)
  Theta <- solve(Sigma)

  if(penalty == "selo"){
    lambda_mat <- selo_deriv(Theta = Theta, lambda = lambda)
    diag(lambda_mat) <- lambda
    }

  fit <- glassoFast::glassoFast(S = Sigma, rho = lambda_mat)
  adj <- ifelse(fit$wi == 0, 0, 1)

  if(threshold){
    adj <- ifelse(abs(cov2cor(fit$wi)) <  1/sqrt(n-p-5), 0, 1)
    } else {
      adj <- ifelse(fit$wi == 0, 0, 1)
      }

  rest <- htf(Sigma, adj)
  Theta <- rest$Theta
  Sigma <- rest$Sigma
  P <- -(cov2cor(Theta) - diag(p))

  returned_object <- list(Theta = Theta,
                          Sigma = Sigma,
                          P = P,
                          fit = fit,
                          adj = adj)

  class(returned_object) <- "ggmncv"

  return(returned_object)
}
