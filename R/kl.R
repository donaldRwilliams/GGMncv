#' @title Kullback-Leibler Divergence
#'
#' @description Compute KL divergence for a multivariate normal distribution.
#'
#' @param true Matrix. The true precision matrix
#'             (inverse of the covariance matrix)
#'
#' @param estimate Matrix. The estimated precision matrix
#'                 (inverse of the covariance matrix)
#'
#' @param stein Logical. Should Stein's loss be computed
#'              (defaults to \code{TRUE})? Note KL divergence is
#'              half of Stein's loss.
#'
#' @return Numeric corresponding to KL divergence.
#'
#' @note A lower value is better, with a score of zero indicating that
#' the estimated precision matrix is identical to the true precision matrix.
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' # nodes
#' p <- 20
#'
#' main <- gen_net(p = p, edge_prob = 0.15)
#'
#' y <- MASS::mvrnorm(250, rep(0, p), main$cors)
#'
#' fit_l1 <- ggmncv(R = cor(y),
#'               n = nrow(y),
#'               penalty = "lasso",
#'               progress = FALSE)
#'
#' # lasso
#' kl_mvn(fit_l1$Theta, solve(main$cors))
#'
#' fit_atan <- ggmncv(R = cor(y),
#'               n = nrow(y),
#'               penalty = "atan",
#'               progress = FALSE)
#'
#' kl_mvn(fit_atan$Theta, solve(main$cors))
#'
#' }
kl_mvn <- function(true, estimate, stein = FALSE) {

  Theta <- true

  hatTheta <- estimate

  p <- ncol(Theta)

  Sigma <- solve(Theta)

  stein_loss <- sum(diag(Sigma %*% hatTheta)) - log(det(Sigma %*% hatTheta)) - p

  if(isFALSE(stein)){

    kl <- stein_loss * 0.5

    return(kl)

  } else {

    return(stein_loss)

  }

}
