#' Ledoit and Wolf Shrinkage Estimator
#'
#' @description Compute the Ledoit and Wolf shrinkage estimator of
#' the covariance matrix \insertCite{ledoit2004well}{GGMncv},
#' which can be used for the \code{initial} inverse covariance matrix
#' in \code{\link{ggmncv}}.
#'
#' @param Y A data matrix (or data.frame) of dimensions \emph{n} by \emph{p}.
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#' @return Inverse correlation matrix.
#'
#' @export
#'
#' @examples
#'
#' # ptsd
#' Y <- ptsd[,1:5]
#'
#' # shrinkage
#' ledoit_wolf(Y)
#'
#' # non-reg
#' solve(cor(Y))
ledoit_wolf <- function(Y, ...){

  if (!requireNamespace("nlshrink", quietly = TRUE)) {
    stop("Please install the '", "nlshrink", "' package.")
  }

  Sigma <- nlshrink::linshrink_cov(as.matrix(Y))

  R <- cov2cor(Sigma)

  Rinv <- solve(R)

  return(Rinv)
}


