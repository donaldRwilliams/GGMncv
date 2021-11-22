#' De-Sparsified Graphical Lasso Estimator
#'
#' @description
#' \loadmathjax
#' Compute the de-sparsified (sometimes called "de-biased") glasso estimator with
#' the approach described in Equation 7 of \insertCite{jankova2015confidence;textual}{GGMncv}.
#' The basic idea is to \emph{undo} \mjseqn{L_1}-regularization, in order
#' to compute \emph{p}-values and confidence intervals
#' (i.e., to make statistical inference).
#'
#' @param object An object of class \code{ggmncv}.
#'
#' @param ... Currently ignored.
#'
#' @return The de-sparsified estimates, including
#'
#' \itemize{
#'
#' \item \code{Theta}:  De-sparsified precision matrix
#'
#' \item \code{P}:  De-sparsified partial correlation matrix
#'
#' }
#'
#'
#' @details
#' According to \insertCite{jankova2015confidence;textual}{GGMncv}, the de-sparisifed estimator,
#' \mjseqn{\hat{\mathrm{\bf T}}}, is defined as
#'
#' \mjseqn{\hat{\mathrm{\bf T}} = 2\hat{\boldsymbol{\Theta}} - \hat{\boldsymbol{\Theta}}\hat{\mathrm{\bf R}}\hat{\boldsymbol{\Theta}},}
#'
#' where \mjseqn{\hat{\boldsymbol{\Theta}}} denotes the graphical lasso estimator of the precision matrix
#' and \mjseqn{\hat{\mathrm{\bf R}}} is the sample correlation matrix. Further details can be
#' found in Section 2 ("Main Results") of \insertCite{jankova2015confidence;textual}{GGMncv}.
#'
#' This approach is built upon earlier work on the de-sparsified lasso estimator
#' \insertCite{javanmard2014confidence,van2014asymptotically,zhang2014confidence}{GGMncv}
#'
#' @note
#' This assumes (reasonably) Gaussian data, and should not to be expected
#' to work for, say, polychoric correlations. Further, all work to date
#' has only looked at the graphical lasso estimator, and not de-sparsifying
#' nonconvex regularization. Accordingly, it is probably best to set
#' \code{penalty = "lasso"} in \code{\link{ggmncv}}.
#'
#' This function only provides the de-sparsified estimator and
#' not \emph{p}-values or confidence intervals (see \code{\link{inference}}).
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # data
#' Y <- GGMncv::Sachs[,1:5]
#'
#' n <- nrow(Y)
#' p <- ncol(Y)
#'
#' # fit model
#' # note: fix lambda, as in the reference
#' fit <- ggmncv(cor(Y), n = nrow(Y),
#'               progress = FALSE,
#'               penalty = "lasso",
#'               lambda = sqrt(log(p)/n))
#'
#' # fit model
#' # note: no regularization
#' fit_non_reg <- ggmncv(cor(Y), n = nrow(Y),
#'                       progress = FALSE,
#'                       penalty = "lasso",
#'                       lambda = 0)
#'
#'
#' # remove (some) bias and sparsity
#' That <- desparsify(fit)
#'
#' # graphical lasso estimator
#' fit$P
#'
#' # de-sparsified estimator
#' That$P
#'
#' # mle
#' fit_non_reg$P
#' @export
desparsify <- function(object, ...){

  if(!is(object, "ggmncv")){
    stop("object must be of class ggmncv")
  }

  p <- ncol(object$Theta)

  # Equation 7 in
  # Jankova, J., & Van De Geer, S. (2015). Confidence intervals for high-dimensional
  # inverse covariance estimation. Electronic Journal of Statistics, 9(1), 1205-1229.
  Theta <- 2 * object$Theta - object$Theta %*% object$R %*% object$Theta

  # partials
  P <- -(cov2cor(Theta) - diag(p))
  returned_object <- list(Theta = Theta, P = P)
  return(returned_object)

}
