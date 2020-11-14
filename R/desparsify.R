#' De-Sparsified Graphical Lasso Estimator
#'
#' @description
#' \loadmathjax
#' Compute the de-sparsified glasso estimator with the approach
#' described in Equation 7 of \insertCite{jankova2015confidence;textual}{GGMncv}.
#'
#' @param object An object of class \code{ggmncv}
#'
#' @param ... Currently ignored
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
#' found in section 2 ("Main Results") of \insertCite{jankova2015confidence;textual}{GGMncv}.
#'
#'
#' @note
#' This assumes (reasonably) Gaussian data.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # data
#' Y <- GGMncv::Sachs
#'
#' # fit model
#' fit <- ggmncv(cor(Y), n = nrow(Y))
#'
#' # remove (some) bias and sparsity
#' That <- desparsify(fit)
#'
#' # de-sparsified partial correlations
#' That$P
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
