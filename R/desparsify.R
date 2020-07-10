#' De-sparsified Graphical Lasso Estimator
#'
#' @description Compute the desparsified glasso estimator with the approach
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
#' \item \code{Theta} De-sparsified precision matrix
#'
#' \item \code{P} De-sparsified partial correlation matrix
#'
#' }
#'
#' @note
#' This assumes the Gaussian data.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd
#'
#' # fit model
#' fit <- GGMncv::GGMncv(cor(Y), n = nrow(Y))
#'
#' desparsify(fit)
#'
#' @export
desparsify <- function(object, ...){
  if(!is(object, "ggmncv")){
    stop("object must be of class ggmncv")
  }
  p <- ncol(object$Theta)
  # Equation 7 in
  # Jankova, J., & Van De Geer, S. (2015). Confidence intervals for high-dimensional
  # inverse covariance estimation. Electronic Journal of Statistics, 9(1), 1205-1229.
  Theta <- 2 * object$Theta - object$Theta%*% object$R %*% object$Theta

  # partials
  P <- -(cov2cor(Theta) - diag(p))
  returned_object <- list(Theta = Theta, P = P)
  return(returned_object)

}
