#' Regression Coefficients from \code{ggmncv} Objects
#'
#' @param object An Object of classs \code{ggmncv}
#' @param ... Currently ignored
#'
#' @return A matrix of regression coefficients
#'
#' @examples
#'
#' \donttest{
#'
#' # data
#' Y <- GGMncv::ptsd
#'
#' # correlations
#' S <- cor(Y)
#'
#' # fit model
#' fit <- GGMncv(S, n = nrow(Y))
#'
#' coefs <- coef(fit)
#'
#' }
#' @export
coef.ggmncv <- function(object,...){

  Theta <- object$Theta
  coefs <- coef_helper(Theta)
  return(coefs)
}
