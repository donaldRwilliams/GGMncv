#' Predict method for \code{ggmncv} Objects
#'
#' @description Predicted values based on a \code{ggmncv} object
#'
#' @param object An object of class \code{ggmncv}
#'
#' @param newdata  An optional data frame in which to look for variables with which to predict.
#'                 If omitted, the fitted values are used.
#'
#' @param ... Currently ignored
#'
#' @return A matrix of predicted values
#'
#' @examples
#' # data
#' Y <- scale(Sachs)
#'
#' # test data
#' Ytest <- Y[1:100,]
#'
#' # training data
#' Ytrain <- Y[101:nrow(Y),]
#'
#' fit <- GGMncv(Ytrain, n = nrow(Ytrain))
#'
#' pred <- predict(fit,
#'                 newdata = Ytest)
#'
#' round(apply((pred - Ytest)^2, 2, mean), 2)
#' @export
predict.ggmncv <- function(object, newdata = NULL,...){

  if(isSymmetric(as.matrix(object$x))){
    stop("data matrix not found")
  }
  if(!is.null(newdata)){
    x_scale <-  newdata
  } else {
    x_scale <- scale(object$x)
  }
  coefs <- unclass(coef(object))
  p <- ncol(x_scale)
  y_hat <- sapply(1:p, function(x)  x_scale[,-x] %*% coefs[x,])
  return(y_hat)
}
