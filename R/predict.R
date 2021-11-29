#' Predict method for \code{ggmncv} Objects
#'
#' @description There is a direct correspondence between the inverse covariance
#' matrix and multiple regression \insertCite{stephens1998,kwan2014regression}{GGMncv}.
#' This readily allows for converting the off diagonal elements to regression coefficients,
#' opening the door to out-of-sample prediction in multiple regression.
#'
#' @param object An object of class \code{\link{ggmncv}}.
#'
#' @param train_data Data used for model fitting (defaults to \code{NULL}).
#'
#' @param newdata  An optional data frame in which to look for variables with which to predict.
#'                 If omitted, the fitted values are used.
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix of predicted values, of dimensions rows
#' (in the training/test data) by the number of nodes (columns).
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
#' fit <- ggmncv(cor(Ytrain), n = nrow(Ytrain),
#'               progress = FALSE)
#'
#' pred <- predict(fit, newdata = Ytest)
#'
#' round(apply((pred - Ytest)^2, 2, mean), 2)
#' @export
predict.ggmncv <- function(object, train_data = NULL,
                           newdata = NULL, ...){

  if(!is.null(newdata)){

    x_scale <-  as.matrix(newdata)

  } else {

    if(is.null(train_data)){
      stop("`train_data` required if `newdata` is NULL")
    }

    object$x <- train_data

    x_scale <- as.matrix(scale(object$x))

  }

  coefs <- as.matrix(unclass(coef(object)))

  p <- ncol(x_scale)

  y_hat <- sapply(1:p, function(x)  x_scale[,-x] %*% coefs[x,])

  return(y_hat)
}
