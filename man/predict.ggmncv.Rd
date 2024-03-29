% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.ggmncv}
\alias{predict.ggmncv}
\title{Predict method for \code{ggmncv} Objects}
\usage{
\method{predict}{ggmncv}(object, train_data = NULL, newdata = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{\link{ggmncv}}.}

\item{train_data}{Data used for model fitting (defaults to \code{NULL}).}

\item{newdata}{An optional data frame in which to look for variables with which to predict.
If omitted, the fitted values are used.}

\item{...}{Currently ignored.}
}
\value{
A matrix of predicted values, of dimensions rows
(in the training/test data) by the number of nodes (columns).
}
\description{
There is a direct correspondence between the inverse covariance
matrix and multiple regression \insertCite{stephens1998,kwan2014regression}{GGMncv}.
This readily allows for converting the off diagonal elements to regression coefficients,
opening the door to out-of-sample prediction in multiple regression.
}
\examples{
# data
Y <- scale(Sachs)

# test data
Ytest <- Y[1:100,]

# training data
Ytrain <- Y[101:nrow(Y),]

fit <- ggmncv(cor(Ytrain), n = nrow(Ytrain),
              progress = FALSE)

pred <- predict(fit, newdata = Ytest)

round(apply((pred - Ytest)^2, 2, mean), 2)
}
\references{
\insertAllCited{}
}
