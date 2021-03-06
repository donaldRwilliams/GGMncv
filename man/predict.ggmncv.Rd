% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.ggmncv}
\alias{predict.ggmncv}
\title{Predict method for \code{ggmncv} Objects}
\usage{
\method{predict}{ggmncv}(object, train_data, newdata = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{ggmncv}}

\item{train_data}{Data used for model fitting.}

\item{newdata}{An optional data frame in which to look for variables with which to predict.
If omitted, the fitted values are used.}

\item{...}{Currently ignored}
}
\value{
A matrix of predicted values
}
\description{
Predicted values based on a \code{ggmncv} object
}
\examples{
# data
Y <- scale(Sachs)

# test data
Ytest <- Y[1:100,]

# training data
Ytrain <- Y[101:nrow(Y),]

fit <- ggmncv(cor(Ytrain), n = nrow(Ytrain))

pred <- predict(fit, train_data = Y,
                newdata = Ytest)

round(apply((pred - Ytest)^2, 2, mean), 2)
}
