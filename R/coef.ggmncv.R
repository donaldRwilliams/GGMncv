#' Regression Coefficients from \code{ggmncv} Objects
#'
#' @description There is a direct correspondence between the inverse covariance
#' matrix and multiple regression \insertCite{stephens1998,kwan2014regression}{GGMncv}.
#' This readily allows for converting the off diagonal elements to regression coefficients,
#' resulting in noncovex penalization for multiple regression modeling.
#'
#'
#' @param object An Object of class \code{ggmncv}.
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#' @return A matrix of regression coefficients.
#'
#' @note
#' The coefficients can be accessed via \code{coefs[1,]},
#' which provides the estimates for predicting the first node.
#'
#' Further, the estimates are essentially computed with both
#' the outcome and predictors scaled to have mean 0 and
#' standard deviation 1.
#'
#'
#' @examples
#'
#' \donttest{
#'
#' # data
#' Y <- GGMncv::ptsd[,1:5]
#'
#' # correlations
#' S <- cor(Y)
#'
#' # fit model
#' fit <- ggmncv(R = S, n = nrow(Y), progress = FALSE)
#'
#' # regression
#' coefs <- coef(fit)
#'
#' coefs
#'
#'
#' # no regularization, resulting in OLS
#'
#' # data
#' # note: scaled for lm()
#' Y <- scale(GGMncv::ptsd[,1:5])
#'
#' # correlations
#' S <- cor(Y)
#'
#' # fit model
#' # note: non reg
#' fit <- ggmncv(R = S, n = nrow(Y), progress = FALSE, lambda = 0)
#'
#' # regression
#' coefs <- coef(fit)
#'
#' # fit lm
#' fit_lm <- lm(Y[,1] ~ 0 + Y[,-1])
#'
#' # ggmncv
#' coefs[1,]
#'
#' # lm
#' as.numeric(coef(fit_lm))
#'
#' }
#'
#' @export
coef.ggmncv <- function(object, ...) {

  # precision matrix
  Theta <- object$Theta

  # inverse to regression
  coefs <- coef_helper(Theta)

  class(coefs) <- c("ggmncv", "coef")

  return(coefs)
}


print_coef <- function(x, ...) {

  p <- nrow(x)

  cat("Estimates:\n\n")

  for (i in seq_len(p)) {

    cat(paste0("node.", i, "\n"))

    nodes_id <-  (1:p)[-i]

    dat <-  as.data.frame(t(x[i, ]))

    colnames(dat) <- paste0("node.", nodes_id)

    print(round(dat, 3), row.names = FALSE)

    cat("---\n\n")

  }
}
