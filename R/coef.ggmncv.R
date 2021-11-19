#' Regression Coefficients from \code{ggmncv} Objects
#'
#' @param object An Object of class \code{ggmncv}.
#'
#' @param ... Currently ignored.
#'
#' @return A matrix of regression coefficients.
#'
#' @note
#' The matrix of coefficients can be accessed via \code{coefs[1,]},
#' which provides the estimates for predicting the first node.
#'
#' Further, the estimates are essentially computed with both
#' the outcome and predictors scaled to have mean 0 and
#' standard deviation 1.
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
    print(dat, row.names = FALSE)
    cat("---\n\n")
  }
}
