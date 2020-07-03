#' Regression Coefficients from \code{ggmncv} Objects
#'
#' @param object An Object of classs \code{ggmncv}
#'
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
coef.ggmncv <- function(object, ...){

  Theta <- object$Theta
  coefs <- coef_helper(Theta)
  class(coefs) <- c("ggmncv", "coef")
  return(coefs)
}


print_coef <- function(x,...){

  p <- ncol(x)
  cat("Estimates:\n\n")
  for(i in 1:p){
    cat(paste("node", i, "\n\n"))
    nodes_id <-  (1:20)[-i]
    dat <-  as.data.frame(t(cs[i,]))
    colnames( dat) <- paste("node", nodes_id)
    print(dat, row.names = FALSE  )
    cat("---\n")
    }
}
