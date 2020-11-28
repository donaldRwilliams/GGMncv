#' Confirm Edges
#'
#' @description Confirmatory testing of edges detected with data-driven
#'              model selection in an independent dataset.
#'
#' @param object An object of classo \code{ggmncv}
#'
#' @param Rnew  Matrix. A correlation matrix of dimensions \emph{p} by \emph{p}.
#'
#' @param method Character string. A correction method for multiple comparison (defaults to \code{fdr}).
#' Can be abbreviated. See \link[stats]{p.adjust}.
#'
#' @param alpha Numeric. Significance level (defaults to \code{0.05}).
#'
#' @return An object of class \code{ggmncv}
#'
#' @export
#'
#' @examples
#' Y <- na.omit(bfi[,1:25])
#'
#' Y_explore <- Y[1:1000,]
#'
#' Y_confirm <- Y[1001:nrow(Y),]
#'
#' fit <- ggmncv(cor(Y_explore),
#'               n = nrow(Y_explore))
#'
#' confirm <- confirm_edges(fit,
#'                          Rnew = cor(Y_confirm),
#'                          method = "fdr",
#'                          alpha = 0.05)
confirm_edges <- function(object, Rnew, method, alpha) {

  if (!is(object, "ggmncv")) {
    stop("must be a ggmncv object.")
  }

  # fitted model
  fit <- object

  # nodes
  p <- ncol(fit$Theta)

  # estimate new R with lasso
  fitnew <-
    ggmncv(
      Rnew,
      n = fit$n,
      penalty = "lasso",
      lambda = sqrt(log(p) / fit$n),
      progress = FALSE
    )

  # debias
  inf <- inference(fitnew)

  confirm_which <- which(fit$adj[upper.tri(diag(p))] == 1)

  # only test those in the object
  ps <- p.adjust(p = inf$uncorrect[upper.tri(diag(p))][confirm_which],
                 method = method)

  confirm_mat <- matrix(0, p, p)

  confirm_mat[upper.tri(diag(p))][confirm_which] <- ifelse(ps < alpha, 1, 0)

  confirm_mat <- symmetric_mat(confirm_mat)

  P <- -(cov2cor(inf$Theta) - diag(p))

  P_confirm <- confirm_mat * P

  returned_object <- list(P = P_confirm, adj = confirm_mat)
  class(returned_object) <- c("ggmncv", "default")
  return(returned_object)
}
