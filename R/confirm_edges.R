#' @title Confirm Edges
#'
#' @description Confirmatory hypothesis testing of edges that were initially
#'              detected with data-driven model selection.
#'
#' @param object An object of class \code{ggmncv}
#'
#' @param Rnew  Matrix. A correlation matrix of dimensions \emph{p} by \emph{p}.
#'
#' @param method Character string. A correction method for multiple comparison
#'               (defaults to \code{fdr}). Can be abbreviated. See \link[stats]{p.adjust}.
#'
#' @param alpha Numeric. Significance level (defaults to \code{0.05}).
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{ggmncv}, including:
#'
#' \itemize{
#'
#' \item{\strong{P}}: Matrix of confirmed edges (partial correlations)
#'
#' \item{\strong{adj}}: Matrix of confirmed edges (adjacency)
#'
#' }
#'
#'
#' @details
#' The basic idea is to merge exploration with confirmation
#' \insertCite{@see for example,  @rodriguez_williams_rast_mulder_2020}{GGMncv}.
#' This is accomplished by testing those edges (in an independent dataset)
#' that were initially detected via data driven model selection.
#'
#' Confirmatory hypothesis testing follows the approach described in
#' \insertCite{jankova2015confidence;textual}{GGMncv}: (1)
#' graphical lasso is computed with lambda fixed to  \mjseqn{\lambda = \sqrt{log(p)/n}},
#' (2) the de-sparsified estimator is computed, and then (3) \emph{p}-values are
#' obtained for the de-sparsified estimator.
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
#'               n = nrow(Y_explore),
#'               progress = FALSE)
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

  returned_object <- list(P = P_confirm,
                          adj = confirm_mat)

  class(returned_object) <- c("ggmncv", "default")
  return(returned_object)
}
