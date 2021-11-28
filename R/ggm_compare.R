#' Compare Edges Between Gaussian Graphical Models
#'
#' @description Establish whether each of the corresponding edges
#'              are significantly different in two groups,
#'              with the de-sparsified estimator of \insertCite{jankova2015confidence}{GGMncv}.
#'
#' @param  object_1 object of class \code{\link{ggmncv}} .
#'
#' @param object_2 An object of class \code{\link{ggmncv}}.
#'
#' @param method Character string. A correction method for
#'               multiple comparisons (defaults to \code{fdr}), which
#'               can be abbreviated. See \link[stats]{p.adjust}.
#'
#' @param alpha Numeric. Significance level (defaults to \code{0.05}).
#'
#' @param ... Currently ignored.
#'
#' @references
#' \insertAllCited{}
#'
#' @note For low-dimensional settings, i.e., when the number of observations
#' far exceeds the number of nodes, this function likely has limited utility and
#' a non regularized approach should be used for comparing edges
#' (see for example \strong{GGMnonreg}).
#'
#' Further, whether the de-sparsified estimator provides nominal error rates
#' remains to be seen, at least across a range of conditions. For example,
#' the simulation results in \insertCite{williams_2021;textual}{GGMncv}
#' demonstrated that the confidence intervals
#' can have (severely) compromised coverage properties (whereas non-regularized methods
#' had coverage at the nominal level).
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{P_diff} De-sparsified partial correlation differences
#'
#' \item \code{adj} Adjacency matrix based on the p-values.
#'
#' \item \code{pval_uncorrected} Uncorrected p-values
#'
#' \item \code{pval_corrected} Corrected p-values
#'
#' \item \code{method} The approach used for multiple comparisons
#'
#' \item \code{alpha} Significance level
#' }
#'
#' @examples
#' # data
#' # note: all edges equal
#' Y1 <- MASS::mvrnorm(250, rep(0, 10), Sigma = diag(10))
#' Y2 <- MASS::mvrnorm(250, rep(0, 10), Sigma = diag(10))
#'
#' # fit models
#' # note: atan penalty by default
#'
#' # group 1
#' fit1 <- ggmncv(cor(Y1), n = nrow(Y1),
#'                progress = FALSE)
#'
#' # group 2
#' fit2 <- ggmncv(cor(Y2), n = nrow(Y2),
#'                progress = FALSE)
#'
#' # compare
#' compare_ggms <- compare_edges(fit1, fit2)
#'
#' compare_ggms
#' @export
compare_edges <- function(object_1,
                          object_2,
                          method = "fdr",
                          alpha = 0.05,
                          ...) {

  # number of nodes
  p <- ncol(object_1$Theta)

  # corrected adjacency matrix
  adj_new <- corrected <- matrix(0, p, p)

  # de-sparsify
  desparse_1 <- inference(object_1)
  desparse_2 <- inference(object_2)

  # difference in precision matrices
  diff <- desparse_1$Theta - desparse_2$Theta

  # standard error of difference
  sd_diff <- sqrt((desparse_1$sds / sqrt(desparse_1$n)) ^ 2 +
                  (desparse_2$sds / sqrt(desparse_2$n)) ^ 2)

  # z score
  z_diff <- diff / sd_diff

  # two sided p-values
  p_values <- 2 * stats::pnorm(abs(z_diff), lower.tail = FALSE)

  # corrected p-values
  corrected_p_values <- stats::p.adjust(p_values[upper.tri(p_values)],
                                        method = method)

  corrected[upper.tri(corrected)] <- corrected_p_values

  # make symmetric
  corrected[lower.tri(corrected)] <- t(corrected)[lower.tri(corrected)]

  # new graph
  adj_new[upper.tri(adj_new)] <- ifelse(corrected_p_values < alpha, 1, 0)
  adj_new[lower.tri(adj_new)] <- t(adj_new)[lower.tri(adj_new)]

  pcor_diff <-
    (-cov2cor(desparse_1$Theta)) - (-cov2cor(desparse_2$Theta))

  returned_object <- list(
    P_diff = pcor_diff,
    adj = adj_new,
    pvals_corrected = corrected,
    p_vals_uncorrected = p_values,
    method = method,
    alpha = alpha
  )

  class(returned_object) <- c("ggmncv",
                              "ggm_compare")
  return(returned_object)

}


print_compare <- function(x, ... ){

  cat("Compare Edges \n")

  cat(paste0(x$method, ": ", x$alpha, "\n"))

  cat("---\n\n")

  adj <- as.data.frame( x$adj)

  colnames(adj) <- 1:ncol(adj)

  print(as.data.frame(adj))
}
