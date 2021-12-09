#' Statistical Inference for Regularized Gaussian Graphical Models
#'
#' @name inference
#'
#' @description Compute \emph{p}-values for each relation based on the
#' de-sparsified glasso estimator \insertCite{jankova2015confidence}{GGMncv}.
#'
#' @param object  An object of class \code{ggmncv}
#'
#' @param method Character string. A correction method for multiple comparison (defaults to \code{fdr}).
#' Can be abbreviated. See \link[stats]{p.adjust}.
#'
#' @param alpha Numeric. Significance level (defaults to \code{0.05}).
#'
#' @param ... Currently ignored.
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{Theta} De-sparsified precision matrix
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
#' @importFrom stats p.adjust pnorm
#'
#' @note
#' This assumes (reasonably) Gaussian data, and should not to be expected
#' to work for, say, polychoric correlations. Further, all work to date
#' has only looked at the graphical lasso estimator, and not de-sparsifying
#' nonconvex regularization. Accordingly, it is probably best to set
#' \code{penalty = "lasso"} in \code{\link{ggmncv}}.
#'
#' Further, whether the de-sparsified estimator provides nominal error rates
#' remains to be seen, at least across a range of conditions. For example,
#' the simulation results in \insertCite{williams_2021;textual}{GGMncv}
#' demonstrated that the confidence intervals
#' can have (severely) compromised coverage properties (whereas non-regularized methods
#' had coverage at the nominal level).
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd[,1:5]
#'
#' # fit model
#' fit <- ggmncv(cor(Y), n = nrow(Y),
#'               progress = FALSE,
#'               penalty = "lasso")
#'
#'
#' # statistical inference
#' inference(fit)
#' @export
inference <- function(object,
                      method = "fdr",
                      alpha = 0.05,
                      ...){
  if(!is(object, "ggmncv")){
    stop("object must be of class ggmncv")
  }
  # columns
  p <- ncol(object$Theta)

  # corrected adjacency matrix
  adj_new <- corrected <- matrix(0, p, p)

  # observations
  n <- object$n

  # precision matrix (sparsified)
  Theta <- object$Theta

  # sd
  sds <- sqrt((tcrossprod(diag(Theta)) + Theta^2))

  # desparsify
  Theta <- desparsify(object)$Theta

  # z stats
  z_stat <- sapply(1:p, function(x){
    Theta[x,]/(sds[x,] /sqrt(n))
  })

  # p values
  p_values <- 2 * stats::pnorm( abs(z_stat), lower.tail = FALSE)

  # corrected p-values
  corrected_p_values <- stats::p.adjust(p_values[upper.tri(p_values)], method = method)
  corrected[upper.tri(corrected)] <- corrected_p_values
  corrected[lower.tri(corrected)] <- t(corrected)[lower.tri(corrected)]

  # new graph
  adj_new[upper.tri(adj_new)] <-  ifelse(corrected_p_values < alpha, 1, 0)
  adj_new[lower.tri(adj_new)] <- t(adj_new)[lower.tri(adj_new)]

  P <- -(cov2cor(Theta) - diag(p))
  P <- P * adj_new
  # return object
  returned_object <- list(Theta = Theta,
                          P = P,
                          adj = adj_new,
                          pval_uncorrect = p_values,
                          p_val_corrected = corrected,
                          method = method,
                          alpha = alpha,
                          sds = sds, n = n)

  class(returned_object) <- c("ggmncv",
                              "inference")
  return(returned_object)
}


#' @rdname inference
#' @examples
#'
#' # alias
#' all.equal(inference(fit), significance_test(fit))
#'
#' @export
significance_test <- inference

print_inference <- function(x, ...){

  cat("Statistical Inference\n")

  cat(paste0(x$method, ": ", x$alpha, "\n"))

  cat("---\n\n")

  adj <- as.data.frame(x$adj)

  colnames(adj) <- 1:ncol(adj)

  print(as.data.frame(adj))
}
