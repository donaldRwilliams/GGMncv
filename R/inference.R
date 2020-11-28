#' Statistical Inference for Gaussian Graphical Models
#'
#' @description Compute p-values for each relation based on the de-sparsified precision matrix
#' \insertCite{jankova2015confidence}{GGMncv}.
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
#' \item \code{uncorrected} Uncorrected p-values
#'
#' \item \code{corrected} Corected p-values
#'
#' \item \code{method} The approach used for multiple comparisons
#'
#' \item \code{alpha} Significance level
#' }
#'
#' @importFrom stats p.adjust pnorm
#'
#' @note
#' This assumes the Gaussian data.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd
#'
#' # fit model
#' fit <- ggmncv(cor(Y), n = nrow(Y))
#'
#'
#' # statistical inference
#' inference(fit)
#'
#'
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
                          uncorrect = p_values,
                          corrected = corrected,
                          method = method,
                          alpha = alpha,
                          sds = sds, n = n)

  class(returned_object) <- c("ggmncv",
                              "inference")
  return(returned_object)
}


print_inference <- function(x, ... ){
  cat("Statistical Inference\n")
  cat(paste0(x$method, ": ", x$alpha, "\n"))
  cat("---\n\n")
  adj <- as.data.frame( x$adj)
  colnames(adj) <- 1:ncol(adj)
  print(as.data.frame(adj))
}
