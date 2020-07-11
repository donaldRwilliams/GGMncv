#' Compare Gaussian Graphical Models
#'
#' @description Compare Gaussian graphical models with the de-sparsified estimator of
#' \insertCite{jankova2015confidence}{GGMncv}.
#'
#' @param object_1 An object of class \code{ggmncv}
#'
#' @param object_2 An object of class \code{ggmncv}
#'
#' @param method Character string. A correction method for multiple comparison (defaults to \code{fdr}).
#' Can be abbreviated. See \link[stats]{p.adjust}.
#'
#' @param alpha Numeric. Significance level (defaults to \code{0.05}).
#'
#' @param ... Currently ignored.
#'
#' \itemize{
#'
#' \item \code{P_diff} De-sparsified partial correlation matrix differences
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
#' @examples
#' # data
#' Y1 <- MASS::mvrnorm(250, rep(0, 10), Sigma = diag(10))
#' Y2 <- MASS::mvrnorm(250, rep(0, 10), Sigma = diag(10))
#'
#' # fit models
#' fit1 <- GGMncv(Y1, n = nrow(Y1))
#' fit2 <- GGMncv(Y2, n = nrow(Y2))
#'
#' # compare
#' compare_ggms <- ggm_compare(fit1, fit2)
#' @export
ggm_compare <- function(object_1, object_2, method = "fdr", alpha = 0.05) {

  # columns
  p <- ncol(object_1$Theta)

  # corrected adjacency matrix
  adj_new <- corrected <- matrix(0, p, p)

  desparse_1 <- inference(object_1)
  desparse_2 <- inference(object_2)

  diff <- desparse_1$Theta - desparse_2$Theta

  sd_diff <-
    sqrt((desparse_1$sds / sqrt(desparse_1$n)) ^ 2 + (desparse_2$sds / sqrt(desparse_2$n)) ^
           2)

  z_diff <- diff / sd_diff


  p_values <- 2 * stats::pnorm(abs(z_diff), lower.tail = FALSE)

  # corrected p-values
  corrected_p_values <-
    stats::p.adjust(p_values[upper.tri(p_values)], method = method)
  corrected[upper.tri(corrected)] <- corrected_p_values
  corrected[lower.tri(corrected)] <-
    t(corrected)[lower.tri(corrected)]

  # new graph
  adj_new[upper.tri(adj_new)] <-
    ifelse(corrected_p_values < alpha, 1, 0)
  adj_new[lower.tri(adj_new)] <- t(adj_new)[lower.tri(adj_new)]

  pcor_diff <-
    (-cov2cor(desparse_1$Theta)) - (-cov2cor(desparse_2$Theta))

  returned_object <- list(
    P_diff = pcor_diff,
    adj = adj_new,
    corrected = corrected,
    uncorrected = p_values,
    method = method,
    alpha = alpha
  )

  class(returned_object) <- c("ggmncv", "ggm_compare")
  return(returned_object)

}


print_compare <- function(x, ... ){
  cat("GGM Compare \n")
  cat(paste0(x$method, ": ", x$alpha, "\n"))
  cat("---\n\n")
  adj <- as.data.frame( x$adj)
  colnames(adj) <- 1:ncol(adj)
  print(as.data.frame(adj))
}
