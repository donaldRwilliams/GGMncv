#' Precision Matrix with Known Graph
#'
#' @description Compute the maximum likelihood estimate of the precision matrix,
#' given a known graphical structure (i.e., an adjacency matrix).
#' This approach was originally described in "The Elements of Statistical Learning"
#' \insertCite{@see pg. 631, @hastie2009elements}{GGMncv}.
#'
#' @param Sigma Covariance matrix
#'
#' @param adj Adjacency matrix that encodes the constraints,
#'            where a zero indicates that element should be zero.
#'
#' @return  A list containing the following:
#'
#' \itemize{
#'
#' \item{\strong{Theta}}: Inverse of the covariance matrix (precision matrix)
#'
#' \item{\strong{Sigma}}: Covariance matrix.
#'
#' \item{\strong{wadj}}: Weighted adjacency matrix, corresponding
#' to the partial correlation network.
#'
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @note
#' The algorithm is written in \code{c++}, and should scale to high dimensions
#' nicely.
#'
#' Note there are a variety of algorithms for this purpose. Simulation
#' studies indicated that this approach is both accurate and computationally
#' efficient \insertCite{@HFT therein, @emmert2019constrained}{GGMncv}
#'
#' @examples
#' \donttest{
#'
#' # data
#' y <- ptsd
#'
#' # fit model
#' fit <- ggmncv(cor(y), n = nrow(y),
#'               penalty = "lasso",
#'               progress = FALSE)
#'
#' # set negatives to zero (sign restriction)
#' adj_new <- ifelse( fit$P <= 0, 0, 1)
#'
#' check_zeros <- TRUE
#'
#' # track trys
#' iter <- 0
#'
#' # iterate until all positive
#' while(check_zeros){
#'   iter <- iter + 1
#'   fit_new <- constrained(cor(y), adj = adj_new)
#'   check_zeros <- any(fit_new$wadj < 0)
#'   adj_new <- ifelse( fit_new$wadj <= 0, 0, 1)
#' }
#'
#'}
#' @export
constrained <- function(Sigma, adj){
  # change to zeros
  # adj <- ifelse(adj == 1, 0, 1)

  # include diagonal!
  diag(adj) <- 1

  # call c++
  fit <- hft_algorithm(Sigma, adj, tol = 1e-10, max_iter = 100)

  Theta <- round(fit$Theta, 3)

  Sigma <- round(fit$Sigma, 3)

  # weighted adjacency matrix
  # partial correlation network
  wadj <- -(cov2cor(Theta) - diag(ncol(adj)))

  returned_object <- list(Theta = Theta,
                          Sigma = Sigma,
                          wadj = wadj)

  return(returned_object)
}
