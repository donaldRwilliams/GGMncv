#' Constrained Precision Matrix
#'
#' @description Compute the maximum likelihood estimate, given certain elements are constrained to zero
#' (e.g., an adjacency matrix). This approach is described in \insertCite{hastie2015statistical;textual}{GGMncv}.
#'
#' @param Sigma Covariance matrix
#'
#' @param adj Matrix with contraints. A zero indicates that element
#'            should be constrained to zero.
#'
#' @return  A list containing the inverse covariance matrix and the covariance matrix.
#'
#' @note The algorithm is written in \code{c++}.
#'
#' @examples
#' \donttest{
#' # data
#' Y <- GGMncv::ptsd[,1:5]
#'
#' # columns
#' p <- ncol(Y)
#'
#' # contstraint matrix
#' constraints <- matrix(0,p,p)
#'
#' # set one value to zero
#' constraints[2,3] <- 1
#' constraints[3,2] <-1
#'
#' # estimate, given constraints
#' fit <- constrained(cor(Y), adj = constraints)
#' Theta <- fit$Theta
#' }
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

  returned_object <- list(Theta = Theta, Sigma = Sigma)

  return(returned_object)
}
