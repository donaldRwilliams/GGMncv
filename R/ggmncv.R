#' GGMncv
#'
#' @param x There are 2 options: either a \code{n} by \code{p} data matrix or a
#'          \code{p} by \code{p} correlation matrix.
#'
#' @param n Numeric. Sample size.
#'
#' @param penalty Character string. Which penalty should be used (defaults to \code{atan}).
#'
#' @param nonreg Logical. Should the precision matrix be reestimated, given the adjacency matrix
#'               (defaults to \code{FALSE})? When set to \code{TRUE}, this provides the \strong{nonregularized},
#'               maximum likelhood estimate with constraints.
#'
#' @param LLA Logical. Should the local linear approximation be used for maxmizing the penalized likelihood ?
#'            The default is \code{FALSE} (see details).
#'
#' @param method Character string. Which correlation coefficient should be computed.
#'               One of "pearson" (default), "spearman", or  "polychoric".
#'
#' @param vip Logical. Should the respective variable inclusion 'probabilities' be
#'            computed (defaults to \code{FALSE})? This is accomplished with a
#'            non-parametric bootstrap (see details)
#'
#' @param vip_iter Numeric. How many number of bootstrap samples (defaults to 1000)?
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ? Note that
#'                 this only applies when \code{vip = TRUE}.
#'
#' @param ... Currently ignored.
#'
#' @importFrom stats cor cov2cor
#' @importFrom psych polychoric
#' @importFrom glassoFast glassoFast
#'
#' @return An object of class \code{ggmncv}, including:
#'
#' \itemize{
#' \item \code{Theta} Inverse covariance matrix
#' \item \code{Sigma} Covariance matrix
#' \item \code{P} Weigthed adjacency matrix
#' \item \code{adj} Adjacency matrix
#' \item \code{lambda} Tuning parameter (i.e., sqrt(log(p)/n))
#' \item \code{fit} glasso fitted model (a list)
#' }
#'
#'
#' @examples
#' # data
#' Y <- BGGM::ptsd
#'
#' # polychoric
#' S <- psych::polychoric(Y)$rho
#'
#' # fit model
#' fit <- GGMncv(S, n = nrow(Y))
#'
#' qgraph::qgraph(fit$P)
#'
#' @export
GGMncv <- function(x, n,
                   penalty = "atan",
                   nonreg = FALSE,
                   LLA = FALSE,
                   method = "pearson",
                   vip = FALSE,
                   vip_iter = 1000,
                   progress = TRUE,
                   ...){

  if(! penalty %in% c("atan", "mcp", "scad", "exp", "selo")){
    stop("penalty not found. \ncurrent options: atan, mcp, scad, exp, selo")
  }

  if (base::isSymmetric(x)) {
    R <- x
  } else {
    if (method == "polychoric") {
    suppressWarnings(
      R <- psych::polychoric(x)$rho
    )
    } else {
      R <- stats::cor(x, method = method)
      }
  }

  # nodes
  p <- ncol(R)

  # identity matrix
  I_p <- diag(p)

  # tuning
  lambda <- sqrt(log(p)/n)

  # inverse covariance matrix
  Theta <- solve(R)

  if(!LLA){
  # lambda matrix
  lambda_mat <-
    eval(parse(text =  paste0(
      penalty, "_deriv(Theta = Theta, lambda = lambda)"
    )))

  diag(lambda_mat) <- lambda

  fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)


  } else {

    Theta <- glassoFast::glassoFast(S = R, rho = lambda)$wi

    lambda_mat <-
      eval(parse(text =  paste0(
        penalty, "_deriv(Theta = Theta, lambda = lambda)"
      )))

    diag(lambda_mat) <- lambda

    fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

   }

  adj <- ifelse(fit$wi == 0, 0, 1)

  if(nonreg) {
    rest <- htf(R, adj)
    Theta <- rest$Theta
    Sigma <- rest$Sigma
    P <- -(stats::cov2cor(Theta) - I_p)
  } else {
    Theta <- fit$wi
    Sigma <- fit$w
    P <- -(stats::cov2cor(Theta) - I_p)
  }

  if(vip){

    if(progress){
      pb <- utils::txtProgressBar(min = 0, max = vip_iter, style = 3)
    }

    vip_results <-sapply(1:vip_iter, function(i){

      Yboot <- Y[sample(1:n, size = n, replace = TRUE),]

      if(method == "polychoric"){

        suppressWarnings(
          R <- psych::polychoric(Yboot)$rho
          )

        } else {

        R <- cor(Yboot, method = method)
      }

      Theta <- solve(R)

      lambda_mat <-
        eval(parse(text =  paste0(
          penalty, "_deriv(Theta = Theta, lambda = lambda)"
        )))

      diag(lambda_mat) <- lambda
      fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)
      adj <- ifelse(fit$wi == 0, 0, 1)

      if(progress){
       utils::setTxtProgressBar(pb, i)
      }

      adj[upper.tri(adj)]
    })

    if(is.null( colnames(Y))){
      cn <- 1:p
    } else {

      cn <- colnames(Y)
    }

    vip_results <-
      data.frame(Relation =  sapply(1:p, function(x)
        paste0(cn, "--", x))[upper.tri(I_p)],
        VIP = rowMeans(vip_results))

  } else {

    vip_results <- NULL
  }

  returned_object <- list(Theta = Theta,
                          Sigma = Sigma,
                          P = P,
                          fit = fit,
                          adj = adj,
                          lambda = lambda,
                          vip_results = vip_results)

  class(returned_object) <- "ggmncv"
  return(returned_object)
}

#' Print \code{ggmncv} Objects
#'
#' @param x An object of class \code{ggmncv}
#' @param ... Currently ignored
#' @export
print.ggmncv <- function(x, ...){
  mat <- round(x$P, 3)
  colnames(mat) <- 1:ncol(x$P)
  rownames(mat) <- 1:ncol(x$P)
  print(mat)
}



#' Plot \code{ggmncv} Objects
#'
#' @description Plot variable inclusion 'probabilities'
#'
#' @param x An object of class \code{ggmncv}
#'
#' @param size Numeric. The size of the points (defaults to 1).
#'
#' @param color Character string. The color of the points (defaults to \code{black})
#'
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object
#'
#' @import ggplot2
#'
#'
#' @examples
#' # data
#' Y <- BGGM::ptsd
#'
#' # polychoric
#' S <- psych::polychoric(Y)$rho
#'
#' # fit model
#' fit <- GGMncv(S, n = nrow(Y),
#'               nonreg = FALSE,
#'               vip = TRUE,
#'               progress = FALSE)
#'
#' # plot VIP
#' plot(fit)
#' @export
plot.ggmncv <- function(x,
                        size = 1,
                        color = "black",
                        ...){
  if(is.null(x$vip_results)){
    stop("variable inclusion 'probabilities' not found (set vip = TRUE)")
  }

  dat <- x$vip_results[order(x$vip_results$VIP),]

  dat$new1 <- factor(dat$Relation,
                     levels = dat$Relation,
                     labels = dat$Relation)

  ggplot(dat,aes(y= new1,
                 x = VIP,
                 group = new1)) +
    geom_point(size = size,
               color = color)  +
    ylab("Relation")

}
