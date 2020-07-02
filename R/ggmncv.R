#' GGMncv
#'
#' @description Estimate Gaussian graphical models with non-convex penalties.
#'
#' @param x There are 2 options: either a \code{n} by \code{p} data matrix or a
#'          \code{p} by \code{p} correlation matrix.
#'
#' @param n Numeric. Sample size.
#'
#' @param lambda Numeric. Tuning parameter governing the degress of penalization. Defaults to
#'               \code{NULL} which results in fixing lambda to \code{sqrt(log(p)/n)}.
#'
#' @param n_lambda Numeric. The number of regularization/thresholding parameters
#'                 (defaults to \code{100}).
#'
#' @param gamma Numeric. Hyperparameter for the penalty function. Defaults to \code{0.1}
#'              which is the recommended value for the default penalty (see details).
#'
#' @param select Logical. Should lambda be selected with BIC (defaults to \code{FALSE})?
#'
#' @param L0_learn Logical. Should lambda be selected based on the non-regularized
#'                 precision matrix (defaults to \code{FALSE}; see details).
#'
#' @param penalty Character string. Which penalty should be used (defaults to \code{atan}).
#'
#' @param refit Logical. Should the precision matrix be refitted, given the adjacency matrix
#'               (defaults to \code{FALSE})? When set to \code{TRUE}, this provides the \strong{nonregularized},
#'               maximum likelhood estimate with constraints.
#'
#' @param LLA Logical. Should the local linear approximation be used for maxmizing the penalized likelihood ?
#'            The default is \code{FALSE} (see details).
#'
#' @param method Character string. Which correlation coefficient should be computed.
#'               One of "pearson" (default), "spearman", or  "polychoric".
#'
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ? Note that
#'                 this only applies when \code{select = TRUE}.
#'
#' @param store Logical. Should all of the fitted models be saved (defaults to \code{NULL}). Note
#'              this only applies when \code{select = TRUE}. and ignored otherwise
#'              (the one model is saved.)
#'
#' @param vip Logical. Should variable inclusion "probabilities" be computed (defaults to \code{FALSE})?
#'
#' @param vip_iter Numeric. How many bootstrap sample for computing \code{vip} (defaults to 1000) ? Note
#'                 also that only the default lambda is implemented (select is not implemented).
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
#' \donttest{
#' # data
#' Y <- GGMncv::ptsd
#'
#' S <- cor(Y)
#'
#' # fit model
#' fit <- GGMncv(S, n = nrow(Y))
#' }
#' @export
GGMncv <- function(x, n,
                   lambda = NULL,
                   n_lambda = 100,
                   gamma = NULL,
                   select = FALSE,
                   L0_learn = FALSE,
                   penalty = "atan",
                   refit = FALSE,
                   LLA = FALSE,
                   method = "pearson",
                   progress = TRUE,
                   store = FALSE,
                   vip = FALSE,
                   vip_iter = 1000,
                   ...){

  if(! penalty %in% c("atan",
                      "mcp",
                      "scad",
                      "exp",
                      "selo",
                      "log")){
    stop("penalty not found. \ncurrent options: atan, mcp, scad, exp, selo, or log")
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

  if(is.null(lambda)){
    # tuning
    lambda_no_select <- sqrt(log(p)/n)
  } else {
      lambda_no_select <- lambda
    }


  if(is.null(gamma)) {
    if(penalty == "scad") {
      gamma <- 3.7
    } else if (penalty == "mcp") {
      gamma <- 3

    } else {
      gamma <- 0.1
      }
    }

  # high dimensional ?
  if(n > p){
    # inverse covariance matrix
    Theta <- solve(R)
  } else {
    Theta <- glassoFast::glassoFast(R,  rho = lambda_no_select)$wi
  }

  # no selection
  if(!select) {

    if (!LLA) {

      # lambda matrix
      lambda_mat <-
        eval(parse(
          text =  paste0(
            penalty,
            "_deriv(Theta = Theta, lambda =  lambda_no_select, gamma = gamma)"
          )
        ))

      diag(lambda_mat) <-  lambda_no_select

      fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

    } else {

      Theta <- glassoFast::glassoFast(S = R, rho =  lambda_no_select)$wi

      lambda_mat <-
        eval(parse(
          text =  paste0(
            penalty,
            "_deriv(Theta = Theta, lambda =  lambda_no_select, gamma = gamma)"
          )
        ))

      diag(lambda_mat) <-  lambda_no_select

      fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

    }

    fitted_models <- NULL
    } else {

    if(is.null(lambda)){
      # take from the huge package:
      # Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012).
      # The huge package for high-dimensional undirected graph estimation in R.
      # The Journal of Machine Learning Research, 13(1), 1059-1062.
      lambda.max <- max(max(R - I_p), -min(R - I_p))
      lambda.min = 0.01 * lambda.max
      lambda <- exp(seq(log(lambda.min), log(lambda.max), length.out = n_lambda))
      }

     n_lambda <- length(lambda)

     if(progress){
       message("selecting lambda")
       pb <- utils::txtProgressBar(min = 0, max = n_lambda, style = 3)
     }

     fits <- lapply(1:n_lambda, function(i){

      if (!LLA) {

        # lambda matrix
        lambda_mat <-
          eval(parse(
            text =  paste0(
              penalty,
              "_deriv(Theta = Theta, lambda = lambda[i], gamma = gamma)"
            )
          ))

        diag(lambda_mat) <- lambda[i]

        fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)
        adj <- ifelse(fit$wi == 0, 0, 1)

      } else {

        Theta <- glassoFast::glassoFast(S = R, rho = lambda[i])$wi

        lambda_mat <-
          eval(parse(
            text =  paste0(
              penalty,
              "_deriv(Theta = Theta, lambda = lambda[i], gamma = gamma)"
            )
          ))

        diag(lambda_mat) <- lambda[i]

        fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)
        adj <- ifelse(fit$wi == 0, 0, 1)
      }


      if(L0_learn){
        Theta <- hft_algorithm(Sigma = R,
                               adj = adj,
                               tol = 0.00001)$Theta

      } else {

        Theta <- fit$wi
      }

      edges <- sum(adj[upper.tri(adj)] != 0)
      log.like <- (n/2) * (log(det(Theta)) - sum(diag(R%*%Theta)))
      bic <- -2 * log.like +  edges * log(n)
      fit$bic <- bic
      fit$lambda <- lambda[i]

      if(progress){
        utils::setTxtProgressBar(pb, i)
      }
      fit
    })

   fit <-  fits[[which.min(  sapply(fits, "[[", "bic"))]]

   if(store){
     fitted_models <- fits

   } else {

     fitted_models <- NULL

   }

   }

  adj <- ifelse(fit$wi == 0, 0, 1)

  if(refit) {
    rest <- hft_algorithm(Sigma = R,
                          adj = adj,
                          tol =  0.00001)
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
      message("\ncomputing vip")
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
      lambda <- sqrt(log(p)/n)
      Theta <- solve(R)
      lambda_mat <-
        eval(parse(text =  paste0(
          penalty, "_deriv(Theta = Theta, lambda = lambda, gamma = gamma)"
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
        paste0(cn, "--", cn[x]))[upper.tri(I_p)],
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
                          vip_results = vip_results,
                          fitted_models = fitted_models)

  class(returned_object) <- c("ggmncv", "default")
  return(returned_object)
}

#' Print \code{ggmncv} Objects
#'
#' @param x An object of class \code{ggmncv}
#' @param ... Currently ignored
#' @export
print.ggmncv <- function(x,...){

  if(is(x, "default")){

  print_ggmncv(x,...)

  }
  if(is(x, "coef")){

    print_coef(x,...)
  }

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
#' @examples
#'
#' \donttest{
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # correlations
#' S <- cor(Y, method = "spearman")
#'
#' # fit model
#' fit <- GGMncv(x = S, n = nrow(Y),
#'               penalty = "atan",
#'               vip = TRUE)
#'
#' # plot
#' plot(fit, size = 4)
#' }
#'
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

print_ggmncv <- function(x, ...){
  mat <- round(x$P, 3)
  colnames(mat) <- 1:ncol(x$P)
  rownames(mat) <- 1:ncol(x$P)
  print(mat)
}
