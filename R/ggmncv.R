#' GGMncv
#'
#' @description
#' \loadmathjax
#' Gaussian graphical models with nonconvex regularization. A survey of these approaches is provided
#' in \insertCite{williams2020beyond;textual}{GGMncv}.
#'
#'
#' @param R Matrix. A correlation matrix of dimensions \emph{p} by \emph{p}.
#'
#' @param n Numeric. The sample size used to compute the information criterion.
#'
#' @param penalty Character string. Which penalty should be used (defaults to \code{"atan"})?
#'
#' @param ic Character string. Which information criterion should be used (defaults to \code{"bic"})?
#'           The options include \code{aic}, \code{ebic} (ebic_gamma defaults to \code{0.5}; see details),
#'           \code{ric}, or any of the generalized information criteria provided in section 5 of
#'           \insertCite{kim2012consistent;textual}{GGMncv}. The options are \code{gic_1}
#'           (i.e., \code{bic}) to \code{gic_6}.
#'
#' @param select Character string. Which tuning parameter should be selected (defaults to \code{"lambda"})?.
#'               The options include \code{"lambda"} (the regularization parameter),
#'               \code{"gamma"} (governs the 'shape'), and \code{"both"}. See details.
#'
#' @param gamma Numeric vector. Hyperparameter for the penalty function. Defaults to 3.7 (\code{SCAD}),
#'              2 (\code{MCP}), 0.5 (\code{adapt}), and 0.01 otherwise with \code{select = "lambda"}.
#'
#' @param lambda Numeric vector. Regularization parameter. Defaults to \code{NULL} that provides default
#'               values with  \code{select = "lambda"} and \code{sqrt(log(p)/n)} with
#'               \code{select = "gamma"}.
#'
#' @param n_lambda Numeric. The number of \mjseqn{\lambda}'s to be evaluated. Defaults to 50.
#'                 This is disregarded if custom values are provided in \code{lambda}.
#'
#' @param n_gamma Numeric. The number of \mjseqn{\gamma}'s to be evaluated. Defaults to 50.
#'                This is disregarded if custom values are provided in \code{lambda}.
#'
#' @param initial Matrix. The initial inverse correlation matrix for computing the penalty
#'                derivative. Defaults to \code{NULL} which uses the inverse of \code{R}.
#'
#' @param LLA Logical. Should the local linear approximation be used (default to \code{FALSE})?
#'
#' @param unreg Logical. Should the models be refitted (or unregularized) with maximum likelihood
#'              (defaults to \code{FALSE})? Setting to \code{TRUE} results in the approach of
#'              \insertCite{Foygel2010;textual}{GGMncv}, but with the regularization path obtained from
#'              nonconvex regularization, as opposed to the \mjseqn{\ell_1}-penalty.
#'
#'
#' @param maxit Numeric. The maximum number of iterations for determining convergence of the LLA
#'              algorithm (defaults to \code{1e4}). Note this can be changed to, say,
#'              \code{2} or \code{3}, which will provide  two and three-step estimators
#'              without convergence check.
#'
#' @param thr Numeric. Threshold for determining convergence of the LLA algorithm
#'            (defaults to \code{1.0e-4}).
#'
#' @param store Logical. Should all of the fitted models be saved (defaults to \code{TRUE})?.
#'
#' @param progress  Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @param ... Additional arguments. Currently gamma in EBIC (\code{ic = "ebic"}) can be set
#'            with \code{ebic_gamma = 1}.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom glassoFast glassoFast
#'
#' @return An object of class \code{ggmncv}, including:
#'
#' \itemize{
#' \item \code{Theta} Inverse covariance matrix
#' \item \code{Sigma} Covariance matrix
#' \item \code{P} Weighted adjacency matrix
#' \item \code{adj} Adjacency matrix
#' \item \code{lambda} Tuning parameter
#' \item \code{fit} glasso fitted model (a list)
#' }
#'
#' @details Several of the penalties are (continuous) approximations to the \mjseqn{\ell_0} penalty,
#' that is, best subset selection. However, the solution does not require enumerating
#' all possible models which results in a computationally efficient solution.
#'
#' \strong{L0 Approximations}
#'
#' \itemize{
#'
#' \item Atan: \code{penalty = "atan"} \insertCite{wang2016variable}{GGMncv}. This is currently the default.
#'
#' \item Seamless \mjseqn{\ell_0}: \code{penalty = "selo"} \insertCite{dicker2013variable}{GGMncv}.
#'
#' \item Exponential: \code{penalty = "exp"}  \insertCite{wang2018variable}{GGMncv}
#'
#' \item Log: \code{penalty = "log"} \insertCite{mazumder2011sparsenet}{GGMncv}.
#'
#' \item Sica: \code{penalty = "sica"}  \insertCite{lv2009unified}{GGMncv}
#'
#' }
#'
#' \strong{Additional penalties}:
#'
#' \itemize{
#'
#' \item SCAD: \code{penalty = "scad"}  \insertCite{fan2001variable}{GGMncv}.
#'
#' \item MCP: \code{penalty = "mcp"} \insertCite{zhang2010nearly}{GGMncv}.
#'
#' \item Adaptive lasso (\code{penalty = "adapt"}): Defaults to  \mjseqn{\gamma = 0.5}
#'  \insertCite{zou2006adaptive}{GGMncv}. Note that for consistency with the
#'  other penalties, \mjseqn{\gamma \rightarrow 0} provides more penalization and
#'  \mjseqn{\gamma = 1} results in \mjseqn{\ell_1} regularization.
#'
#' \item Lasso:  \code{penalty = "lasso"}  \insertCite{tibshirani1996regression}{GGMncv}.
#'
#' }
#'
#' \strong{Gamma} (\mjseqn{\gamma}):
#'
#' The \code{gamma} argument corresponds to additional hyperparameter for each penalty.
#' The defaults are set to the recommended values from the respective papers.
#'
#' \strong{LLA}
#'
#' The local linear approximate is noncovex penalties was described in
#' \insertCite{fan2009network}{GGMncv}. This is essentially a iteratively reweighted (g)lasso.
#' Note that by default \code{LLA = FALSE}. This is due to the work
#' of \insertCite{zou2008one;textual}{GGMncv}, which suggested that, so long as the starting
#' values are good enough, then a one-step estimator is sufficient. In the case of low-dimensional data,
#' the sample based inverse covariance matrix is used to compute the penalty.
#' This is expected to work well, assuming that \mjseqn{n} is sufficiently larger than  \mjseqn{p}.
#'
#' \strong{EBIC}
#'
#' When setting \code{ic = "ebic"}  the hyperparameter that determines the additional penalty to BIC is
#' passed via the \code{...} argument. This must be specified as \code{ebic_gamma = 1}. The  default is
#' \code{0.5}.
#'
#' @importFrom stats cor cov2cor
#'
#' @export
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' S <- cor(Y)
#'
#' # fit model
#' fit <- ggmncv(S, n = nrow(Y))
#'
#' # plot
#' plot(get_graph(fit))
ggmncv <- function(R,
                   n,
                   penalty = "atan",
                   ic = "bic",
                   select = "lambda",
                   gamma = NULL,
                   lambda = NULL,
                   n_lambda = 50,
                   n_gamma = 50,
                   initial = NULL,
                   LLA = FALSE,
                   unreg = FALSE,
                   maxit = 1e4,
                   thr = 1.0e-4,
                   store = TRUE,
                   progress = TRUE,
                   ...) {

  if (!penalty %in% c("atan",
                      "mcp",
                      "scad",
                      "exp",
                      "selo",
                      "log",
                      "lasso",
                      "sica",
                      "lq",
                      "adapt")) {
    stop("penalty not found. \ncurrent options: atan, mcp, scad, exp, selo, or log")
  }
  if (select == "lambda") {

    # nodes
    p <- ncol(R)

    # identity matrix
    I_p <- diag(p)

    if (is.null(initial)) {
      Theta <- solve(R)

    } else {

      Theta <- initial

    }
    if (is.null(gamma)) {
      if (penalty == "scad") {
        gamma <- 3.7

      } else if (penalty == "mcp") {
        gamma <- 2

      } else if (penalty == "adapt") {
        gamma <- 0.5

      } else {
        gamma <- 0.01

      }
    }

    if (is.null(lambda)) {
      # take from the huge package:
      # Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012).
      # The huge package for high-dimensional undirected graph estimation in R.
      # The Journal of Machine Learning Research, 13(1), 1059-1062.
      lambda.max <- max(max(R - I_p),-min(R - I_p))
      lambda.min = 0.01 * lambda.max
      lambda <-
        exp(seq(log(lambda.min), log(lambda.max), length.out = n_lambda))
    }

    n_lambda <- length(lambda)

    if (progress) {
      message("selecting lambda")
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_lambda,
                                  style = 3)
    }

    iterations <- 0

    fits <- lapply(1:n_lambda, function(i) {
      if (!LLA) {
        # lambda matrix
        lambda_mat <-
          eval(parse(
            text =  paste0(
              penalty,
              "_deriv(Theta = Theta, lambda = lambda[i], gamma = gamma)"
            )
          ))

        fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

        Theta <- fit$wi

        adj <- ifelse(fit$wi == 0, 0, 1)

      } else {
        Theta_new <- glassoFast::glassoFast(S = R, rho = lambda[i])$wi

        convergence <- 1

        iterations <- 0

        while (convergence > thr & iterations < maxit) {

          Theta_old <- Theta_new

          lambda_mat <-
            eval(parse(
              text =  paste0(
                penalty,
                "_deriv(Theta = Theta_old, lambda = lambda[i], gamma = gamma)"
              )
            ))

          fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

          Theta_new <- fit$wi

          Theta <- fit$wi

          iterations <- iterations + 1

          convergence <- mean(abs(Theta_new -  Theta_old))

          fit$iterations <- iterations

        }

        adj <- ifelse(fit$wi == 0, 0, 1)

      }

      if(unreg){

        refit <-  constrained(R, adj)
        fit$wi <- refit$Theta
        fit$w  <- refit$Sigma
        Theta  <- refit$Theta
      }

      edges <- sum(adj[upper.tri(adj)] != 0)

      fit$ic <- gic_helper(
        Theta = Theta,
        R = R,
        edges = edges,
        n = n,
        p = p,
        type = ic, ...
      )

      fit$lambda <- lambda[i]
      fit$gamma <- gamma

      if (progress) {
        utils::setTxtProgressBar(pb, i)
      }

      fit

    })

    if (store) {
      fitted_models <- fits

    } else {
      fitted_models <- NULL
    }

    which_min <- which.min(sapply(fits, "[[", "ic"))

    fit <-  fits[[which_min]]

    Theta <- fit$wi

    Sigma <- fit$w

    adj <- ifelse(fit$wi == 0, 0, 1)

    P <- -(stats::cov2cor(Theta) - I_p)

    returned_object <- list(
      Theta = Theta,
      Sigma = Sigma,
      P = P,
      fit = fit,
      adj = adj,
      lambda = lambda,
      lambda_min = lambda[which_min],
      fitted_models = fitted_models,
      penalty = penalty,
      n = n,
      iterations =  iterations,
      select = select,
      R = R
    )

    class(returned_object) <- c("ggmncv",
                                "default")

    }  else if(select == "gamma") {

    R <- cor(Y)

    initial <- NULL

    if (is.null(initial)) {

      Theta <- solve(R)

    } else {

      Theta <- initial

    }


    # nodes
    p <- ncol(R)

    # identity matrix
    I_p <- diag(p)

    if(is.null(gamma)) {

      if (penalty == "scad") {

        gamma <- seq(2.001, 5, n_gamma)

      } else if (penalty == "mcp") {

        gamma <- seq(1.001, 4, n_gamma)

      } else if (penalty == "adapt") {

        gamma <- seq(0.1, 1, n_gamma)

      } else {

        gamma <- seq(0.001, 0.1, length.out =  n_gamma)

      }
    }


    check_error <- check_gamma(penalty, gamma)

    n_gamma <- length(gamma)

    if (progress) {
      message("selecting gamma")
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_gamma,
                                  style = 3)
    }

    iterations <- 0

    lambda <- sqrt(log(p) / n)

    fits <- lapply(1:n_gamma, function(i) {

      if (!LLA) {
        # lambda matrix
        lambda_mat <-
          eval(parse(
            text =  paste0(
              penalty,
              "_deriv(Theta = Theta, lambda = lambda, gamma = gamma[i])"
            )
          ))

        fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

        Theta <- fit$wi

        adj <- ifelse(fit$wi == 0, 0, 1)

      } else {

        Theta_new <- glassoFast::glassoFast(S = R, rho = lambda)$wi

        convergence <- 1

        iterations <- 0

        while (convergence > thr & iterations < maxit) {

          Theta_old <- Theta_new

          lambda_mat <-
            eval(parse(
              text =  paste0(
                penalty,
                "_deriv(Theta = Theta_old, lambda = lambda, gamma = gamma[i])"
              )
            ))

          fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

          Theta_new <- fit$wi

          Theta <- fit$wi

          iterations <- iterations + 1

          convergence <- mean(abs(Theta_new -  Theta_old))

          fit$iterations <- iterations

        }

        adj <- ifelse(fit$wi == 0, 0, 1)

      }

      if(unreg){

        refit <-  constrained(R, adj)
        fit$wi <- refit$Theta
        fit$w  <- refit$Sigma
      }

      edges <- sum(adj[upper.tri(adj)] != 0)

      fit$ic <- gic_helper(
        Theta = Theta,
        R = R,
        edges = edges,
        n = n,
        p = p,
        type = ic, ...
      )

      fit$gamma <- gamma[i]
      fit$lambda <- lambda

      if (progress) {
        utils::setTxtProgressBar(pb, i)
      }

      fit

    })

    if (store) {
      fitted_models <- fits

    } else {
      fitted_models <- NULL
    }

    which_min <- which.min(sapply(fits, "[[", "ic"))

    fit <-  fits[[which_min]]

    Theta <- fit$wi

    Sigma <- fit$w

    adj <- ifelse(fit$wi == 0, 0, 1)

    P <- -(stats::cov2cor(Theta) - I_p)

    returned_object <- list(
      Theta = Theta,
      Sigma = Sigma,
      P = P,
      fit = fit,
      adj = adj,
      lambda = lambda,
      gamma_min = gamma[which_min],
      fitted_models = fitted_models,
      penalty = penalty,
      n = n,
      iterations =  iterations,
      select = select,
      R = R
    )

    class(returned_object) <- c("ggmncv",
                                "default")

  } else if (select == "both") {

    # nodes
    p <- ncol(R)

    # identity matrix
    I_p <- diag(p)



    if (is.null(initial)) {
      Theta <- solve(R)

    } else {

      Theta <- initial

    }

    if (is.null(lambda)) {
      # take from the huge package:
      # Zhao, T., Liu, H., Roeder, K., Lafferty, J., & Wasserman, L. (2012).
      # The huge package for high-dimensional undirected graph estimation in R.
      # The Journal of Machine Learning Research, 13(1), 1059-1062.
      lambda.max <- max(max(R - I_p),-min(R - I_p))
      lambda.min = 0.01 * lambda.max
      lambda <-
        exp(seq(log(lambda.min), log(lambda.max), length.out = n_lambda))
    }

    n_lambda <- length(lambda)

    if(is.null(gamma)) {

      if (penalty == "scad") {

        gamma <- seq(2.001, 5, n_gamma)

      } else if (penalty == "mcp") {

        gamma <- seq(1.001, 4, n_gamma)

      } else if (penalty == "adapt") {

        gamma <- seq(0.1, 1, n_gamma)

      } else {

        gamma <- seq(0.001, 0.1, length.out =  n_gamma)

      }
    }

    if (progress) {
      message("selecting lambda and gamma")
      pb <- utils::txtProgressBar(min = 0,
                                  max = n_lambda,
                                  style = 3)
    }
    iterations <- 0

    fits_all <-  lapply(1:n_lambda, function(x) {


      fits <- lapply(1:n_gamma, function(i) {

        if (!LLA) {
          # lambda matrix
          lambda_mat <-
            eval(parse(
              text =  paste0(
                penalty,
                "_deriv(Theta = Theta, lambda = lambda[x], gamma = gamma[i])"
              )
            ))

          fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

          Theta <- fit$wi

          adj <- ifelse(fit$wi == 0, 0, 1)

        } else {

          Theta_new <- glassoFast::glassoFast(S = R, rho = lambda[x])$wi

          convergence <- 1

          iterations <- 0

          while (convergence > thr & iterations < maxit) {
            Theta_old <- Theta_new

            lambda_mat <-
              eval(parse(
                text =  paste0(
                  penalty,
                  "_deriv(Theta = Theta_old, lambda = lambda[x], gamma = gamma[i])"
                )
              ))

            fit <- glassoFast::glassoFast(S = R, rho = lambda_mat)

            Theta_new <- fit$wi

            Theta <- fit$wi

            iterations <- iterations + 1

            convergence <- mean(abs(Theta_new -  Theta_old))

            fit$iterations <- iterations

          }

          adj <- ifelse(fit$wi == 0, 0, 1)

        }


        if(isTRUE(unreg)){

          refit <-  constrained(Sigma = R, adj = adj)
          fit$wi <- refit$Theta
          fit$w  <- refit$Sigma
          Theta <- refit$Theta
        }

        edges <- sum(adj[upper.tri(adj)] != 0)

        fit$ic <- gic_helper(
          Theta = Theta,
          R = R,
          edges = edges,
          n = n,
          p = p,
          type = ic
        )

        fit$gamma <- gamma[i]
        fit$lambda <- lambda[x]


        fit

      })

      if (progress) {
        utils::setTxtProgressBar(pb, x)
      }

      fits
    })


    unnest <- fits_all[[1]]
    for(i in 2:n_lambda){
      unnest <- c(fits_all[[i]], unnest)
    }

    if (store) {
      fitted_models <- unnest

    } else {
      fitted_models <- NULL
    }

    which_min <- which.min(sapply(unnest, "[[" , "ic"))

    fit <-  unnest[[which_min]]

    Theta <- fit$wi

    Sigma <- fit$w

    adj <- ifelse(fit$wi == 0, 0, 1)

    P <- -(stats::cov2cor(Theta) - I_p)

    returned_object <- list(
      Theta = Theta,
      Sigma = Sigma,
      P = P,
      fit = fit,
      adj = adj,
      lambda = lambda,
      gamma_min = fit$gamma,
      lambda_min = fit$lambda,
      fitted_models = fitted_models,
      penalty = penalty,
      n = n,
      iterations =  iterations,
      select = select,
      R = R
    )
  } else {
    stop("select must be 'lambda', 'gamma', or 'both'.")
  }
  return(returned_object)
}


#' Print \code{ggmncv} Objects
#'
#' @param x An object of class \code{ggmncv}
#'
#' @param ... Currently ignored
#'
#' @importFrom methods is
#'
#' @export
print.ggmncv <- function(x,...){

  if(methods::is(x, "default")){

    print_ggmncv(x,...)

  }
  if(methods::is(x, "coef")){

    print_coef(x,...)
  }

  if(methods::is(x, "inference")){
    print_inference(x, ...)
  }

  if(methods::is(x, "ggm_compare")){
    print_compare(x,...)
  }

}



#' Plot \code{ggmncv} Objects
#'
#' @description Plot the solution path for the partial correlations.
#'
#' @param x An object of class \code{ggmncv}
#'
#' @param size Numeric. The size of the points (\code{eip}).
#'
#' @param alpha Numeric. The transparency of the lines. Only for the solution path options.
#'
#' @param ... Currently ignored.
#'
#' @return A \code{ggplot} object
#'
#' @importFrom  ggplot2 aes ggplot  geom_point ylab facet_grid geom_line
#' geom_vline geom_hline xlab ylab ggtitle theme element_text
#'
#' @importFrom reshape melt
#'
#' @examples
#'
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # correlations
#' S <- cor(Y, method = "spearman")
#'
#' # fit model
#' fit <- ggmncv(R = S, n = nrow(Y))
#'
#' # plot
#' plot(fit)
#'
#' @export
plot.ggmncv <- function(x, size = 1, alpha = 0.5, ...){

  if(x$select != "lambda"){
    stop("select must be 'lambda'.")
  }

  if(is.null(x$fitted_models)){
    stop("models not stored.")
  }

  n_lambda <-  length(x$lambda)

  if(n_lambda == 0) {
    stop("solution path not found. must set 'select = TRUE'.")
  }

  p <- ncol(x$Theta)
  which_min <- which(x$lambda  == x$lambda_min)
  lambda_min <- x$lambda[which_min]

  Theta_std <- t(sapply(1:n_lambda, function(i)
    cov2cor(x$fitted_models[[i]]$wi)[upper.tri(diag(p))]))

  non_zero <- sum(Theta_std[which_min,] !=0)
  dat_res <-  reshape::melt(Theta_std)
  dat_res$X1 <- round(x$lambda, 3)
  dat_res$penalty <- x$penalty

  plt <- ggplot(dat_res, aes(y = -value,
                             x = X1,
                             group = as.factor(X2),
                             color = as.factor(X2))) +
    facet_grid(~ penalty) +

    geom_line(show.legend = FALSE,
              alpha = alpha,
              size = size) +

    geom_hline(yintercept = 0,  color = "black") +
    xlab(expression(lambda)) +
    ylab(expression(hat(rho))) +
    ggtitle(paste0(non_zero, " edges (",
                   round(non_zero /  p*(p-1)*.5),
                   "% connectivity)")) +
    theme(axis.title  = element_text(size = 12),
          strip.text = element_text(size = 12))

  return(plt)
}



print_ggmncv <- function(x, ...){
  mat <- round(x$P, 3)
  colnames(mat) <- 1:ncol(x$P)
  rownames(mat) <- 1:ncol(x$P)
  print(mat)
}
