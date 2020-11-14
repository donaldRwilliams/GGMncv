rm(list = ls())

check_gamma <- function(penalty, gamma){

  if (penalty == "scad") {
    if (any(gamma <= 3)) {
      stop("gamma must be greater than 3")
    }
  } else if (penalty == "mcp") {
    if (any(gamma <= 1)) {
      stop("gamma must be greater than 1")
    }

  } else {
    if (any(gamma < 0)) {
      stop("gamma must be positive")
    }
  }
}

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
                   ordinal = FALSE,
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
        exp(seq(log(lambda.min), log(lambda.max * 2), length.out = n_lambda))
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
      iterations =  iterations
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
    iterations =  iterations
  )

  class(returned_object) <- c("ggmncv",
                              "default")

  } else if (select == "all") {

    initial <- NULL
    lambda <- NULL
    n_lambda <- 50
    gamma <- NULL
    n_gamma <- 50
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

# lapply(1:2, function(x) fits_all[[x]])
#
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
    iterations =  iterations
  )
  } else {
    stop("select must be 'lambda', 'gamma', or 'all'")
  }
return(returned_object)
}




ptsd_pcors <- corpcor::cor2pcor(  cor(BGGM::ptsd) )
ptsd_pcors <- ifelse(abs(ptsd_pcors) < 0.05, 0, ptsd_pcors)
ptsd_cors <- corpcor::pcor2cor(ptsd_pcors)

# Y <- MASS::mvrnorm(200, rep(0,20), ptsd_cors)
Y <- BGGM::gen_ordinal(20000, 20, levels = 4, cor_mat = ptsd_cors)

R <- psych::polychoric(Y)$rho


fit1 <- ggmncv(R = cor(Y, method = "spearman"),
              n = nrow(Y),
              penalty = "atan",
              # n_gamma = 50,
              select = "lambda",
              # maxit = 2,
              # unreg = TRUE,
              store = TRUE)

test <- lapply(fits_all, unlist)

unnest <- fits_all[[1]]

system.time({
for(i in 2:n_lambda){

  unnest <- c(fits_all[[i]], unnest)

}})

t <- lapply(1:50, function(x) fits_all[[x]] )
t[[50]]

# fit2 <- ggmncv(R = cor(Y),
#                n = nrow(Y),ic = "ebic",
#                penalty = "atan",
#                select = "gamma",
#                store = TRUE,
#                unreg = FALSE, ebic_gamma = 1)

BDgraph::compare(ifelse(ptsd_pcors == 0, 0, 1), fit1$adj)




BGGM:::KL(solve(ptsd_cors), hatTheta = fit1$Theta)
BGGM:::KL(solve(ptsd_cors), hatTheta = fit2$Theta)




