#' @title Network Comparison Test
#'
#' @description A re-implementation and extension of the permutation based
#' network comparison test introduced in \insertCite{van2017comparing;textual}{GGMncv}.
#' Such extensions include scaling to networks with many nodes and the option to
#' use custom test-statistics.
#'
#' @param Y_g1 A matrix (or data.frame) of dimensions \emph{n} by \emph{p},
#'             corresponding to the first dataset (\emph{p} must be the same
#'             for \code{Y_g1} and \code{Y_g2}).
#'
#' @param Y_g2 A matrix of dimensions \emph{n} by \emph{p}, corresponding to the
#'             second dataset (\emph{p} must be the same for \code{Y_g1} and \code{Y_g2}).
#'
#' @param iter Numeric. Number of (Monte Carlo) permutations (defaults to \code{1000}).
#'
#' @param desparsify Logical. Should the de-sparsified glasso estimator be
#'                   computed (defaults to \code{TRUE})? This is much faster,
#'                   as the tuning parameter is fixed to
#'                   \mjseqn{\lambda = \sqrt{log(p)/n}}.
#'
#' @param method character string. Which correlation coefficient (or covariance)
#'               is to be computed. One of "pearson" (default), "kendall",
#'               or "spearman".
#'
#' @param FUN A function or list of functions (defaults to \code{NULL}),
#'            specifying custom test-statistics. See \strong{Examples}.
#'
#'
#' @param cores Numeric. Number of cores to use when executing the permutations in
#'              parallel (defaults to \code{1}).
#'
#' @param progress Logical. Should a progress bar be included
#'                 (defaults to \code{TRUE})?
#'
#' @param update_progress How many times should the progress bar be updated
#'                        (defaults to \code{4})? Note that setting this to a
#'                        large value should result in the worse performance,
#'                        due to additional overhead communicating among the
#'                        parallel processes.
#'
#' @param ... Additional arguments passed to \code{\link{ggmncv}}.

#'
#' @details
#'
#' \strong{User-Defined Functions}
#'
#' These functions must have two arguments, corresponding
#' to the partial correlation network for each group. An
#' example is provided below.
#'
#'
#' For user-defined functions (\code{FUN}), absolute values are used
#' to compute the p-value, assuming more than one value is returned
#' (e.g., centrality). This is done to mimic the \code{R} package
#' \strong{NCT}.
#'
#' A fail-safe method to ensure the p-value is computed correctly is
#' to access the permutations and observed values from the \code{nct}
#' object.
#'
#' Finally, comparing edges is not implemented. The most straightforward
#' way to do this is with \code{\link{compare_edges}}, which
#' uses the de-sparsified estimator.
#'
#' @note
#'
#' In \insertCite{van2017comparing;textual}{GGMncv}, it was suggested that
#' these are tests of \emph{invariance}. To avoid confusion, that
#' terminology is not used in \strong{GGMncv}. This is because
#' these tests assume invariance or the null is \emph{true}, and thus
#' can only be used to detect differences. Hence, it would be incorrect
#' to suggest networks are the same, or evidence for invariance,
#' by merely failing to reject the null hypothesis
#' \insertCite{williams_null}{GGMncv}.
#'
#' For the defaults, Jensen-Shannon divergence is a symmetrized version
#' of Kullback-Leibler divergence (the average of both directions).
#'
#' \strong{Computational Speed}
#'
#' This implementation has two key features that should make it
#' scale to larger networks: (1) parallel computation and (2) the
#' \code{R} package \strong{glassoFast} is used under the hood
#' (as opposed to \strong{glasso}). CPU (time) comparisons are
#' provided in \insertCite{sustik2012glassofast;textual}{GGMncv}.
#'
#' \strong{Non-regularized}
#'
#' Non-regularized can be implemented by setting \code{lambda = 0}. Note
#' this is provided to \code{\link{ggmncv}} via \code{...}.
#'
#' @references
#' \insertAllCited{}
#'
#' @return A list of class \code{nct}, including the following
#'
#' \itemize{
#'
#' \item \code{glstr_pvalue}: Global strength p-value.
#'
#' \item \code{sse_pvalue}: Sum of square error p-value.
#'
#' \item \code{jsd_pvalue}: Jensen-Shannon divergence p-value.
#'
#' \item \code{max_pvalue}: Maximum difference p-value.
#'
#'  \item \code{glstr_obs}: Global strength observed.
#'
#' \item \code{sse_obs}: Sum of square error observed.
#'
#' \item \code{jsd_obs}: Jensen-Shannon divergence observed.
#'
#' \item \code{max_obs}: Maximum difference observed.
#'
#' \item \code{glstr_perm}: Global strength permutations.
#'
#' \item \code{sse_perm}: Sum of square error permutations.
#'
#' \item \code{jsd_perm}: Jensen-Shannon divergence permutations.
#'
#' \item \code{max_perm}: Maximum difference permutations.
#'
#'
#' }
#'
#' For user-defined functions, i.e., those provided to \code{FUN},
#' the function name is pasted to \code{_pvalue}, \code{_obs}, and
#' \code{_perm}.
#'
#' @export
#'
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapply
#'
#' @examples
#' \donttest{
#' # generate network
#' main <- gen_net(p = 10)
#'
#' # assume groups are equal
#' y1 <- MASS::mvrnorm(n = 500,
#'                     mu = rep(0, 10),
#'                     Sigma = main$cors)
#'
#' y2 <- MASS::mvrnorm(n = 500,
#'                     mu = rep(0, 10),
#'                     Sigma = main$cors)
#'
#' compare_ggms <- nct(y1, y2, iter = 500,
#'                     progress = FALSE)
#'
#' compare_ggms
#'
#' # custom function
#' # note: x & y are partial correlation networks
#'
#' # correlation
#' Correlation <- function(x, y){
#' cor(x[upper.tri(x)], y[upper.tri(y)])
#' }
#'
#' compare_ggms <- nct(y1, y2,iter = 100,
#'                     FUN = Correlation,
#'                     progress = FALSE)
#'
#' compare_ggms
#'
#' # correlation and strength
#'
#' Strength <- function(x, y){
#' NetworkToolbox::strength(x) - NetworkToolbox::strength(y)
#' }
#'
#' compare_ggms <- nct(y1, y2, iter = 100,
#'                     FUN = list(Correlation = Correlation,
#'                                Strength = Strength),
#'                     progress = FALSE)
#'
#' compare_ggms
#' }
nct <- function(Y_g1, Y_g2,
                iter = 1000,
                desparsify = TRUE,
                method = "pearson",
                FUN = NULL,
                cores = 1,
                progress = TRUE,
                update_progress = 4,
                ...){

  p_g1 <- ncol(Y_g1)

  p_g2 <- ncol(Y_g2)

  if(p_g1 != p_g2){
    stop("Yg1 and Yg2 must have the same number of nodes (columns)")
  }

  I_p <- diag(p_g1)

  n_g1 <- nrow(Y_g1)

  n_g2 <- nrow(Y_g2)

  n_total <- n_g1 + n_g2

  if(desparsify){

    fit_g1 <-
      desparsify(
        ggmncv(
          R = cor(Y_g1, method = method),
          n = n_g1, lambda = sqrt(log(p_g1)/n_g1),
          progress = FALSE
        )
      )

    fit_g2 <-
      desparsify(
        ggmncv(
          R = cor(Y_g2, method = method),
          n = n_g2, lambda = sqrt(log(p_g2)/n_g2),
          progress = FALSE
        )
      )

  } else {

    fit_g1 <- ggmncv(R = cor(Y_g1, method = method),
                     n = n_g1, ...,
                     progress = FALSE)

    fit_g2 <- ggmncv(R = cor(Y_g2, method = method),
                     n = n_g2, ...,
                     progress = FALSE)
  }

  pcor_g1 <- fit_g1$P
  pcor_g2 <- fit_g2$P

  # observed: global strength
  glstr_obs_g1 <- sum(abs(pcor_g1[upper.tri(I_p)]))
  glstr_obs_g2 <- sum(abs(pcor_g2[upper.tri(I_p)]))
  glstr_obs <- abs(glstr_obs_g1 - glstr_obs_g2)

  # observed: SSE
  sse_obs <- sum((pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)])^2)

  # observed: jensen shannon distance
  jsd_obs <- jsd(fit_g1$Theta, fit_g2$Theta)

  # observed: max diff
  max_obs <- max(abs(pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)]))

  obs <- list(glstr_obs = glstr_obs,
              sse_obs = sse_obs,
              jsd_obs = jsd_obs,
              max_obs = max_obs)

  n_seq <- seq_len(n_total)

  stacked_data <- rbind(Y_g1, Y_g2)

  if(!is.null(FUN)){

    if(length(FUN) == 1 & is(FUN, "list")){
      stop("one function in `FUN` cannot be a list.")
    }

    if(is(FUN, "function")){

      obs$fun_obs <- list(FUN(pcor_g1, pcor_g2))

    } else if (is(FUN, "list")){

      obs$fun_obs <- lapply(FUN, function(f) f(pcor_g1, pcor_g2))

    } else {

      stop("custom not supported. must be a function or list of functions")

    }

  }

  if(cores == 1){

    if(progress){

      pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)

    }

    iter_results <- lapply(X = 1:iter, function(x){

      perm_g1 <- sample(1:n_total, size = n_g1, replace = FALSE)

      Y_g1_perm <- stacked_data[perm_g1,]

      Y_g2_perm <- stacked_data[n_seq[-perm_g1],]

      if(desparsify){

        fit_g1 <-
          desparsify(
            ggmncv(
              R = cor(Y_g1_perm, method = method),
              n = n_g1, lambda = sqrt(log(p_g1)/n_g1),
              progress = FALSE
            )
          )

        fit_g2 <-
          desparsify(
            ggmncv(
              R = cor(Y_g2_perm, method = method),
              n = n_g2, lambda = sqrt(log(p_g2)/n_g2),
              progress = FALSE
            )
          )

      } else {

        fit_g1 <- ggmncv(R = cor(Y_g1_perm, method = method),
                         n = n_g1, ...,
                         progress = FALSE)

        fit_g2 <- ggmncv(R = cor(Y_g2_perm, method = method),
                         n = n_g2, ...,
                         progress = FALSE)
      }

      pcor_g1 <- fit_g1$P
      pcor_g2 <- fit_g2$P

      # observed: global strength
      glstr_perm_g1 <- sum(abs(pcor_g1[upper.tri(I_p)]))
      glstr_perm_g2 <- sum(abs(pcor_g2[upper.tri(I_p)]))
      glstr_diff_perm <- abs(glstr_perm_g1 - glstr_perm_g2)

      # observed: SSE
      sse_perm <- sum((pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)])^2)

      # observed: jensen shannon distance
      jsd_perm <- jsd(fit_g1$Theta, fit_g2$Theta)

      # observed: max diff
      max_diff_perm <- max(abs(pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)]))

      returned_obj <- list(glstr_perm = glstr_diff_perm,
                           sse_perm = sse_perm,
                           jsd_perm = jsd_perm,
                           max_perm = max_diff_perm)


      if(!is.null(FUN)){

        if(is(FUN, "function")){

          returned_obj$fun_perm <- FUN(pcor_g1, pcor_g2)

        } else if (is(FUN, "list")){

          returned_obj$fun_perm <- lapply(FUN, function(f) f(pcor_g1, pcor_g2))

        } else {

          stop("custom not supported. must be a function or list of functions")

        }

      }

      if(progress){

        utils::setTxtProgressBar(pb, x)

      }

      return(returned_obj)

    })

  } else {

    cl <- parallel::makeCluster(cores)

    if(progress){

      pbapply::pboptions(nout = update_progress)

      iter_results <- pbapply::pblapply(X = 1:iter, function(x){

        perm_g1 <- sample(1:n_total, size = n_g1, replace = FALSE)

        Y_g1_perm <- stacked_data[perm_g1,]

        Y_g2_perm <- stacked_data[n_seq[-perm_g1],]

        if(desparsify){

          fit_g1 <-
            desparsify(
              ggmncv(
                R = cor(Y_g1_perm, method = method),
                n = n_g1, lambda = sqrt(log(p_g1)/n_g1),
                progress = FALSE
              )
            )

          fit_g2 <-
            desparsify(
              ggmncv(
                R = cor(Y_g2_perm, method = method),
                n = n_g2, lambda = sqrt(log(p_g2)/n_g2),
                progress = FALSE
              )
            )

        } else {

          fit_g1 <- ggmncv(R = cor(Y_g1_perm, method = method),
                           n = n_g1, ...,
                           progress = FALSE)

          fit_g2 <- ggmncv(R = cor(Y_g2_perm, method = method),
                           n = n_g2, ...,
                           progress = FALSE)
        }

        pcor_g1 <- fit_g1$P
        pcor_g2 <- fit_g2$P

        # observed: global strength
        glstr_perm_g1 <- sum(abs(pcor_g1[upper.tri(I_p)]))
        glstr_perm_g2 <- sum(abs(pcor_g2[upper.tri(I_p)]))
        glstr_diff_perm <- abs(glstr_perm_g1 - glstr_perm_g2)

        # observed: SSE
        sse_perm <- sum((pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)])^2)

        # observed: jensen shannon distance
        jsd_perm <- jsd(fit_g1$Theta, fit_g2$Theta)

        # observed: max diff
        max_diff_perm <- max(abs(pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)]))

        returned_obj <- list(glstr_perm = glstr_diff_perm,
                             sse_perm = sse_perm,
                             jsd_perm = jsd_perm,
                             max_perm = max_diff_perm)


        if(!is.null(FUN)){

          if(is(FUN, "function")){

            returned_obj$fun_perm <- FUN(pcor_g1, pcor_g2)

          } else if (is(FUN, "list")){

            returned_obj$fun_perm <- lapply(FUN, function(f) f(pcor_g1, pcor_g2))

          } else {

            stop("custom not supported. must be a function or list of functions")

          }

        }

        return(returned_obj)

      }, cl = cl)

    } else {

      iter_results <- parallel::parLapply(X = 1:iter, function(x){

        perm_g1 <- sample(1:n_total, size = n_g1, replace = FALSE)

        Y_g1_perm <- stacked_data[perm_g1,]

        Y_g2_perm <- stacked_data[n_seq[-perm_g1],]

        if(desparsify){

          fit_g1 <-
            desparsify(
              ggmncv(
                R = cor(Y_g1_perm, method = method),
                n = n_g1, lambda = sqrt(log(p_g1)/n_g1),
                progress = FALSE
              )
            )

          fit_g2 <-
            desparsify(
              ggmncv(
                R = cor(Y_g2_perm, method = method),
                n = n_g2, lambda = sqrt(log(p_g2)/n_g2),
                progress = FALSE
              )
            )

        } else {

          fit_g1 <- ggmncv(R = cor(Y_g1_perm, method = method),
                           n = n_g1, ...,
                           progress = FALSE)

          fit_g2 <- ggmncv(R = cor(Y_g2_perm, method = method),
                           n = n_g2, ...,
                           progress = FALSE)
        }

        pcor_g1 <- fit_g1$P
        pcor_g2 <- fit_g2$P

        # observed: global strength
        glstr_perm_g1 <- sum(abs(pcor_g1[upper.tri(I_p)]))
        glstr_perm_g2 <- sum(abs(pcor_g2[upper.tri(I_p)]))
        glstr_diff_perm <- abs(glstr_perm_g1 - glstr_perm_g2)

        # observed: SSE
        sse_perm <- sum((pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)])^2)

        # observed: jensen shannon distance
        jsd_perm <- jsd(fit_g1$Theta, fit_g2$Theta)

        # observed: max diff
        max_diff_perm <- max(abs(pcor_g1[upper.tri(I_p)] - pcor_g2[upper.tri(I_p)]))

        returned_obj <- list(glstr_perm = glstr_diff_perm,
                             sse_perm = sse_perm,
                             jsd_perm = jsd_perm,
                             max_perm = max_diff_perm)

        if(!is.null(FUN)){

          if(is(FUN, "function")){

            returned_obj$fun_perm <- FUN(pcor_g1, pcor_g2)

          } else if (is(FUN, "list")){

            returned_obj$fun_perm <- lapply(FUN, function(f) f(pcor_g1, pcor_g2))

          } else {

            stop("custom not supported. must be a function or list of functions")

          }

        }

        return(returned_obj)

      }, cl = cl)

    }

    stopCluster(cl)

  }

  custom_results <- list()

  perm_list <- list()

  if(!is.null(FUN)){

    if (length(FUN) == 1) {

      if (length(obs$fun_obs[[1]]) == 1) {

        perm_list[[1]] <-
          sapply(iter_results, "[[", "fun_perm", simplify = TRUE)

        custom_results[[1]] <-
          mean(as.numeric(perm_list[[1]]) >= obs$fun_obs[[1]])

        names(perm_list) <-
          paste0(as.character(substitute(FUN)), "_perm")

        names(custom_results) <-
          paste0(as.character(substitute(FUN)), "_pvalue")

        names(obs$fun_obs) <-
          paste0(as.character(substitute(FUN)), "_obs")

      } else {

        perm_list[[1]] <-  t(sapply(iter_results, "[[", "fun_perm",
                                  simplify = TRUE))

        custom_results[[1]]  <- colMeans(apply(perm_list[[1]],
                                               MARGIN = 2, function(x) {
                    abs(x) >= abs(obs$fun_obs[[1]])
                  }))

        names(custom_results) <-
          paste0(as.character(substitute(FUN)), "_pvalue")

        names(perm_list) <-
          paste0(as.character(substitute(FUN)), "_perm")

        names(obs$fun_obs) <-
          paste0(as.character(substitute(FUN)), "_obs")

      } # end more than 1 test stat

    } else {

      for(i in seq_along(FUN)){

        res_i <-
          do.call(rbind,
                  t(sapply(iter_results, "[[", "fun_perm",
                           simplify = TRUE))[,i])

        perm_list[[i]] <- res_i

        if (length(obs$fun_obs[[i]]) > 1) {

          custom_results[[i]]  <- colMeans(apply(res_i,
                                                 MARGIN = 2, function(x) {
                                                   abs(x) >= abs(obs$fun_obs[[i]])
                                                 }))

        } else {

          custom_results[[i]] <- mean(res_i >= obs$fun_obs[[i]])

        }

      }

      names(custom_results) <- paste0(names(FUN), "_pvalue")

      names(perm_list) <- paste(names(FUN), "_perm")

      names(obs$fun_obs) <- paste0(names(obs$fun_obs), "_obs")

    }

  }

  default_names <- names(iter_results[[1]])[1:4]

  perm_list_default <- lapply(default_names, function(x) {
    results_x <- sapply(iter_results, "[[", x)
    return(results_x)
  })

  default_results <- lapply(1:4, function(x) {
    results_x <- mean(perm_list_default[[x]] >= obs[[x]])
    return(results_x)
  })

  names(default_results) <- paste0(sub(pattern = "\\_.*",
                                       replacement = "",
                                       x = default_names), "_pvalue")

  names(perm_list_default) <- default_names

  returned_object <- c(default_results,
                       custom_results,
                       obs[1:4],
                       obs$fun_obs,
                       perm_list_default,
                       perm_list)

  class(returned_object) <- "nct"

  return(returned_object)
}

#' Print \code{nct} Objects
#'
#' @param x An object of class \code{nct}
#' @param ... Currently ignored.
#' @export
print.nct <- function(x, ...){

  check_defaults <- length(grep("_pvalue", names(x)))
  cat("Network Comparsion Test\n")
  cat("(GGMncv Edition)\n")
  cat("----\n")
  cat("Maximum Difference\n")
  cat("p-value:", x$max_pvalue, "\n")
  cat("----\n")
  cat("Global Strength\n")
  cat("p-value:", x$glstr_pvalue, "\n")
  cat("----\n")
  cat("Sum of Squared Error\n")
  cat("p-value:", x$sse_pvalue, "\n")
  cat("----\n")
  cat("Jensen-Shannon divergence\n")
  cat("p-value:", x$jsd_pvalue, "\n")
  cat("----\n")

  if (check_defaults > 4) {
    cat("note: compute p-values manually for custom tests. see vignettes.")
  }

}


