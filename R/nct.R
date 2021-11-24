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
#' @param ... Additional arguments passed to \code{\link{ggmncv}}.
#'
#' @references
#' \insertAllCited{}
#'
#' @return
#'
#' @export
#'
#' @importFrom pbapply pboptions
#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
nct <- function(Y_g1, Y_g2,
                iter = 1000,
                desparsify = TRUE,
                method = "pearson",
                FUN = NULL,
                cores = 1,
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

    if(is(FUN, "function")){

      obs$fun_obs <- FUN(pcor_g1, pcor_g2)

    } else if (is(FUN, "list")){

      obs$fun_obs <- lapply(FUN, function(f) f(pcor_g1, pcor_g2))

    } else {

      stop("custom not supported. must be a function or list of functions")

    }

  }

  cl <- makeCluster(cores)

  pboptions(nout = 10)

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

  stopCluster(cl)

  custom_results <- list()
  perm_list <- list()

  if(!is.null(FUN)){

    if(length(FUN) == 1){

      if(length(obs$fun_obs) == 1){
        perm_list[[1]] <- sapply(iter_results, "[[", "fun_perm", simplify = TRUE)
        names(perm_list) <- paste0(as.character(substitute(FUN)), "_perm")
        custom_results[[1]] <- mean(as.numeric(perm_list[[1]]) >= obs$fun_obs)
        names(custom_results) <- paste0(as.character(substitute(FUN)), "_pvalue")
        obs$fun_obs <- list(obs$fun_obs)
        names(obs$fun_obs) <- paste0(as.character(substitute(FUN)), "_obs")

      }



      #names(perm_list) <- paste0(as.character(substitute(FUN)), "_perm")

      #custom_results[[1]] <- mean(perm_list[[1]] >= obs$fun_obs[[i]])

    } else {

      for(i in seq_len(length(FUN))) {

        res_x <- do.call(rbind, sapply(iter_results, "[[", "fun_perm",
                                       simplify = TRUE)[i,])

        perm_list[[i]] <- res_x

        if(ncol(res_x) > 1){

          custom_results[[i]]  <- colMeans(t(apply(res_x, MARGIN = 1, function(x) {
            abs(x) >= abs(obs$fun_obs[[i]])
          })))

        } else {

          custom_results[[i]] <- mean(res_x >= obs$fun_obs[[i]])

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
                       perm_list_default,
                       perm_list,
                       obs[1:4],
                       obs$fun_obs)

  return(returned_object)
}


