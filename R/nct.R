#' Title
#'
#' @param Y_g1
#' @param Y_g2
#' @param iter
#' @param desparsify
#' @param method
#' @param FUN
#' @param cores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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

  cl <- parallel::makeCluster(cores)

  pbapply::pboptions(nout = 10)

  parallel::clusterExport(cl,varlist = "jsd")

  iter_results <- pblapply(X = 1:iter, function(x){

    perm_g1 <- sample(1:n_total, size = n_g1, replace = FALSE)

    Y_g1_perm <- stacked_data[perm_g1,]

    Y_g2_perm <- stacked_data[n_seq[-perm_g1],]

    if(desparsify){

      fit_g1 <-
        GGMncv::desparsify(
          GGMncv::ggmncv(
            R = cor(Y_g1_perm, method = method),
            n = n_g1, lambda = sqrt(log(p_g1)/n_g1),
            progress = FALSE
          )
        )

      fit_g2 <-
        GGMncv::desparsify(
          GGMncv::ggmncv(
            R = cor(Y_g2_perm, method = method),
            n = n_g2, lambda = sqrt(log(p_g2)/n_g2),
            progress = FALSE
          )
        )

    } else {

      fit_g1 <- GGMncv::ggmncv(R = cor(Y_g1_perm, method = method),
                       n = n_g1, ...,
                       progress = FALSE)

      fit_g2 <- GGMncv::ggmncv(R = cor(Y_g2_perm, method = method),
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

  parallel::stopCluster(cl)

  custom_results <- list()

  if(!is.null(FUN)){

    if(length(FUN) == 1){

      custom_results[[1]] <- t(sapply(iter_results, "[[", "fun_perm", simplify = TRUE))

    } else {

      for(i in seq_len(length(FUN))) {

        res_x <- do.call(rbind, sapply(iter_results, "[[", "fun_perm", simplify = TRUE)[i,])

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

    }

    names(obs$fun_obs) <- paste0(names(obs$fun_obs), "_obs")

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

  ret <- c(default_results,
    custom_results,
    perm_list_default,
    obs[1:4],
    obs$fun_obs)

  return(ret)
}


