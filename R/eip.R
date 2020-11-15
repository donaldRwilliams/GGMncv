#' Bootstrapped Edge Inclusion 'Probabilities'
#'
#'
#' @description
#' \loadmathjax
#' Compute edge inclusion 'probabilities' with a non-parametric bootstrap.
#'
#' @param Y Matrix. A matrix of dimensions \emph{n} by \emph{p}.
#'
#' @param method Character string. Which correlation coefficient (or covariance) is to be
#'        computed. One of "pearson" (default), "kendall", or "spearman."
#'        Defaults to \code{pearson}.
#'
#' @param samples Numeric. How many boostrap samples (defaults to \code{500}).
#'
#' @param penalty Character string. Which penalty should be used (defaults to \code{"atan"})?
#'
#' @param ic Character string. Which information criterion should be used (defaults to \code{"bic"})?
#'           The options include \code{aic}, \code{ebic} (ebic_gamma defaults to \code{0.5}; see details),
#'           \code{ric}, or any of the generalized information criteria provided in section 5 of
#'           \insertCite{kim2012consistent;textual}{GGMncv}. The options are \code{gic_1}
#'           (i.e., \code{bic}) to \code{gic_6}.
#'
#' @param gamma Numeric. Hyperparameter for the penalty function. Defaults to 3.7 (\code{SCAD}),
#'              2 (\code{MCP}), 0.5 (\code{adapt}), and 0.01 otherwise with \code{select = "lambda"}.
#'
#' @param lambda Numeric vector. Regularization parameter. Defaults to \code{NULL} that provides default
#'               values with  \code{select = "lambda"} and  \code{sqrt(log(p)/n)} with
#'               \code{select = "gamma"}.
#'
#' @param n_lambda Numeric. The number of \mjseqn{\lambda}'s to be evaluated. Defaults to 50.
#'                 This is disregarded if custom values are provided in \code{lambda}.
#'
#' @param n_gamma Numeric. The number of \mjseqn{\gamma}'s to be evaluated. Defaults to 50.
#'                This is disregarded if custom values are provided in \code{lambda}.
#'
#' @param unreg Logical. Should the models be refitted (or unregularized) with maximum likelihood
#'              (defaults to \code{FALSE})? Setting to \code{TRUE} results in the approach of
#'              \insertCite{Foygel2010;textual}{GGMncv}, but with the regularization path obtained from
#'              nonconvex regularization, as opposed to the \mjseqn{\ell_1}-penalty.
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE}) ?
#'
#' @param ... Additional arguments. Currently gamma in EBIC (\code{ic = "ebic"}) can be set
#'            with \code{ebic_gamma = 1}.
#'
#' @return An object of class \code{eip}
#' @export
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # compute eip's
#' boot_samps <- boot_eip(Y, samples  = 10)
boot_eip <- function(Y,
                     method = "pearson",
                     samples = 50,
                     penalty = "atan",
                     ic = "bic",
                     gamma = NULL,
                     lambda = NULL,
                     n_lambda = 50,
                     n_gamma = 50,
                     unreg = FALSE,
                     progress = TRUE,...){

  n <- nrow(Y)

  p <- ncol(Y)

  I_p <- diag(p)

  if(progress){
    message("\ncomputing eip's")
    pb <- utils::txtProgressBar(min = 0, max = samples, style = 3)
  }

  boot_samps <-  sapply(1:samples, function(i) {

    Yboot <- Y[sample(1:n, size = n, replace = TRUE),]

    R <- cor(Yboot, method = method)

    fit <- ggmncv(R = R,
                      n = n,
                      penalty = penalty,
                      ic = ic,
                      select = "lambda",
                      gamma = gamma,
                      lambda = lambda,
                      n_lambda = n_lambda,
                      n_gamma = n_gamma,
                      unreg = unreg,
                      store = FALSE,
                      progress = FALSE)

    adj <- fit$adj

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

  eip_results <-
    data.frame(Relation =  sapply(1:p, function(x)
      paste0(cn, "--", cn[x]))[upper.tri(I_p)],
      EIP = rowMeans(boot_samps))

  returned_object <- list(eip_results = eip_results)
  class(returned_object) <- c("eip")

  return(returned_object)
}



#' Plot Edge Inclusion 'Probabilities'
#'
#' @param x An object of class \code{eip}
#'
#' @param color Character string. Color for \code{geom_point}.
#'
#' @param size Numeric. Size of \code{geom_point}.
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#'
#' @export
#'
#' @examples
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # compute eip's
#' boot_samps <- boot_eip(Y, B = 10)
#'
#'
#' plot(boot_samps)
plot.eip <- function(x, color = "black", size = 1,...){

  dat <- x$eip_results[order(x$eip_results$EIP),]

  dat$new1 <- factor(dat$Relation,
                     levels = dat$Relation,
                     labels = dat$Relation)

  plt <- ggplot(dat,aes(y= new1,
                        x = EIP,
                        group = new1)) +
    geom_point(size = size,
               color = color)  +
    ylab("Relation") +
    theme(axis.title  = element_text(size = 12),
          strip.text = element_text(size = 12))

  return(plt)

}

#' Print \code{eip} Objects
#'
#' @param x An object of class \code{eip}
#'
#' @param ... Currently ignored.
#'
#' @export
print.eip <- function(x, ...){
  cat("Edge Inclusion 'Probabilities':\n\n")
  print(data.frame(Relation = x$eip_results$Relation,
                   EIP = x$eip_results$EIP),
        row.names = F)
  cat("-----")
}
