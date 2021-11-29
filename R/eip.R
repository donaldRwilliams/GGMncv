#' @title Bootstrapped Edge Inclusion 'Probabilities'
#'
#' @description
#' \loadmathjax
#' Compute the number of times each edge was selected
#' when performing a non-parametric bootstrap
#' \insertCite{@see Figure 6.7, @hastie2009elements}{GGMncv}.
#'
#' @param Y A matrix of dimensions \emph{n} by \emph{p}.
#'
#' @param method Character string. Which correlation coefficient (or covariance)
#'               is to be computed. One of "pearson" (default), "kendall",
#'               or "spearman".
#'
#' @param samples Numeric. How many bootstrap samples (defaults to \code{500})?
#'
#'
#' @param progress Logical. Should a progress bar be included (defaults to \code{TRUE})?
#'
#' @param ... Additional arguments passed to \code{\link{ggmncv}}.
#'
#' @references
#' \insertAllCited{}
#'
#' @return An object of class \code{eip} that includes the "probabilities" in a
#' data frame.
#'
#' @note Although \insertCite{hastie2009elements;textual}{GGMncv} suggests
#' this approach provides probabilities, to avoid confusion with Bayesian inference,
#' these are better thought of as "probabilities" (or better yet proportions).
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#' # data (ptsd symptoms)
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # compute eip's
#' boot_samps <- boot_eip(Y, samples  = 100, progress = FALSE)
#'
#' boot_samps
#'}
boot_eip <- function(Y,
                     method = "pearson",
                     samples = 500,
                     progress = TRUE, ...){

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

    fit <- ggmncv(R = R, n = n, progress = FALSE, ...)

    adj <- fit$adj

    if(progress){
      utils::setTxtProgressBar(pb, i)
    }

    adj[upper.tri(adj)]

  })

  if (is.null(colnames(Y))) {
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
#' \donttest{
#'
#' # data
#' Y <- GGMncv::ptsd[,1:10]
#'
#' # compute eip's
#' boot_samps <- boot_eip(Y, B = 10, progress = FALSE)
#'
#'
#' plot(boot_samps)
#'
#' }
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

#' Print the Head of \code{eip} Objects
#'
#' @param x An object of class \code{eip}
#'
#' @param n Numeric. Number of rows to print.
#'
#' @param ... Currently ignored.
#'
#' @export
head.eip <- function(x, n = 5,...){
  print(x$eip_results[1:n,], row.names = FALSE )
}
