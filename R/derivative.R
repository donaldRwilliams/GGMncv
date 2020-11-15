#' Penalty Derivative
#'
#' @description Compute the derivative for a nonconvex penalty.
#'
#' @param theta Numeric vector. Values for which the derivative is computed.
#'
#' @param penalty Character string. Character string. Which penalty should be
#'                used (defaults to \code{"atan"})? See \code{\link[GGMncv]{ggmncv}} for the
#'                available penalties.
#'
#' @param lambda Numeric.  Regularization parameter (defaults to \code{1}).
#'
#' @param gamma Numeric vector. Hyperparameter(s) for the penalty function.
#'
#' @return A list of class \code{penalty_derivative}
#' @export
#'
#' @examples
#' deriv <- penalty_derivative(theta =  seq(-5,5,length.out = 10000),
#'                             lambda = 1,
#'                             gamma = c(0.01, 0.05, 0.1))
penalty_derivative <- function(theta = seq(-5,5,length.out = 100000),
                               penalty = "atan",
                               lambda = 1,
                               gamma = c(0.01, 0.05)){

  pen_deriv <- lapply(1:length(gamma), function(x){
    deriv_mat <-
      eval(parse(
        text =  paste0(
          penalty,
          "_deriv(Theta = as.matrix(theta), lambda = lambda, gamma = gamma[x])"
        )
      ))
    data.frame(deriv = deriv_mat,
               thetas = abs(theta),
               gamma = gamma[x],
               penalty = penalty)
  })
  deriv <- do.call(rbind.data.frame, pen_deriv)
  returned_object <- list(deriv = deriv, lambda = lambda)
  class(returned_object) <- "penalty_derivative"
  return(returned_object)
}



#' Plot \code{penalty_derivative} Objects
#'
#' @param x An object of class \code{penalty_derivative}.
#'
#' @param size Numeric. Line size in \code{geom_line}.
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#' @export
#'
#' @importFrom ggplot2 scale_color_discrete scale_y_continuous
#'
#' @examples
#' \donttest{
#' pen_deriv <- penalty_derivative(theta =  seq(-5,5,length.out = 10000),
#'                             lambda = 1,
#'                             gamma = c(0.01, 0.05, 0.1))
#' plot(pen_deriv)
#' }
plot.penalty_derivative <- function(x, size = 1, ...) {

  plt <- ggplot(x$deriv,
                aes(
                  x = thetas,
                  y = deriv,
                  color = as.factor(gamma),
                  group = gamma
                ))  +
    geom_line(size = size) +
    ylab(expression(italic(p * "'")[lambda][gamma] ~ "(" * theta * ")")) +
    xlab(expression(theta)) +
    scale_color_discrete(name = expression(gamma)) +
    scale_y_continuous(limits = c(0, x$lambda))

  return(plt)
}
