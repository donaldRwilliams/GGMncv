#' Penalty Derivative
#'
#' @name penalty_derivative
#'
#' @description Compute the derivative for a nonconvex penalty.
#'
#' @param theta Numeric vector. Values for which the derivative is computed.
#'
#' @param penalty Character string. Which penalty should be
#'                used (defaults to \code{"atan"})?
#'                See \code{\link[GGMncv]{ggmncv}} for the
#'                available penalties.
#'
#' @param lambda Numeric.  Regularization parameter (defaults to \code{1}).
#'
#' @param gamma Numeric vector. Hyperparameter(s) for the penalty function
#'
#' @references
#' \insertAllCited{}
#'
#' @return A list of class \code{penalty_derivative}, including the following:
#'
#' \itemize{
#'
#' \item \code{deriv}: Data frame including the derivative, theta, gamma,
#' and the chosen penalty.
#'
#' \item \code{lambda}: Regularization parameter.
#'
#' }
#'
#' @note  Some care is required for specifying \code{gamma}. For example,
#' the default value for \code{scad} is 3.7 and it \emph{must} be some
#' value greater than 2 \insertCite{fan2001variable}{GGMncv}. The
#' default values in \strong{GGMncv} are set to recommended values in the
#' respective papers.
#'
#' @export
#'
#' @examples
#' deriv <- penalty_derivative(theta =  seq(-5,5,length.out = 10000),
#'                             lambda = 1,
#'                             gamma = c(0.01, 0.05, 0.1))
#'
#' head(deriv$deriv)
penalty_derivative <- function(theta = seq(-5,5, length.out = 100000),
                               penalty = "atan",
                               lambda = 1,
                               gamma = c(0.01, 0.05)){



  if(!is.matrix(theta)){

    pen_deriv <- lapply(1:length(gamma), function(x) {
      deriv_mat <-
        eval(parse(
          text =  paste0(
            penalty,
            "_deriv(Theta = as.matrix(theta), lambda = lambda, gamma = gamma[x])"
          )
        ))

      res_x <- data.frame(
        deriv = deriv_mat,
        thetas = abs(theta),
        gamma = gamma[x],
        penalty = penalty
      )

      return(res_x)

})


  deriv <- do.call(rbind.data.frame, pen_deriv)

  returned_object <- list(deriv = deriv, lambda = lambda)

  } else if (is.matrix(theta)){

    if(length(gamma) != 1){
      stop("only one value for gamma is allowed when theta is a matrix.")
    }

    deriv_mat <-
      eval(parse(
        text =  paste0(
          penalty,
          "_deriv(Theta = as.matrix(theta), lambda = lambda, gamma = gamma)"
        )
      ))

    returned_object <- list(deriv = deriv_mat)

  } else {

    stop("theta not supported. must be numeric or matrix.")

  }

  class(returned_object) <- "penalty_derivative"

  return(returned_object)
}

#' Plot \code{penalty_derivative} Objects
#'
#' @param x An object of class \code{\link{penalty_derivative}}.
#'
#' @param size Numeric. Line size in \code{geom_line}.
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#'
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

  if (is.data.frame(x$deriv)) {

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

  } else {

    stop("x not supported. must be a data frame")

  }

  return(plt)
}
