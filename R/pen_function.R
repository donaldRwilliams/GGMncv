#' Penalty Function
#'
#' @name penalty_function
#'
#' @description Compute the penalty function for nonconvex penalties.
#'
#' @inheritParams penalty_derivative
#'
#' @references
#' \insertAllCited{}
#'
#' @return A list of class \code{penalty_function}, including the following:
#'
#' \itemize{
#'
#' \item \code{deriv}: Data frame including the penalty function,
#' theta, gamma, and the chosen penalty.
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
#' func <- penalty_function(theta =  seq(-5,5,length.out = 10000),
#'                             lambda = 1,
#'                             gamma = c(0.01, 0.05, 0.1))
#'
#' head(func)
penalty_function <- function(theta = seq(-5,5,length.out = 100000),
                             penalty = "atan",
                             lambda = 1,
                             gamma = c(0.01, 0.05)){

  pen_func <- lapply(1:length(gamma), function(x){

    pen_mat <-
      eval(parse(
        text =  paste0(
          penalty,
          "_pen(Theta = as.matrix(theta), lambda = lambda, gamma = gamma[x])"
        )
      ))

    res_x <- data.frame(
      pen = pen_mat,
      thetas = theta,
      gamma = gamma[x],
      penalty = penalty
    )

    return(res_x)

  })

  pen <- do.call(rbind.data.frame, pen_func)

  returned_object <- list(pen = pen)

  class(returned_object) <- "penalty_function"

  return(returned_object)

}


#' Plot \code{penalty_function} Objects
#'
#' @param x An object of class\code{\link{penalty_function}}.
#'
#' @param size Numeric. Line size in \code{geom_line}.
#'
#' @param ... Currently ignored.
#'
#' @return An object of class \code{ggplot}
#'
#' @export
#'
#' @examples
#' \donttest{
#' func <- penalty_function(theta =  seq(-5,5,length.out = 10000),
#'                             lambda = 1,
#'                             gamma = c(0.01, 0.05, 0.1))
#' plot(func)
#' }
plot.penalty_function <- function(x, size = 1, ...) {
  plt <- ggplot(x$pen,
                aes(
                  x = thetas,
                  y = pen,
                  color = as.factor(gamma),
                  group = gamma
                ))  +
    geom_line(size = size) +
    ylab(expression(italic(p)[lambda][gamma] ~ "(" * theta * ")")) +
    xlab(expression(theta)) +
    scale_color_discrete(name = expression(gamma))
  return(plt)
}
