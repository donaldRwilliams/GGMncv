#' Binary Classification
#'
#' @param estimate Matrix. Estimated graph (adjacency matrix)
#'
#' @param true Matrix. True graph (adjacency matrix)
#'
#' @param model_name Character string. Name of the method or penalty
#' (defaults to \code{NULL})
#'
#'
#' @return A data frame containing specificity (1 - false positive rate),
#' sensitivity (true positive rate), precision (1 - false discovery rate),
#' f1_score, and mcc (Matthews correlation coefficient).
#'
#' @export
#'
#' @examples
#' \donttest{
#' p <- 20
#' n <- 500
#'
#' true_net <- gen_net(p = p, edge_prob = 0.25)
#'
#' y <- MASS::mvrnorm(n = n,
#'                    mu = rep(0, p),
#'                    Sigma = true_net$cors)
#'
#' # default
#' fit_atan <- ggmncv(R = cor(y),
#'                    n = nrow(y),
#'                    penalty = "atan",
#'                    progress = FALSE)
#'
#' # lasso
#' fit_l1 <- ggmncv(R = cor(y),
#'                  n = nrow(y),
#'                  penalty = "lasso",
#'                  progress = FALSE)
#'
#' # atan scores
#' score_binary(estimate = true_net$adj,
#'              true = fit_atan$adj,
#'              model_name = "atan")
#'
#' score_binary(estimate = fit_l1$adj,
#'              true = true_net$adj,
#'              model_name = "lasso")
#'}
score_binary <- function (estimate, true, model_name = NULL) {

  True <- as.matrix(true)

  Estimate <- as.matrix(estimate)

  TN <-
    ifelse(True[upper.tri(True)] == 0 &
             Estimate[upper.tri(Estimate)] ==
             0, 1, 0)
  TN <- sum(TN)

  FP <-
    ifelse(True[upper.tri(True)] == 0 &
             Estimate[upper.tri(Estimate)] !=
             0, 1, 0)

  FP <- sum(FP)

  TP <-
    ifelse(True[upper.tri(True)] != 0 &
             Estimate[upper.tri(Estimate)] !=
             0, 1, 0)

  TP <- sum(TP)

  FN <-
    ifelse(True[upper.tri(True)] != 0 &
             Estimate[upper.tri(Estimate)] ==
             0, 1, 0)
  FN <- sum(FN)

  Specificity <- TN / (TN + FP)

  Sensitivity <- TP / (TP + FN)

  Precision <- TP / (TP + FP)

  Recall <- TP / (TP + FN)

  F1_score <- 2 * ((Precision * Recall) / (Precision + Recall))

  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  results <- c(Specificity,
               Sensitivity,
               Precision,
               Recall,
               F1_score,
               MCC)

  results_name <- c("specificity",
                    "sensitivity",
                    "precision",
                    "recall",
                    "f1_score",
                    "mcc")

  results <-
    cbind.data.frame(measure = results_name,
                     score = results)

  if(!is.null(model_name)){
    results$model_name <- model_name
  }
  return(results)
}
