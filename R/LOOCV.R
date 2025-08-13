#' Calculates the LOOCV for the model fit
#' Needs more work
#'
#' @noRd

LOOCV <- function(object) {
  if (!inherits(object, "MuFiMeshGP")) {
    stop("The object is not of class \"MuFiMeshGP\" \n")
  }
  X <- object$X
  t <- object$t
  Y <- object$Y
  Ki <- object$Ki
  l <- object$l
  if (object$used_args$trend.type == "SK") {
    beta0 <- object$used_args$mean.known
    pred <- Y - (Ki %*% (Y - beta0)) / diag(Ki)
    sd2 <- 1 / diag(Ki)
    rmse <- mean((pred - Y)^2)
  }
  if (object$used_args$trend.type == "OK") {
    beta0 <- c(object$estiP$beta)
    pred <- Y - (Ki %*% (Y - beta0)) / diag(Ki)
    sd2 <- 1 / diag(Ki) - rowSums(Ki)^2 / (diag(Ki)^2 * sum(Ki))
    rmse <- mean((pred - Y)^2)
  }
  if (object$used_args$trend.type == "UK") {
    beta <- object$estiP$beta
    myF <- cbind(rep(1, nrow(X)), object$regF(X, t))
    pred <- Y - (Ki %*% (Y - myF %*% beta)) / diag(Ki)
    sd2 <- NULL
    rmse <- mean((pred - Y)^2)
  }
  return(list(pred = pred, sd2 = sd2, rmse = rmse))
}

#' Calculates the CRPS given the predictions (needs to be updated)
#'
#' @noRd

CRPS <- function(y, pred) {
  mu <- pred$mean
  sd <- pred$sd
  crps <- sd *
    (2 *
      dnorm((y - mu) / sd) +
      (y - mu) / sd * (2 * pnorm((y - mu) / sd) - 1) -
      1 / sqrt(pi))
  return(crps)
}
