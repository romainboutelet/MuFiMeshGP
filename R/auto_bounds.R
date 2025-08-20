#' Creates automatic bounds given design for decay parameter
#'
#' @param X design matrix
#'
#' @return a list with \code{lower} and \code{upper} bounds
#'
#' @noRd

auto_bounds <- function(
  X,
  min_cor = 0.01,
  max_cor = 0.5,
  covtype = "Gaussian",
  p = 0.05
) {
  Xsc <- find_reps(X, rep(1, nrow(X)), rescale = T)
  dists <- isoDist1_cpp(Xsc$X0)
  repr_low_dist <- quantile(
    x = dists[lower.tri(dists)],
    probs = p,
    names = FALSE
  )
  repr_lar_dist <- quantile(
    x = dists[lower.tri(dists)],
    probs = 1 -
      p,
    names = FALSE
  )
  if (covtype == "Gaussian") {
    theta_min <- -repr_low_dist / log(min_cor)
    theta_max <- -repr_lar_dist / log(max_cor)
    return(list(
      lower = theta_min *
        (Xsc$inputBounds[2, ] -
          Xsc$inputBounds[1, ])^2,
      upper = theta_max * (Xsc$inputBounds[2, ] - Xsc$inputBounds[1, ])^2
    ))
  } else {
    stop("Non Gaussian kernel is not implemented yet")
    #   tmpfun <- function(theta, repr_dist, covtype, value) {
    #     cov_auto(
    #       matrix(sqrt(repr_dist / ncol(X)), ncol = ncol(X)),
    #       matrix(0, ncol = ncol(X)),
    #       type = covtype,
    #       theta = theta
    #     ) -
    #       value
    #   }
    #   theta_min <- try(
    #     uniroot(
    #       tmpfun,
    #       interval = c(sqrt(.Machine$double.eps), 100),
    #       covtype = covtype,
    #       value = min_cor,
    #       repr_dist = repr_low_dist,
    #       tol = sqrt(.Machine$double.eps)
    #     )$root
    #   )
    #   if (is(theta_min, "try-error")) {
    #     warning(
    #       "The automatic selection of lengthscales bounds was not successful. Perhaps provide lower and upper values."
    #     )
    #     theta_min <- 0.01
    #   }
    #   theta_max <- try(
    #     uniroot(
    #       tmpfun,
    #       interval = c(sqrt(.Machine$double.eps), 100),
    #       covtype = covtype,
    #       value = max_cor,
    #       repr_dist = repr_lar_dist,
    #       tol = sqrt(.Machine$double.eps)
    #     )$root,
    #     silent = TRUE
    #   )
    #   if (is(theta_max, "try-error")) {
    #     theta_max <- 5
    #   }
    #   return(list(
    #     lower = theta_min *
    #       (Xsc$inputBounds[2, ] -
    #         Xsc$inputBounds[1, ]),
    #     upper = max(1, theta_max) *
    #       (Xsc$inputBounds[2, ] - Xsc$inputBounds[1, ])
    #   ))
  }
}

#' Creates covariance function for auto_bounds
#'
#' @noRd

cov_auto <- function(
  X1,
  X2 = NULL,
  theta,
  type = c("Gaussian", "Matern5_2", "Matern3_2")
) {
  type <- match.arg(type)
  if (type == "Gaussian")
    return(cov_gen(
      x1 = X1,
      x2 = X2,
      t1 = rep(0, nrow(X1)),
      t2 = rep(0, nrow(X2)),
      phi1sq = 1 / theta^2,
      phi2sq = rep(1, length(theta)),
      sigma1sq = 1,
      sigma2sq = 0
    ))
  if (type == "Matern5_2") stop("Not implemented yet")
  # return(cov_Matern5_2(X1 = X1, X2 = X2, theta = theta))
  if (type == "Matern3_2") stop("Not implemented yet")
  # return(cov_Matern3_2(X1 = X1, X2 = X2, theta = theta))
}

#' Finds repetitive design points
#'
#' @noRd

## All the functions in this file are based on code from the package hetGP,
## licensed under the LGPL-3 license

find_reps <- function(
  X,
  Z,
  return.Zlist = TRUE,
  rescale = FALSE,
  normalize = FALSE,
  inputBounds = NULL
) {
  if (is.null(dim(X))) X <- matrix(X, ncol = 1)
  if (nrow(X) == 1) {
    if (return.Zlist)
      return(list(X0 = X, Z0 = Z, mult = 1, Z = Z, Zlist = list(Z)))
    return(list(X0 = X, Z0 = Z, mult = 1, Z = Z))
  }
  if (rescale) {
    if (is.null(inputBounds)) inputBounds <- apply(X, 2, range)
    X <- (X -
      matrix(
        inputBounds[1, ],
        nrow = nrow(X),
        ncol = ncol(X),
        byrow = TRUE
      )) %*%
      diag(1 / (inputBounds[2, ] - inputBounds[1, ]), ncol(X))
  }
  outputStats <- NULL
  if (normalize) {
    outputStats <- c(mean(Z), var(Z))
    Z <- (Z - outputStats[1]) / sqrt(outputStats[2])
  }
  X0 <- unique(X)
  if (nrow(X) == nrow(X0)) {
    if (return.Zlist)
      return(list(
        X0 = X,
        Z0 = Z,
        mult = rep(1, length(Z)),
        Z = Z,
        Zlist = as.list(Z),
        inputBounds = inputBounds,
        outputStats = outputStats
      ))
    return(list(
      X0 = X,
      Z0 = Z,
      mult = rep(1, length(Z)),
      Z = Z,
      inputBounds = inputBounds,
      outputStats = outputStats
    ))
  }
  corresp <- find_corres(X0, X)
  Zlist <- split(Z, corresp)
  mult <- as.numeric(unlist(lapply(Zlist, length)))
  if (return.Zlist)
    return(list(
      X0 = X0,
      Z0 = unlist(lapply(Zlist, mean)),
      mult = mult,
      Z = unlist(Zlist),
      Zlist = Zlist,
      inputBounds = inputBounds,
      outputStats = outputStats
    ))
  return(list(
    X0 = X0,
    Z0 = unlist(lapply(Zlist, mean)),
    mult = mult,
    Z = unlist(Zlist),
    inputBounds = inputBounds,
    outputStats = outputStats
  ))
}
