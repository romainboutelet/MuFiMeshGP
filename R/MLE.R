#' Binder function for MLE
#'
#' @noRd

MLE <- function(
  trend.type,
  par,
  X,
  t,
  Y,
  l,
  covtype,
  myregF,
  mean.known,
  H.known,
  iso = FALSE,
  nugget = NULL
) {
  if (covtype == "Matern5_2") {
    stop("'Matern5_2' is not available at the moment")
  } else if (covtype == "Matern3_2") {
    stop("'Matern3_2' is not available at the moment")
  }
  if (trend.type == "SK") {
    return(MLE_SK(
      par = par,
      X = X,
      t = t,
      Y = Y,
      l = l,
      covtype = covtype,
      mean.known = mean.known,
      H.known = H.known,
      iso = iso,
      nugget = nugget
    ))
  } else if (trend.type == "OK") {
    return(RMLE_OK(
      par = par,
      X = X,
      t = t,
      Y = Y,
      l = l,
      covtype = covtype,
      myregF = myregF,
      H.known = H.known,
      iso = iso,
      nugget = nugget
    ))
  } else if (trend.type == "UK") {
    return(RMLE_UK(
      par = par,
      X = X,
      t = t,
      Y = Y,
      l = l,
      covtype = covtype,
      myregF = myregF,
      H.known = H.known,
      iso = iso,
      nugget = nugget
    ))
  }
}


#' Calculates the MLE
#'
#' @noRd

MLE_SK <- function(
  par,
  X,
  t,
  Y,
  l,
  covtype,
  mean.known,
  H.known,
  iso = FALSE,
  nugget = NULL
) {
  d <- ncol(X)
  n <- nrow(X)
  if (is.null(H.known)) {
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
      H <- par[2 * d + 2]
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
      H <- par[4]
    }
  } else {
    H <- H.known
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
    }
  }
  myK <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  myKi <- chol2inv(chol(myK))
  myval <- 1 /
    2 *
    as.double(determinant(myK)$modulus) +
    n / 2 * log(crossprod(Y - mean.known, myKi %*% (Y - mean.known)))
  return(c(myval))
}

#' Calculates the RMLE
#'
#' @noRd

RMLE_OK <- function(
  par,
  X,
  t,
  Y,
  l,
  covtype,
  myregF,
  H.known,
  iso = FALSE,
  nugget = NULL
) {
  d <- ncol(X)
  p <- 1
  n <- nrow(X)
  if (is.null(H.known)) {
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
      H <- par[2 * d + 2]
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
      H <- par[4]
    }
  } else {
    H <- H.known
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
    }
  }
  # Create contrasts
  Z <- crossprod(
    diag(x = 1, nrow = n) -
      tcrossprod(
        myregF %*%
          chol2inv(chol(crossprod(myregF, myregF))),
        myregF
      ),
    Y
  )
  myK <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  myKi <- chol2inv(chol(myK))
  myKiF <- myKi %*% myregF
  myW <- crossprod(myregF, myKiF)
  myWi <- chol2inv(chol(myW))
  myval <- 1 /
    2 *
    as.double(determinant(myK)$modulus) +
    1 / 2 * as.double(determinant(myW)$modulus) +
    (n - p) /
      2 *
      log(crossprod(
        Z,
        (myKi -
          myKiF %*%
            tcrossprod(myWi, myKiF)) %*%
          Z
      ))
  return(c(myval))
}

#' Calculates the RMLE
#'
#' @noRd

RMLE_UK <- function(
  par,
  X,
  t,
  Y,
  l,
  covtype,
  myregF,
  H.known,
  iso = FALSE,
  nugget = NULL
) {
  d <- ncol(X)
  n <- nrow(X)
  p <- ncol(myregF)
  if (is.null(H.known)) {
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
      H <- par[2 * d + 2]
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
      H <- par[4]
    }
  } else {
    H <- H.known
    if (!iso) {
      phi1sq <- exp(par[1:d])
      phi2sq <- exp(par[(d + 1):(2 * d)])
      sigma1sq <- 1
      sigma2sq <- exp(par[2 * d + 1])
    } else {
      phi1sq <- exp(par[1])
      phi2sq <- exp(par[2])
      sigma1sq <- 1
      sigma2sq <- exp(par[3])
    }
  }
  # Create contrasts
  Z <- crossprod(
    diag(x = 1, nrow = n) -
      tcrossprod(
        myregF %*%
          chol2inv(chol(crossprod(myregF, myregF))),
        myregF
      ),
    Y
  )
  myK <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  myKi <- solve(myK)
  myKiF <- myKi %*% myregF
  myW <- crossprod(myregF, myKiF)
  myWi <- chol2inv(chol(myW))
  myval <- 1 /
    2 *
    as.double(determinant(myK)$modulus) +
    1 / 2 * as.double(determinant(myW)$modulus) +
    (n - p) /
      2 *
      log(crossprod(
        Z,
        (myKi -
          myKiF %*%
            tcrossprod(myWi, myKiF)) %*%
          Z
      ))
  return(c(myval))
}
