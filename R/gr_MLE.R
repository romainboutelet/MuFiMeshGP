#' Binder function for gr_MLE
#'
#' @noRd

gr_MLE <- function(
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
    return(gr_MLE_SK(
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
    return(gr_RMLE_OK(
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
    return(gr_RMLE_UK(
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


#' Calculates the gradient of the MLE
#'
#' @noRd

gr_MLE_SK <- function(
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
  K1 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = 0,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K2 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = 0,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K <- K1 + K2
  Ki <- chol2inv(chol(K))
  DX <- distance.MuFiMeshGP(X)
  if (covtype == "Gaussian") {
    if (!iso && d > 1) {
      dK_phi1sq <- array(dim = c(n, n, d))
      dK_phi2sq <- array(dim = c(n, n, d))
      for (i in 1:d) {
        dK_phi1sq[,, i] <- (-DX[,, i]) * K1
        dK_phi2sq[,, i] <- (-DX[,, i]) * K2
      }
    } else {
      DX <- drop(DX)
      dK_phi1sq <- -DX * K1
      dK_phi2sq <- -DX * K2
    }
  }
  tmp <- crossprod(Ki, Y - mean.known)
  if (!iso && d > 1) {
    gr_nlogl_phi1sq <- rep(0, d)
    gr_nlogl_phi2sq <- rep(0, d)
    for (i in 1:d) {
      gr_nlogl_phi1sq[i] <- sum(fast_diag_cpp(Ki, dK_phi1sq[,, i])) -
        crossprod(tmp, crossprod(dK_phi1sq[,, i], tmp)) /
          crossprod(Y - mean.known, tmp)
      gr_nlogl_phi2sq[i] <- sum(fast_diag_cpp(Ki, dK_phi2sq[,, i])) -
        crossprod(tmp, crossprod(dK_phi2sq[,, i], tmp)) /
          crossprod(Y - mean.known, tmp)
    }
  } else {
    gr_nlogl_phi1sq <- sum(fast_diag_cpp(Ki, dK_phi1sq)) -
      crossprod(tmp, crossprod(dK_phi1sq, tmp)) /
        crossprod(Y - mean.known, tmp)
    gr_nlogl_phi2sq <- sum(fast_diag_cpp(Ki, dK_phi2sq)) -
      crossprod(tmp, crossprod(dK_phi2sq, tmp)) /
        crossprod(Y - mean.known, tmp)
  }
  dK_sigma2sq <- K2 / sigma2sq
  gr_nlogl_sigma2sq <- 1 /
    2 *
    (sum(fast_diag_cpp(Ki, dK_sigma2sq)) -
      n *
        crossprod(tmp, crossprod(dK_sigma2sq, tmp)) /
        crossprod(Y - mean.known, tmp))
  if (is.null(H.known)) {
    dK_H <- K2 * dK_H_cpp(t, H, l)
    gr_nlogl_H <- 1 /
      2 *
      (sum(fast_diag_cpp(Ki, dK_H)) -
        n *
          crossprod(tmp, crossprod(dK_H, tmp)) /
          crossprod(Y - mean.known, tmp))
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq,
      gr_nlogl_H
    ))
  } else {
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq
    ))
  }
}

#' Calculates the gradient of the RMLE
#'
#' @noRd

gr_RMLE_OK <- function(
  par,
  X,
  t,
  Y,
  myregF,
  l,
  covtype,
  H.known,
  iso = FALSE,
  nugget = NULL
) {
  d <- ncol(X)
  n <- nrow(X)
  p <- 1
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
  K1 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = 0,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K2 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = 0,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K <- K1 + K2
  Z <- crossprod(
    diag(x = 1, nrow = n) -
      tcrossprod(
        myregF %*%
          chol2inv(chol(crossprod(myregF, myregF))),
        myregF
      ),
    Y
  )
  Ki <- chol2inv(chol(K))
  KiF <- Ki %*% myregF
  W <- c(crossprod(myregF, KiF))
  Wi <- 1 / W
  KiFWi <- KiF * Wi
  tmp <- crossprod(Ki, Z)
  num <- crossprod(Z, (Ki - Wi * tcrossprod(KiF, KiF)) %*% Z)
  DX <- distance.MuFiMeshGP(X)
  if (!iso && d > 1) {
    if (covtype == "Gaussian") {
      dK_phi1sq <- array(dim = c(n, n, d))
      dK_phi2sq <- array(dim = c(n, n, d))
      for (i in 1:d) {
        dK_phi1sq[,, i] <- (-DX[,, i]) * K1
        dK_phi2sq[,, i] <- (-DX[,, i]) * K2
      }
    }
    dW_phi1sq <- rep(0, d)
    dW_phi2sq <- rep(0, d)
    for (i in 1:d) {
      dW_phi1sq[i] <- -crossprod(KiF, dK_phi1sq[,, i]) %*% KiF
      dW_phi2sq[i] <- -crossprod(KiF, dK_phi2sq[,, i]) %*% KiF
    }
    gr_nlogl_phi1sq <- rep(0, d)
    gr_nlogl_phi2sq <- rep(0, d)
    for (i in 1:d) {
      dnum_phi1sq <- 2 *
        crossprod(
          tmp,
          tcrossprod(
            dK_phi1sq[,, i] %*%
              KiFWi,
            KiF
          )
        ) %*%
          Z +
        crossprod(Z, KiFWi) %*% dW_phi1sq[i] %*% crossprod(KiFWi, Z) -
        crossprod(tmp, crossprod(dK_phi1sq[,, i], tmp))
      gr_nlogl_phi1sq[i] <- 1 /
        2 *
        (sum(fast_diag_cpp(Ki, dK_phi1sq[,, i])) +
          Wi * dW_phi1sq[i] +
          (n - p) * dnum_phi1sq / num)
      dnum_phi2sq <- 2 *
        crossprod(
          tmp,
          tcrossprod(
            dK_phi2sq[,, i] %*%
              KiFWi,
            KiF
          )
        ) %*%
          Z +
        crossprod(Z, KiFWi) %*% dW_phi2sq[i] %*% crossprod(KiFWi, Z) -
        crossprod(tmp, crossprod(dK_phi2sq[,, i], tmp))
      gr_nlogl_phi2sq[i] <- 1 /
        2 *
        (sum(fast_diag_cpp(Ki, dK_phi2sq[,, i])) +
          Wi * dW_phi2sq[i] +
          (n - p) * dnum_phi2sq / num)
    }
  } else {
    DX <- drop(DX)
    if (covtype == "Gaussian") {
      dK_phi1sq <- -DX * K1
      dK_phi2sq <- -DX * K2
      dW_phi1sq <- c(-crossprod(KiF, dK_phi1sq) %*% KiF)
      dW_phi2sq <- c(-crossprod(KiF, dK_phi2sq) %*% KiF)
    }
    dnum_phi1sq <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_phi1sq %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_phi1sq %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_phi1sq, tmp))
    gr_nlogl_phi1sq <- 1 /
      2 *
      (sum(fast_diag_cpp(Ki, dK_phi1sq)) +
        Wi * dW_phi1sq +
        (n - p) * dnum_phi1sq / num)
    dnum_phi2sq <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_phi2sq %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_phi2sq %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_phi2sq, tmp))
    gr_nlogl_phi2sq <- 1 /
      2 *
      (fast_trace_cpp(Ki, dK_phi2sq) +
        Wi * dW_phi2sq +
        (n - p) * dnum_phi2sq / num)
  }
  dK_sigma1sq <- K1 / sigma1sq
  dW_sigma1sq <- -c(crossprod(KiF, dK_sigma1sq) %*% KiF)
  dnum_sigma1sq <- 2 *
    crossprod(
      tmp,
      tcrossprod(
        dK_sigma1sq %*%
          KiFWi,
        KiF
      )
    ) %*%
      Z +
    crossprod(Z, KiFWi) %*% dW_sigma1sq %*% crossprod(KiFWi, Z) -
    crossprod(tmp, crossprod(dK_sigma1sq, tmp))
  gr_nlogl_sigma1sq <- 1 /
    2 *
    (sum(fast_diag_cpp(Ki, dK_sigma1sq)) +
      Wi * dW_sigma1sq +
      (n - p) * dnum_sigma1sq / num)
  dK_sigma2sq <- K2 / sigma2sq
  dW_sigma2sq <- -c(crossprod(KiF, dK_sigma2sq) %*% KiF)
  dnum_sigma2sq <- 2 *
    crossprod(
      tmp,
      tcrossprod(
        dK_sigma2sq %*%
          KiFWi,
        KiF
      )
    ) %*%
      Z +
    crossprod(Z, KiFWi) %*% dW_sigma2sq %*% crossprod(KiFWi, Z) -
    crossprod(tmp, crossprod(dK_sigma2sq, tmp))
  gr_nlogl_sigma2sq <- 1 /
    2 *
    (sum(fast_diag_cpp(Ki, dK_sigma2sq)) +
      Wi * dW_sigma2sq +
      (n - p) * dnum_sigma2sq / num)
  if (is.null(H.known)) {
    dK_H <- K2 * dK_H_cpp(t, H, l)
    dW_H <- -c(crossprod(KiF, dK_H) %*% KiF)
    dnum_H <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_H %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_H %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_H, tmp))
    gr_nlogl_H <- 1 /
      2 *
      (sum(fast_diag_cpp(Ki, dK_H)) +
        Wi * dW_H +
        (n - p) * dnum_H / num)
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq,
      gr_nlogl_H
    ))
  } else {
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq
    ))
  }
}

#' Calculates the gradient of the RMLE
#'
#' @noRd

gr_RMLE_UK <- function(
  par,
  X,
  t,
  Y,
  myregF,
  l,
  covtype,
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
  K1 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = 0,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K2 <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = 0,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  K <- K1 + K2
  Z <- crossprod(
    diag(x = 1, nrow = n) -
      tcrossprod(
        myregF %*%
          chol2inv(chol(crossprod(myregF, myregF))),
        myregF
      ),
    Y
  )
  Ki <- chol2inv(chol(K))
  KiF <- Ki %*% myregF
  W <- crossprod(myregF, KiF)
  Wi <- chol2inv(chol(W))
  KiFWi <- KiF %*% Wi
  tmp <- crossprod(Ki, Z)
  num <- crossprod(Z, (Ki - KiF %*% tcrossprod(Wi, KiF)) %*% Z)
  DX <- distance.MuFiMeshGP(X)
  if (!iso && d > 1) {
    if (covtype == "Gaussian") {
      dK_phi1sq <- array(dim = c(n, n, d))
      dK_phi2sq <- array(dim = c(n, n, d))
      for (i in 1:d) {
        dK_phi1sq[,, i] <- (-DX[,, i]) * K1
        dK_phi2sq[,, i] <- (-DX[,, i]) * K2
      }
    }
    dW_phi1sq <- array(dim = c(p, p, d))
    dW_phi2sq <- array(dim = c(p, p, d))
    for (i in 1:d) {
      dW_phi1sq[,, i] <- -crossprod(KiF, dK_phi1sq[,, i]) %*% KiF
      dW_phi2sq[,, i] <- -crossprod(KiF, dK_phi2sq[,, i]) %*% KiF
    }
    gr_nlogl_phi1sq <- rep(0, d)
    gr_nlogl_phi2sq <- rep(0, d)
    for (i in 1:d) {
      dnum_phi1sq <- 2 *
        crossprod(
          tmp,
          tcrossprod(
            dK_phi1sq[,, i] %*%
              KiFWi,
            KiF
          )
        ) %*%
          Z +
        crossprod(Z, KiFWi) %*% dW_phi1sq[,, i] %*% crossprod(KiFWi, Z) -
        crossprod(tmp, crossprod(dK_phi1sq[,, i], tmp))
      gr_nlogl_phi1sq[i] <- 1 /
        2 *
        (sum(fast_diag_cpp(Ki, dK_phi1sq[,, i])) +
          sum(fast_diag_cpp(Wi, dW_phi1sq[,, i])) +
          (n - p) * dnum_phi1sq / num)
      dnum_phi2sq <- 2 *
        crossprod(
          tmp,
          tcrossprod(
            dK_phi2sq[,, i] %*%
              KiFWi,
            KiF
          )
        ) %*%
          Z +
        crossprod(Z, KiFWi) %*% dW_phi2sq[,, i] %*% crossprod(KiFWi, Z) -
        crossprod(tmp, crossprod(dK_phi2sq[,, i], tmp))
      gr_nlogl_phi2sq[i] <- 1 /
        2 *
        (sum(fast_diag_cpp(Ki, dK_phi2sq[,, i])) +
          sum(fast_diag_cpp(Wi, dW_phi2sq[,, i])) +
          (n - p) * dnum_phi2sq / num)
    }
  } else {
    DX <- drop(DX)
    if (covtype == "Gaussian") {
      dK_phi1sq <- -DX * K1
      dK_phi2sq <- -DX * K2
    }
    dW_phi1sq <- -crossprod(KiF, dK_phi1sq) %*% KiF
    dW_phi2sq <- -crossprod(KiF, dK_phi2sq) %*% KiF
    dnum_phi1sq <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_phi1sq %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_phi1sq %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_phi1sq, tmp))
    gr_nlogl_phi1sq <- 1 /
      2 *
      (sum(fast_diag_cpp(Ki, dK_phi1sq)) +
        sum(fast_diag_cpp(Wi, dW_phi1sq)) +
        (n - p) * dnum_phi1sq / num)
    dnum_phi2sq <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_phi2sq %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_phi2sq %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_phi2sq, tmp))
    gr_nlogl_phi2sq <- 1 /
      2 *
      (sum(fast_diag_cpp(Ki, dK_phi2sq)) +
        sum(fast_diag_cpp(Wi, dW_phi2sq)) +
        (n - p) * dnum_phi2sq / num)
  }
  dK_sigma1sq <- K1 / sigma1sq
  dW_sigma1sq <- -crossprod(KiF, dK_sigma1sq) %*% KiF
  dnum_sigma1sq <- 2 *
    crossprod(
      tmp,
      tcrossprod(
        dK_sigma1sq %*%
          KiFWi,
        KiF
      )
    ) %*%
      Z +
    crossprod(Z, KiFWi) %*% dW_sigma1sq %*% crossprod(KiFWi, Z) -
    crossprod(tmp, crossprod(dK_sigma1sq, tmp))
  gr_nlogl_sigma1sq <- 1 /
    2 *
    (sum(fast_diag_cpp(Ki, dK_sigma1sq)) +
      sum(fast_diag_cpp(Wi, dW_sigma1sq)) +
      (n - p) * dnum_sigma1sq / num)
  dK_sigma2sq <- K2 / sigma2sq
  dW_sigma2sq <- -crossprod(KiF, dK_sigma2sq) %*% KiF
  dnum_sigma2sq <- 2 *
    crossprod(
      tmp,
      tcrossprod(
        dK_sigma2sq %*%
          KiFWi,
        KiF
      )
    ) %*%
      Z +
    crossprod(Z, KiFWi) %*% dW_sigma2sq %*% crossprod(KiFWi, Z) -
    crossprod(tmp, crossprod(dK_sigma2sq, tmp))
  gr_nlogl_sigma2sq <- 1 /
    2 *
    (sum(fast_diag_cpp(Ki, dK_sigma2sq)) +
      sum(fast_diag_cpp(Wi, dW_sigma2sq)) +
      (n - p) * dnum_sigma2sq / num)
  if (is.null(H.known)) {
    dK_H <- K2 * dK_H_cpp(t, H, l)
    dW_H <- -crossprod(KiF, dK_H) %*% KiF
    dnum_H <- 2 *
      crossprod(
        tmp,
        tcrossprod(
          dK_H %*%
            KiFWi,
          KiF
        )
      ) %*%
        Z +
      crossprod(Z, KiFWi) %*% dW_H %*% crossprod(KiFWi, Z) -
      crossprod(tmp, crossprod(dK_H, tmp))
    gr_nlogl_H <- 1 /
      2 *
      (fast_trace_cpp(Ki, dK_H) +
        fast_trace_cpp(Wi, dW_H) +
        (n - p) * dnum_H / num)
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq,
      gr_nlogl_H
    ))
  } else {
    return(c(
      phi1sq * gr_nlogl_phi1sq,
      phi2sq * gr_nlogl_phi2sq,
      sigma2sq * gr_nlogl_sigma2sq
    ))
  }
}
