#' Calculate W for the IMSPE
#'
#' @noRd

Wij <- function(x1, x2 = NULL, phi1sq, sigma1sq, covtype) {
  if (covtype == "Gaussian") {
    if (is.null(x2)) {
      return(Wijs1_gauss_cpp(x = x1, phi1sq = phi1sq, sigma1sq = sigma1sq))
    } else
      return(Wijs2_gauss_cpp(
        x1 = x1,
        x2 = x2,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq
      ))
  }
}

#' Calculate H for the IMSPE
#'
#' @noRd

Hij <- function(
  mu,
  phi1sq,
  sigma1sq,
  covtype,
  trend.type,
  trend.dim,
  trend.pol,
  interaction
) {
  d <- ncol(mu)
  n <- nrow(mu)
  if (covtype == "Gaussian") {
    Phi1 <- matrix(phi1sq, nrow = n, ncol = d, byrow = T)
    I1 <- sqrt(pi / Phi1) *
      (pnorm(sqrt(2 * Phi1) * mu) -
        pnorm(sqrt(2 * Phi1) * (mu - 1)))
    h1 <- sigma1sq * matrix(apply(I1, 1, 'prod'), ncol = 1)
    if (trend.type == "OK") return(h1)
    if (trend.dim == "fidelity") {
      if (is.null(interaction)) {
        return(cbind(h1, matrix(rep(0, n), ncol = 1)))
      } else if (interaction == "linear") {
        return(cbind(h1, matrix(rep(0, (d + 1) * n), ncol = d + 1)))
      } else if (interaction == "quadratic") {
        return(cbind(
          h1,
          matrix(
            rep(0, (2 * d + 1 + choose(d, 2)) * n),
            ncol = 2 * d + 1 + choose(d, 2)
          )
        ))
      }
    } else if (trend.dim == "input") {
      if (is.null(trend.pol)) {
        if (interaction == "linear") {
          return(cbind(h1, matrix(rep(0, d * n), ncol = d)))
        } else if (interaction == "quadratic") {
          return(cbind(
            h1,
            matrix(
              rep(0, (2 * d + choose(d, 2)) * n),
              ncol = 2 * d + choose(d, 2)
            )
          ))
        }
      }
    }
    I2 <- (exp(-Phi1 * mu^2) - exp(-Phi1 * (1 - mu)^2))
    I3 <- (mu - 1) * exp(-Phi1 * (1 - mu)^2) - mu * exp(-Phi1 * mu^2) + I1

    h2 <- matrix(NA, nrow = n, ncol = d)
    for (k in 1:d)
      if (d > 1)
        h2[, k] <- (I2[, k] / (2 * Phi1[, k]) + mu[, k] * I1[, k]) *
          matrix(apply(I1[, -k, drop = F], 1, 'prod'), ncol = 1) else
        h2[, k] <- (I2[, k] / (2 * Phi1[, k]) + mu[, k] * I1[, k])
    h2 <- sigma1sq * h2
    if (trend.dim == "input") {
      if (trend.pol == "linear") {
        if (is.null(interaction)) {
          return(cbind(h1, 2 * h2 - matrix(rep(h1, d), ncol = d)))
        } else if (interaction == "linear") {
          return(cbind(
            h1,
            2 * h2 - matrix(rep(h1, d), ncol = d),
            matrix(rep(0, d * n), ncol = d)
          ))
        } else if (interaction == "quadratic") {
          return(cbind(
            h1,
            2 * h2 - matrix(rep(h1, d), ncol = d),
            matrix(
              rep(0, (2 * d + choose(d, 2)) * n),
              ncol = 2 * d + choose(d, 2)
            )
          ))
        }
      }
    } else if (trend.dim == "both") {
      if (trend.pol == "linear") {
        if (is.null(interaction)) {
          return(cbind(h1, 2 * h2 - matrix(rep(h1, d), ncol = d)))
        } else if (interaction == "linear") {
          return(cbind(
            h1,
            2 * h2 - matrix(rep(h1, d), ncol = d),
            matrix(rep(0, (d + 1) * n), ncol = d + 1)
          ))
        } else if (interaction == "quadratic") {
          return(cbind(
            h1,
            2 * h2 - matrix(rep(h1, d), ncol = d),
            matrix(
              rep(0, (2 * d + 1 + choose(d, 2)) * n),
              ncol = 2 * d + 1 + choose(d, 2)
            )
          ))
        }
      }
    }
    h3 <- matrix(NA, nrow = n, ncol = d)
    if (d > 1) {
      for (k in 1:d)
        h3[, k] <- (I3[, k] /
          (2 * Phi1[, k]) +
          mu[, k] * I2[, k] / Phi1[, k] +
          mu[, k]^2 * I1[, k]) *
          matrix(apply(I1[, -k, drop = F], 1, 'prod'), ncol = 1)
    } else {
      for (k in 1:d)
        h3[, k] <- I3[, k] /
          (2 * Phi1[, k]) +
          mu[, k] * I2[, k] / Phi1[, k] +
          mu[, k]^2 * I1[, k]
    }
    h3 <- sigma1sq * h3
    if (d > 1) {
      h4 <- matrix(NA, nrow = n, ncol = choose(d, 2))
      ind <- matrix(NA, nrow = choose(d, 2), ncol = 2)
      count <- 1
      for (i in 1:(d - 1)) {
        ind[count:(count + d - i - 1), ] <- cbind(rep(i, (d - i)), (i + 1):(d))
        count <- count + d - i
      }
      for (i in 1:choose(d, 2)) {
        h4[, i] <- (I2[, ind[i, 1]] /
          (2 * Phi1[, ind[i, 1]]) +
          mu[, ind[i, 1]] * I1[, ind[i, 1]]) *
          (I2[, ind[i, 2]] /
            (2 * Phi1[, ind[i, 2]]) +
            mu[, ind[i, 2]] * I1[, ind[i, 2]]) *
          matrix(apply(I1[, -ind[i, ], drop = F], 1, 'prod'), ncol = 1)
      }
      h4 <- sigma1sq * h4
    }
    if (trend.dim == "input") {
      if (trend.pol == "quadratic") {
        if (is.null(interaction)) {
          if (d > 1) {
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4
            ))
          } else {
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d)
            ))
          }
        } else if (interaction == "linear") {
          if (d > 1)
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4,
              matrix(rep(0, d * n), ncol = d)
            )) else
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              matrix(rep(0, d * n), ncol = d)
            ))
        } else if (interaction == "quadratic") {
          if (d > 1)
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4,
              matrix(
                rep(0, (2 * d + choose(d, 2)) * n),
                ncol = 2 * d + choose(d, 2)
              )
            )) else
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              matrix(rep(0, (2 * d) * n), ncol = 2 * d)
            ))
        }
      }
    } else if (trend.dim == "both") {
      if (trend.pol == "quadratic") {
        if (is.null(interaction)) {
          if (d > 1) {
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4,
              matrix(rep(0, n), ncol = 1)
            ))
          } else {
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              matrix(rep(0, n), ncol = 1)
            ))
          }
        } else if (interaction == "linear") {
          if (d > 1)
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4,
              matrix(rep(0, (d + 1) * n), ncol = d + 1)
            )) else
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              matrix(rep(0, (d + 1) * n), ncol = d + 1)
            ))
        } else if (interaction == "quadratic") {
          if (d > 1)
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              h4,
              matrix(
                rep(0, (2 * d + 1 + choose(d, 2)) * n),
                ncol = 2 * d + 1 + choose(d, 2)
              )
            )) else
            return(cbind(
              h1,
              2 * h2 - matrix(rep(h1, d), ncol = d),
              6 * h3 - 6 * h2 + matrix(rep(h1, d), ncol = d),
              matrix(rep(0, (2 * d + 1) * n), ncol = 2 * d + 1)
            ))
        }
      }
    }
  }
}


#' Calculate G for the IMSPE
#'
#' @noRd

Gij <- function(d, trend.type, trend.dim, trend.pol, interaction) {
  if (trend.type == "OK") {
    G <- matrix(data = 1, nrow = 1, ncol = 1)
  } else if (trend.type == "UK") {
    if (trend.dim == "fidelity") {
      if (is.null(interaction)) {
        G <- matrix(data = 0, nrow = 2, ncol = 2)
        G[1, 1] <- 1
      } else if (interaction == "linear") {
        G <- matrix(data = 0, nrow = d + 2, ncol = d + 2)
        G[1, 1] <- 1
      } else if (interaction == "quadratic") {
        G <- matrix(
          data = 0,
          nrow = 2 * d + choose(d, 2) + 2,
          ncol = 2 * d + choose(d, 2) + 2
        )
        G[1, 1] <- 1
      }
    } else if (trend.dim == "input") {
      if (is.null(trend.pol)) {
        if (interaction == "linear") {
          G <- matrix(data = 0, nrow = d + 1, ncol = d + 1)
        } else if (interaction == "quadratic") {
          G <- matrix(
            data = 0,
            nrow = 2 * d + choose(d, 2) + 1,
            ncol = 2 * d + choose(d, 2) + 1
          )
        }
        G[1, 1] <- 1
      } else if (trend.pol == "linear") {
        if (is.null(interaction)) {
          G <- matrix(data = 0, nrow = d + 1, ncol = d + 1)
        } else if (interaction == "linear") {
          G <- matrix(data = 0, nrow = 2 * d + 1, ncol = 2 * d + 1)
        } else if (interaction == "quadratic") {
          G <- matrix(
            data = 0,
            nrow = 3 * d + choose(d, 2) + 1,
            ncol = 3 * d + choose(d, 2) + 1
          )
        }
        G[1, 1] <- 1
        if (d > 1) diag(G[2:(d + 1), 2:(d + 1)]) <- 1 / 3 else
          G[2:(d + 1), 2:(d + 1)] <- 1 / 3
      } else if (trend.pol == "quadratic") {
        if (is.null(interaction)) {
          G <- matrix(
            data = 0,
            nrow = 2 * d + choose(d, 2) + 1,
            ncol = 2 * d + choose(d, 2) + 1
          )
        } else if (interaction == "linear") {
          G <- matrix(
            data = 0,
            nrow = 3 * d + choose(d, 2) + 1,
            ncol = 3 * d + choose(d, 2) + 1
          )
        } else if (interaction == "quadratic") {
          G <- matrix(
            data = 0,
            nrow = 4 * d + 2 * choose(d, 2) + 1,
            ncol = 4 * d + 2 * choose(d, 2) + 1
          )
        }
        G[1, 1] <- 1
        if (d > 1) diag(G[2:(d + 1), 2:(d + 1)]) <- 1 / 3 else
          G[2:(d + 1), 2:(d + 1)] <- 1 / 3
        if (d > 1)
          diag(G[(d + 2):(2 * d + 1), (d + 2):(2 * d + 1)]) <- 1 / 5 else
          G[(d + 2):(2 * d + 1), (d + 2):(2 * d + 1)] <- 1 / 5
        if (d > 1) {
          if (d == 2) {
            G[
              (2 * d + 2):(2 * d + 1 + choose(d, 2)),
              (2 * d + 2):(2 * d + 1 + choose(d, 2))
            ] <- 1 / 9
          } else {
            diag(G[
              (2 * d + 2):(2 * d + 1 + choose(d, 2)),
              (2 * d + 2):(2 * d + 1 + choose(d, 2))
            ]) <- 1 / 9
          }
        }
      }
    } else if (trend.dim == "both") {
      if (trend.pol == "linear") {
        if (is.null(interaction)) {
          G <- matrix(data = 0, nrow = d + 2, ncol = d + 2)
        } else if (interaction == "linear") {
          G <- matrix(data = 0, nrow = 2 * d + 2, ncol = 2 * d + 2)
        } else if (interaction == "quadratic") {
          G <- matrix(
            data = 0,
            nrow = 3 * d + choose(d, 2) + 2,
            ncol = 3 * d + choose(d, 2) + 2
          )
        }
        G[1, 1] <- 1
        if (d > 1) diag(G[2:(d + 1), 2:(d + 1)]) <- 1 / 3 else
          G[2:(d + 1), 2:(d + 1)] <- 1 / 3
      } else if (trend.pol == "quadratic") {
        if (is.null(interaction)) {
          G <- matrix(
            data = 0,
            nrow = 2 * d + choose(d, 2) + 2,
            ncol = 2 * d + choose(d, 2) + 2
          )
        } else if (interaction == "linear") {
          G <- matrix(
            data = 0,
            nrow = 3 * d + choose(d, 2) + 2,
            ncol = 3 * d + choose(d, 2) + 2
          )
        } else if (interaction == "quadratic") {
          G <- matrix(
            data = 0,
            nrow = 4 * d + 2 * choose(d, 2) + 2,
            ncol = 4 * d + 2 * choose(d, 2) + 2
          )
        }
        G[1, 1] <- 1
        if (d > 1) diag(G[2:(d + 1), 2:(d + 1)]) <- 1 / 3 else
          G[2:(d + 1), 2:(d + 1)] <- 1 / 3
        if (d > 1)
          diag(G[(d + 2):(2 * d + 1), (d + 2):(2 * d + 1)]) <- 1 / 5 else
          G[(d + 2):(2 * d + 1), (d + 2):(2 * d + 1)] <- 1 / 5
        if (d > 1) {
          if (d == 2) {
            G[
              (2 * d + 2):(2 * d + 1 + choose(d, 2)),
              (2 * d + 2):(2 * d + 1 + choose(d, 2))
            ] <- 1 / 9
          } else {
            diag(G[
              (2 * d + 2):(2 * d + 1 + choose(d, 2)),
              (2 * d + 2):(2 * d + 1 + choose(d, 2))
            ]) <- 1 / 9
          }
        }
      }
    }
  }
  return(G)
}


#' To be updated
#'
#' @noRd

dx_kn <- function(X, x, k1n, k2n, phi1sq, phi2sq, covtype) {
  if (covtype == "Gaussian")
    return(dx_kn_gauss_cpp(X, x, k1n, k2n, phi1sq, phi2sq))
}

#' To be updated
#'
#' @noRd

dt_sigmasq <- function(dt_kn, Ki, kn, dt_k, d) {
  if (d > 1)
    return(dt_sigmasq_cpp(dt_kn = dt_kn, Ki = Ki, kn = kn, dt_k = dt_k)) else
    return(dt_k - 2 * crossprod(dt_kn, crossprod(Ki, kn)))
}

#' To be updated
#'
#' @noRd

dx_wn <- function(x1, x2, wn, phi1sq, sigma1sq, covtype) {
  if (covtype == "Gaussian") {
    return(
      sigma1sq *
        c1_gauss_cpp(X = x1, x = x2, sigma = sqrt(1 / phi1sq), W = wn)
    )
  }
}

#' To be updated
#'
#' @noRd

dx_w <- function(x, w, phi1sq, sigma1sq, covtype) {
  if (covtype == "Gaussian") {
    return(sigma1sq * c2_gauss_cpp(x = x, t = sqrt(1 / phi1sq), w = w))
  }
}
