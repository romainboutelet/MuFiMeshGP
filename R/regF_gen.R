#' Creates the regression function for the mean
#'
#' @description Creates the regression function for the GP mean.
#'
#' @param trend.dim which dimension should the trend follow: \code{"input"},
#' \code{"fidelity"}, or \code{"both"}.
#' @param trend.pol Which polynomial degree should the input mean trend have:
#' \code{"linear"} or \code{"quadratic"}.
#' @param interaction polynomial degree of the interaction between input trend
#' and fidelity trend of: \code{NULL}, \code{"linear"}, or \code{"quadratic"}.
#' \code{"linear"} or \code{"quadratic"}.
#' @param l convergence rate parameter, usually \code{l = 4}.
#' @param d input space dimension
#'
#' @return a function
#'
#' @export

regF_gen <- function(trend.dim, trend.pol, interaction, l, d) {
  if (trend.dim == "input") {
    if (is.null(trend.pol)) {
      if (is.null(interaction)) {
        stop("OK must be used if trend.pol and interaction are NULL.")
      } else if (interaction == "linear") {
        regF <- function(x, t) return((2 * x - 1) * t^(l / 2))
      } else if (interaction == "quadratic") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2),
              spat.int(x) * t^(l / 2)
            )) else
            return(cbind(
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2)
            ))
        }
      } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
    } else if (trend.pol == "linear") {
      if (is.null(interaction)) {
        regF <- function(x, t) return(2 * x - 1)
      } else if (interaction == "linear") {
        regF <- function(x, t) return(cbind(2 * x - 1, (2 * x - 1) * t^(l / 2)))
      } else if (interaction == "quadratic") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2),
              spat.int(x) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2)
            ))
        }
      } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
    } else if (trend.pol == "quadratic") {
      if (is.null(interaction)) {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(2 * x - 1, 6 * x^2 - 6 * x + 1, spat.int(x))) else
            return(cbind(2 * x - 1, 6 * x^2 - 6 * x + 1))
        }
      } else if (interaction == "linear") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              spat.int(x),
              (2 * x - 1) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              (2 * x - 1) * t^(l / 2)
            ))
        }
      } else if (interaction == "quadratic") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              spat.int(x),
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2),
              spat.int(x) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2)
            ))
        }
      } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
    } else stop("trend.pol must be: NULL, 'linear', or 'quadratic'")
  } else if (trend.dim == "fidelity") {
    if (is.null(interaction)) {
      regF <- function(x, t) return(cbind(t^(l / 2)))
    } else if (interaction == "linear") {
      regF <- function(x, t) return(cbind(t^(l / 2), (2 * x - 1) * t^(l / 2)))
    } else if (interaction == "quadratic") {
      regF <- function(x, t) {
        d <- ncol(x)
        if (d >= 2)
          return(cbind(
            t^(l / 2),
            (2 * x - 1) * t^(l / 2),
            (6 * x^2 - 6 * x + 1) * t^(l / 2),
            spat.int(x) * t^(l / 2)
          )) else
          return(cbind(
            t^(l / 2),
            (2 * x - 1) * t^(l / 2),
            (6 * x^2 - 6 * x + 1) * t^(l / 2)
          ))
      }
    } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
  } else if (trend.dim == "both") {
    if (is.null(trend.pol))
      stop(
        "Please use trend.dim = 'fidelity' instead of trend.dim = 'both' and trend.pol = NULL"
      )
    if (trend.pol == "linear") {
      if (is.null(interaction)) {
        regF <- function(x, t) return(cbind(2 * x - 1, t^(l / 2)))
      } else if (interaction == "linear") {
        regF <- function(x, t)
          return(cbind(2 * x - 1, t^(l / 2), (2 * x - 1) * t^(l / 2)))
      } else if (interaction == "quadratic") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              t^(l / 2),
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2),
              spat.int(x) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              t^(l / 2),
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2)
            ))
        }
      } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
    } else if (trend.pol == "quadratic") {
      if (is.null(interaction)) {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              spat.int(x),
              t^(l / 2)
            )) else return(cbind(2 * x - 1, 6 * x^2 - 6 * x + 1, t^(l / 2)))
        }
      } else if (interaction == "linear") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              spat.int(x),
              t^(l / 2),
              (2 * x - 1) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              t^(l / 2),
              (2 * x - 1) * t^(l / 2)
            ))
        }
      } else if (interaction == "quadratic") {
        regF <- function(x, t) {
          d <- ncol(x)
          if (d >= 2)
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              spat.int(x),
              t^(l / 2),
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2),
              spat.int(x) * t^(l / 2)
            )) else
            return(cbind(
              2 * x - 1,
              6 * x^2 - 6 * x + 1,
              t^(l / 2),
              (2 * x - 1) * t^(l / 2),
              (6 * x^2 - 6 * x + 1) * t^(l / 2)
            ))
        }
      } else stop("interaction must be: NULL, 'linear', or 'quadratic'")
    } else
      stop("trend.pol must be: 'linear', or 'quadratic' for trend.dim = 'both'")
  }
  return(regF)
}

#' Calculates the interation between fidelity and input
#'
#' @noRd

spat.int <- function(x) {
  spat_int <- matrix(NA, nrow(x), ncol = choose(ncol(x), 2))
  ind <- matrix(NA, nrow = ncol(spat_int), ncol = 2)
  count <- 1
  for (i in 1:(ncol(x) - 1)) {
    ind[count:(count + ncol(x) - i - 1), ] <- cbind(
      rep(i, (ncol(x) - i)),
      (i + 1):(ncol(x))
    )
    count <- count + ncol(x) - i
  }
  for (i in 1:ncol(spat_int)) {
    spat_int[, i] <- (2 * x[, ind[i, 1]] - 1) * (2 * x[, ind[i, 2]] - 1)
  }
  return(spat_int)
}
