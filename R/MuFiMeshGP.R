#' Prediction of the MuFiMeshGP emulator for any fidelity level.
#'
#' @description The function computes the posterior mean and standard deviation of the
#' MuFiMeshGP model.
#'
#' @seealso \code{\link{MuFiMeshGP}} for the model.
#'
#' @details From the model fitted by \code{\link{MuFiMeshGP}} or \code{\link{update.MuFiMeshGP}}
#'  the posterior mean and standard deviation are calculated for any input
#'  location and fidelity level.
#'  For details, see Boutelet and Sung (2025, <arXiv:2503.23158>).
#'
#' @param X matrix of input locations. Each row represents a sample.
#' @param t vector of fidelity levels. Each element is a sample and is connected
#'  to the corresponding row in \code{X}.
#' @param Y vector of response values.
#' @param covtype covariance kernel type, only 'Gaussian' is available for now,
#'  'Matern5_2' or 'Matern3_2' will be available soon (see \code{\link{cov_gen}}).
#' @param trend.type,trend.dim,trend.pol,interaction define the mean function form of
#'  the Gaussian process. \code{trend.type} can be: "SK" in which case
#'  \code{mean.known} needs to be specified as a scalar; "OK" in which case the constant
#'  mean will be evaluated through MLE; "UK" in which case \code{trend.dim}
#'  specifies whether the trend will be along the input space (\code{"input"}),
#'  the fidelity space (\code{"fidelity"}), or both (\code{"both"}).
#'  If \code{trend.dim} is \code{"input"} or \code{"both"}, the user can use the
#'  \code{trend.pol} to specify if the trend on the input space alone should be
#'  \code{"linear"} or \code{"quadratic}.
#'  Finally, if \code{trend.dim} is \code{"both"}, then an \code{interaction}
#'  term specify the polynomial order (\code{"linear} or \code{"quadratic"}) of
#'   the input space trend that is multiplied to the fidelity space trend.
#'   See \code{\link{regF_gen}} for further details.
#' @param mean.known Specifies the mean if \code{"SK"} as \code{trend.type},
#'  scalar.
#' @param H.known allow the user to specify the value of H as
#'  \code{H.known}, a scalar in (0,1).
#' @param gradient whether or not the gradient of the log-likelihood shouldbe
#'  used in the parameter estimation.
#' @param init Where should the parameter estimation start from, a vector.
#' @param single_fidelity can be used as \code{TRUE} to use \code{MuFiMeshGP}
#' as a single fidelity Gaussian Process. This will set \code{sigma2sq} as 0.
#' @param param.bounds a list with two arguments(\code{lower} and \code{upper})
#'  describing the bounds used for MLE optimization of \code{phi1sq} and \code{phi2sq}.
#' Each argument should be a vector of length \code{ncol(X)}.
#' If \code{NULL} the bounds of \code{phi1sq} and \code{phi2sq}
#' are specified automatically from the design matrix.
#' @param iso whether the covariance function will be isotropic (\code{TRUE}
#' or \code{FALSE})
#' @param l rate of convergence of the system (see Details), scalar.
#' @param nugget (optional) for controlling numerical error.
#' @param ncores (optional) number of cores for parallelization.
#'
#'
#' @return a list which is given the S3 class "MuFiMeshGP"
#' @importFrom methods is
#' @importFrom stats dnorm optim pnorm quantile var
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @useDynLib MuFiMeshGP, .registration = TRUE
#' @export
#' @examples
#' # Example code
#'
#' f <- function(x, t){
#'   x <- c(x)
#'   return(exp(-1.4*x)*cos(3.5*pi*x)+sin(40*x)/10*t^2)
#' }
#'
#' set.seed(1)
#' X <- matrix(runif(15,0,1), ncol = 1)
#' tt <- runif(15,0.5,2)
#'
#' Y <- f(c(X), tt)
#'
#' fit.mufimeshgp <- MuFiMeshGP(X, tt, Y)
#'
#' xx <- matrix(seq(0,1,0.01), ncol = 1)
#' ftrue <- f(xx, 0)
#'
#' # predict
#' pred.mufimeshgp <- predict(fit.mufimeshgp, xx, rep(0,101))
#'
#' mu <- pred.mufimeshgp$mean
#' s <- pred.mufimeshgp$sd
#' lower <- mu + qnorm(0.025)*s
#' upper <- mu + qnorm(0.975)*s
#'
#' # plot
#'
#' par(mfrow = c(1,1))
#' plot(xx, ftrue, "l", ylim = c(-1,1.3), ylab = "y", xlab = "x")
#' lines(c(xx), mu, col = "blue")
#' lines(c(xx), lower, col = "blue", lty = 2)
#' lines(c(xx), upper, col = "blue", lty = 2)
#' points(c(X), Y, col = "red")
#'
#'
#' ### RMSE ###
#' print(sqrt(mean((ftrue - mu))^2))
#'

MuFiMeshGP <- function(
  X,
  t,
  Y,
  covtype = "Gaussian",
  trend.type = "OK",
  trend.dim = "input",
  trend.pol = "quadratic",
  interaction = NULL,
  mean.known = NULL,
  H.known = NULL,
  gradient = TRUE,
  init = NULL,
  single_fidelity = F,
  param.bounds = NULL,
  iso = FALSE,
  l = 4,
  nugget = 1e-6,
  ncores = 1
) {
  eps <- sqrt(.Machine$double.eps)
  n <- nrow(X)
  d <- ncol(X)
  myregF <- NULL
  regF <- NULL
  Gijs <- NULL
  if (is.null(H.known)) myH <- 0.5 else myH <- H.known
  if (trend.type == "SK" & is.null(mean.known)) {
    stop("mean.known must be specified for SK")
  } else if (trend.type == "SK" & !is.null(mean.known)) {
    p <- 0
  } else if (trend.type == "OK") {
    myregF <- matrix(1, nrow = n, ncol = 1)
    trend.pol <- NULL
    interaction <- NULL
    p <- 1
    Gijs <- Gij(
      d = ncol(X),
      trend.type = trend.type,
      trend.dim = trend.dim,
      trend.pol = trend.pol,
      interaction = interaction
    )
  } else if (trend.type == "UK") {
    if (is.null(trend.dim)) stop("UK needs trend.dim to be non NULL")
    if (trend.dim %in% c("input", "both") & is.null(trend.pol))
      stop("trend.pol needs to be non NULL if trend.dim is 'input' or 'both'")
    regF <- regF_gen(
      trend.dim = trend.dim,
      trend.pol = trend.pol,
      interaction = interaction,
      l = l
    )
    myregF <- cbind(rep(1, n), regF(X, t))
    p <- dim(myregF)[2]
    Gijs <- Gij(
      d = ncol(X),
      trend.type = trend.type,
      trend.dim = trend.dim,
      trend.pol = trend.pol,
      interaction = interaction
    )
  }
  if (is.null(param.bounds)) {
    if (length(unique(t)) == 1) {
      bounds <- auto_bounds(X)
      single_fidelity <- T
    } else bounds <- auto_bounds(cbind(X, t))
    lower_phisq <- as.numeric(1 / bounds$upper[1:d])
    upper_phisq <- as.numeric(1 / bounds$lower[1:d])
  } else {
    tryCatch(
      expr = {
        lower_phisq <- param.bounds$lower
        upper_phisq <- param.bounds$upper
      },
      error = function(e) {
        print("param.bounds must be a list (see manual).")
      }
    )
  }
  sq.int <- max(Y) - min(Y)
  if (is.null(init)) {
    myphi1sq <- sqrt(lower_phisq * upper_phisq)
    myphi2sq <- myphi1sq
    mytausq <- 1 / mean(t)^l
  } else {
    if (!iso && d > 1) {
      myphi1sq <- init$phi1sq
      myphi2sq <- init$phi2sq
    } else {
      myphi1sq <- init$phi1sq[1]
      myphi2sq <- init$phi2sq[1]
    }
    mytausq <- init$sigma2sq / init$sigma1sq
    if (is.null(H.known)) myH <- init$H
  }
  envtmp <- environment()
  if (!single_fidelity) {
    fn <- function(par, env) {
      loglik <- MLE(
        trend.type = trend.type,
        par = par,
        X = X,
        t = t,
        Y = Y,
        l = l,
        covtype = covtype,
        myregF = myregF,
        mean.known = mean.known,
        H.known = H.known,
        iso = iso,
        nugget = nugget
      )
      if (!is.null(env) && !is.na(loglik)) {
        if (is.null(env$min_loglik) || loglik > env$min_loglik) {
          env$min_loglik <- loglik
          env$arg_min <- par
        }
      }
      return(loglik)
    }
    if (gradient) {
      gr <- function(par, env) {
        gr_loglik <- gr_MLE(
          trend.type = trend.type,
          par = par,
          X = X,
          t = t,
          Y = Y,
          l = l,
          covtype = covtype,
          myregF = myregF,
          mean.known = mean.known,
          H.known = H.known,
          iso = iso,
          nugget = nugget
        )
        return(gr_loglik)
      }
    } else gr <- NULL
    lower.params <- log(c(lower_phisq, lower_phisq, 1e-4 / max(t)^l))
    upper.params <- log(c(upper_phisq, upper_phisq, 1e4 / max(t)^l))
  } else {
    fn <- function(par, env) {
      loglik <- MLE(
        trend.type = trend.type,
        par = c(par, rep(1, d), 0, 0.5),
        X = X,
        t = t,
        Y = Y,
        l = l,
        covtype = covtype,
        myregF = myregF,
        mean.known,
        H.known = H.known,
        iso = iso,
        nugget = nugget
      )
      if (!is.null(env) && !is.na(loglik)) {
        if (is.null(env$min_loglik) || loglik > env$min_loglik) {
          env$min_loglik <- loglik
          env$arg_min <- par
        }
      }
      return(loglik)
    }
    if (gradient) {
      gr <- function(par, env) {
        gr_loglik <- gr_MLE(
          trend.type = trend.type,
          par = c(par, rep(1, d), 0, 0.5),
          X = X,
          t = t,
          Y = Y,
          l = l,
          covtype = covtype,
          myregF = myregF,
          mean.known = mean.known,
          H.known = H.known,
          iso = iso,
          nugget = nugget
        )
        return(gr_loglik)
      }
    } else gr <- NULL
    lower.params <- log(lower_phisq)
    upper.params <- log(upper_phisq)
  }
  if (is.null(H.known)) {
    if (!single_fidelity) {
      par <- c(log(c(myphi1sq, myphi2sq, mytausq)), myH)
      lower.par <- c(lower.params, eps)
      upper.par <- c(upper.params, 1 - eps)
    } else {
      par <- log(myphi1sq)
      lower.par <- lower.params
      upper.par <- upper.params
    }
  } else if (!is.null(H.known)) {
    par <- log(c(myphi1sq, myphi2sq, mytausq))
    lower.par <- lower.params
    upper.par <- upper.params
  }
  res <- try(optim(
    par = par,
    fn = fn,
    method = "L-BFGS-B",
    lower = lower.par,
    upper = upper.par,
    gr = gr,
    env = envtmp
  ))
  if (is(res, "try-error")) {
    res <- list(
      par = envtmp$arg_min,
      value = envtmp$min_loglik,
      counts = NA,
      message = "Optimization stopped due to NAs,
                  use best value so far"
    )
  }
  pars <- res$par
  if (!single_fidelity) {
    if (!iso) {
      phi1sq <- exp(pars[1:d])
      phi2sq <- exp(pars[(d + 1):(2 * d)])
      tausq <- exp(pars[2 * d + 1])
      if (is.null(H.known)) H <- pars[2 * d + 2] else H <- H.known
    } else {
      phi1sq <- exp(rep(pars[1], d))
      phi2sq <- exp(rep(pars[2], d))
      tausq <- exp(pars[3])
      if (is.null(H.known)) H <- pars[4] else H <- H.known
    }
  } else {
    if (!iso) {
      phi1sq <- exp(pars[1:d])
      phi2sq <- rep(0, d)
      tausq <- 0
      H <- 0
    } else {
      phi1sq <- exp(rep(pars[1], d))
      phi2sq <- rep(0, d)
      tausq <- 0
      H <- 0
    }
  }
  if (!is.null(H.known)) H <- H.known
  myK <- cov_gen(
    x1 = X,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = 1,
    sigma2sq = tausq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = nugget
  )
  Ki <- chol2inv(chol(myK))
  if (trend.type %in% c("OK", "UK")) {
    myKiF <- Ki %*% myregF
    myW <- crossprod(myregF, myKiF)
    myWi <- chol2inv(chol(myW))
    beta <- tcrossprod(myWi, myKiF) %*% Y
  }
  if (trend.type == "SK") tmp <- (Y - mean.known) else
    tmp <- Y - myregF %*% beta
  sigma1sq <- c(crossprod(tmp, Ki) %*% tmp) / (n - p)
  sigma2sq <- tausq * sigma1sq
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
  Ki <- chol2inv(chol(myK))
  if (is.null(H.known)) {
    if (trend.type == "SK") {
      estiP = list(
        phi1sq = phi1sq,
        phi2sq = phi2sq,
        sigma1sq = sigma1sq,
        sigma2sq = sigma2sq,
        H = H
      )
    } else if (trend.type %in% c("OK", "UK")) {
      estiP = list(
        phi1sq = phi1sq,
        phi2sq = phi2sq,
        sigma1sq = sigma1sq,
        sigma2sq = sigma2sq,
        beta = beta,
        H = H
      )
    }
  } else {
    if (trend.type == "SK") {
      estiP = list(
        phi1sq = phi1sq,
        phi2sq = phi2sq,
        sigma1sq = sigma1sq,
        sigma2sq = sigma2sq
      )
    } else if (trend.type %in% c("OK", "UK")) {
      estiP = list(
        phi1sq = phi1sq,
        phi2sq = phi2sq,
        sigma1sq = sigma1sq,
        sigma2sq = sigma2sq,
        beta = beta
      )
    }
  }
  MuFimodel <- list(
    X = X,
    t = t,
    Y = Y,
    estiP = estiP,
    regF = regF,
    Ki = Ki,
    used_args = list(
      covtype = covtype,
      trend.type = trend.type,
      trend.dim = trend.dim,
      trend.pol = trend.pol,
      interaction = interaction,
      H.known = H.known,
      param.bounds = param.bounds,
      mean.known = mean.known,
      l = l,
      iso = iso,
      nugget = nugget,
      Gijs = Gijs,
      ncores = ncores,
      gradient = gradient,
      init = init,
      single_fidelity = single_fidelity
    )
  )
  class(MuFimodel) <- "MuFiMeshGP"
  return(MuFimodel)
}
