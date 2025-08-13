#' Updates the \code{MuFiMeshGP} model fit with new observations
#'
#' @description The function updates the current \code{MuFiMeshGP} model.
#'
#' @seealso \code{\link{MuFiMeshGP}} for initializing the model.
#'
#' @details From the model fitted by \code{\link{MuFiMeshGP}} or \code{\link{update.MuFiMeshGP}}
#'  the posterior mean and standard deviation are calculated for any input
#'  location and fidelity level.
#'  For details, see Boutelet and Sung (2025, <arXiv:2503.23158>).
#'
#' @param object an object of class \code{MuFiMeshGP}.
#' @param x matrix of new input locations.
#' @param t new tunable parameter, a scalar.
#' @param y observation corresponding to input location \code{x} and tunable
#'  parameter \code{t}.
#' @param param.estim if \code{TRUE}, the hyper-parameters are estimated by running it
#'  through \code{\link{MuFiMeshGP}}. If \code{FALSE}, the hyper-parameters from
#'  \code{object} are used to update the \code{MuFiMeshGP} model fit.
#' @param init See \code{\link{MuFiMeshGP}}.
#' @param ... no other argument.
#'
#' @return a list which is given the S3 class "MuFiMeshGP"
#' @rdname update.MuFiMeshGP
#' @title update.MuFiMeshGP
#' @method update MuFiMeshGP
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
#' par(mfrow = c(2,1))
#' plot(xx, ftrue, "l", ylim = c(-1,1.3), ylab = "y", xlab = "x")
#' lines(c(xx), mu, col = "blue")
#' lines(c(xx), lower, col = "blue", lty = 2)
#' lines(c(xx), upper, col = "blue", lty = 2)
#' points(c(X), Y, col = "red")
#'
#' ### RMSE ###
#' print(sqrt(mean((ftrue - mu))^2))
#'
#' best <- IMSPE_AL(fit.mufimeshgp, 0.5, 2, function(t) return(1 / t^2))
#' new.Y <- f(best$x, best$t)
#' fit.mufimeshgp <- update(fit.mufimeshgp, best$x, best$t, new.Y)
#'
#' pred.mufimeshgp <- predict(fit.mufimeshgp, xx, rep(0, 101))
#' mu <- pred.mufimeshgp$mean
#' s <- pred.mufimeshgp$sd
#' lower <- mu + qnorm(0.025)*s
#' upper <- mu + qnorm(0.975)*s
#'
#' plot(xx, ftrue, "l", ylim = c(-1,1.3), ylab = "y", xlab = "x")
#' lines(c(xx), mu, col = "blue")
#' lines(c(xx), lower, col = "blue", lty = 2)
#' lines(c(xx), upper, col = "blue", lty = 2)
#' points(c(X), Y, col = "red")
#' points(c(best$x), new.Y, col = "green")
#'
#' ### RMSE ###
#' print(sqrt(mean((ftrue - mu))^2))

update.MuFiMeshGP <- function(
  object,
  x,
  t,
  y,
  param.estim = TRUE,
  init = NULL,
  ...
) {
  if (!inherits(object, "MuFiMeshGP")) {
    stop("The object is not of class \"MuFiMeshGP\" \n")
  }
  d <- ncol(object$X)
  x <- matrix(x, ncol = d)
  X <- rbind(object$X, x)
  t <- c(object$t, t)
  Y <- c(object$Y, y)
  regF <- object$regF
  estiP <- object$estiP
  covtype <- object$used_args$covtype
  trend.type <- object$used_args$trend.type
  trend.dim <- object$used_args$trend.dim
  trend.pol <- object$used_args$trend.pol
  interaction <- object$used_args$interaction
  H.known <- object$used_args$H.known
  param.known <- object$used_args$param.known
  param.bounds <- object$used_args$param.bounds
  mean.known <- object$used_args$mean.known
  l <- object$used_args$l
  iso <- object$used_args$iso
  nugget <- object$used_args$nugget
  ncores <- object$used_args$ncores
  gradient <- object$used_args$gradient
  single_fidelity <- object$used_args$single_fidelity
  if (is.null(H.known)) H <- object$estiP$H else H <- H.known
  if (any(is.na(y))) {
    estiP <- object$estiP
    Ki <- chol2inv(chol(cov_gen(
      x1 = X,
      t1 = t,
      phi1sq = estiP$phi1sq,
      phi2sq = estiP$phi2sq,
      sigma1sq = estiP$sigma1sq,
      sigma2sq = estiP$sigma2sq,
      l = l,
      covtype = covtype,
      H = H,
      iso = iso,
      nugget = nugget
    )))
    used_args <- object$used_args
    return(list(
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
        param.known = param.known,
        param.bounds = param.bounds,
        mean.known = mean.known,
        l = l,
        iso = iso,
        nugget = nugget,
        ncores = ncores,
        gradient = gradient,
        init = init,
        single_fidelity = single_fidelity
      )
    ))
  }
  if (param.estim) {
    return(MuFiMeshGP(
      X = X,
      t = t,
      Y = Y,
      covtype = covtype,
      trend.type = trend.type,
      trend.dim = trend.dim,
      trend.pol = trend.pol,
      interaction = interaction,
      H.known = H.known,
      param.known = param.known,
      param.bounds = param.bounds,
      mean.known = mean.known,
      l = l,
      iso = iso,
      nugget = nugget,
      ncores = ncores,
      gradient = gradient,
      init = init,
      single_fidelity = single_fidelity
    ))
  } else {
    estiP <- object$estiP
    kn <- cov_gen(
      x1 = object$X,
      x2 = x,
      t1 = object$t,
      t2 = t,
      phi1sq = estiP$phi1sq,
      phi2sq = estiP$phi2sq,
      sigma1sq = estiP$sigma1sq,
      sigma2sq = estiP$sigma2sq,
      l = l,
      covtype = covtype,
      H = H,
      iso = iso,
      nugget = nugget
    )
    k <- cov_gen(
      x1 = x,
      t1 = t,
      phi1sq = estiP$phi1sq,
      phi2sq = estiP$phi2sq,
      sigma1sq = estiP$sigma1sq,
      sigma2sq = estiP$sigma2sq,
      l = l,
      covtype = covtype,
      H = H,
      iso = iso,
      nugget = nugget
    )
    sigma2 <- c(k - crossprod(kn, object$Ki %*% kn))
    myc <- -object$Ki %*% kn / sigma2
    Ki <- rbind(
      cbind(object$Ki + sigma2 * tcrossprod(myc, myc), myc),
      cbind(t(myc), 1 / sigma2)
    )
    updated <- list(
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
        param.known = param.known,
        param.bounds = param.bounds,
        mean.known = mean.known,
        l = l,
        iso = iso,
        nugget = nugget,
        ncores = ncores,
        gradient = gradient,
        init = init,
        single_fidelity = single_fidelity
      )
    )
    class(updated) <- "MuFiMeshGP"
    return(updated)
  }
}
