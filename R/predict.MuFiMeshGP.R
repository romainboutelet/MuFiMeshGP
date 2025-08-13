#' Prediction of the MuFiMeshGP emulator for any fidelity level.
#'
#' @description The function computes the posterior mean and standard deviation of the
#' \code{MuFiMeshGP} model.
#'
#' @seealso \code{\link{MuFiMeshGP}} for the model
#'
#' @details From the object fitted by \code{\link{MuFiMeshGP}} or \code{\link{update.MuFiMeshGP}}
#'  the posterior mean and standard deviation are calculated for any input
#'  location and fidelity level.
#'  For details, see Boutelet and Sung (2025, <arXiv:2503.23158>).
#'
#' @param object an object of class \code{MuFiMeshGP}.
#' @param x matrix of new input locations to predict.
#' @param t vector or new fidelity levels to use for predictions.
#' @param ... no other argument.
#'
#' @return
#' \itemize{
#'   \item \code{mean}: vector of predictive posterior mean.
#'   \item \code{sd}: vector of predictive posterior standard deviation.
#' }
#'
#' @rdname predict.MuFiMeshGP
#' @title predict.MuFiMeshGP
#' @method predict MuFiMeshGP
#' @export
#' @examples
#' # Example code
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

predict.MuFiMeshGP <- function(object, x, t, ...) {
  if (!inherits(object, "MuFiMeshGP")) {
    stop("The object is not of class \"MuFiMeshGP\" \n")
  }
  X <- object$X
  Y <- object$Y
  tt <- object$t
  phi1sq <- object$estiP$phi1sq
  phi2sq <- object$estiP$phi2sq
  sigma1sq <- object$estiP$sigma1sq
  sigma2sq <- object$estiP$sigma2sq
  l <- object$l
  nugget <- object$used_args$nugget
  if (is.null(object$used_args$H.known)) H <- object$estiP$H else
    H <- object$used_args$H.known
  if (object$used_args$trend.type == "SK")
    mean.known <- object$used_args$mean.known
  nrowX <- nrow(X)
  d <- ncol(X)
  Ki <- object$Ki
  kn <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = object$used_args$l,
    covtype = object$used_args$covtype,
    H = H,
    iso = object$used_args$iso,
    nugget = nugget
  )
  k <- diag(cov_gen(
    x1 = x,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = object$used_args$l,
    covtype = object$used_args$covtype,
    H = H,
    iso = object$used_args$iso,
    nugget = 0
  ))
  if (object$used_args$trend.type == "SK") {
    myCalm11 <- Y - mean.known
    mean <- mean.known + crossprod(kn, crossprod(Ki, myCalm11))
    sd <- sqrt(ifelse(
      k - colSums(kn * crossprod(Ki, kn)) >= 0,
      k - colSums(kn * crossprod(Ki, kn)),
      0
    ))
  }
  if (object$used_args$trend.type == "OK") {
    beta <- c(object$estiP$beta)
    myregF <- matrix(1, nrow = nrowX, ncol = 1)
    myregf <- matrix(rep(1, nrow(x)), nrow = 1)
    myCalm11 <- Y - myregF %*% beta
    mean <- crossprod(myregf, beta) + crossprod(kn, crossprod(Ki, myCalm11))
    Var1 <- ifelse(
      k - colSums(kn * crossprod(Ki, kn)) >= 0,
      k - colSums(kn * crossprod(Ki, kn)),
      0
    )
    gamma <- myregf - crossprod(myregF, crossprod(Ki, kn))
    M <- crossprod(myregF, crossprod(Ki, myregF))
    Var2 <- ifelse(
      colSums(gamma * solve(M, gamma)) >= 0,
      colSums(gamma * solve(M, gamma)),
      0
    )
    sd <- sqrt(Var1 + Var2)
  }
  if (object$used_args$trend.type == "UK") {
    beta <- object$estiP$beta
    myregF <- cbind(rep(1, nrowX), object$regF(X, tt))
    myregf <- t(cbind(rep(1, nrow(x)), object$regF(x, t)))
    myCalm11 <- Y - myregF %*% beta
    mean <- crossprod(myregf, beta) + crossprod(kn, crossprod(Ki, myCalm11))
    Var1 <- ifelse(
      k - colSums(kn * crossprod(Ki, kn)) >= 0,
      k - colSums(kn * crossprod(Ki, kn)),
      0
    )
    gamma <- myregf - crossprod(myregF, crossprod(Ki, kn))
    M <- crossprod(myregF, crossprod(Ki, myregF))
    Var2 <- ifelse(
      colSums(gamma * solve(M, gamma)) >= 0,
      colSums(gamma * solve(M, gamma)),
      0
    )
    sd <- sqrt(Var1 + Var2)
  }
  return(list(mean = mean, sd = sd))
}
