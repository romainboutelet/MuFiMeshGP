#' IMSPE optimal design point search
#'
#' @description Search for the best next design point according to the IMSPE
#' criterion given a current \code{MuFiMeshGP} model fit.
#'
#' @param object Current \code{MuFiMeshGP} model fit.
#' @param t.min,t.max Lower and upper bounds on the fidelity space for the search.
#' @param cost.func Function that maps the tunable parameter \code{t} to the
#'  corresponding cost running a simulation at that fidelity level. For example,
#'  \code{function(t) 1/t^2}.
#' @param cost.new (optional) Cost of running a new simulation at a new fidelity
#'  level, scalar.
#' @param gr whether the gradient should be used in the optimization of
#'  the IMSPE. (Not recommended due to numerical errors)
#' @param gr_cost.func If \code{grad} is \code{TRUE}, the user needs to specify
#' the gradient of the cost function, as a function.
#' @param DesCand Design candidates to evaluate from.
#' @param Wijs,Hijs (optional) Matrices from previous IMSPE search to obtain
#' faster computation through matrix decomposition.
#' @param control list of arguments udes for the optimization.
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
#'

IMSPE_AL <- function(
  object,
  t.min,
  t.max,
  cost.func,
  cost.new = 0,
  gr = FALSE,
  gr_cost.func = NULL,
  DesCand = NULL,
  Wijs = NULL,
  Hijs = NULL,
  control = list(
    multi.start.n = 20,
    maxit = 20,
    DesStart = NULL,
    seed = NULL,
    ncores = 1
  )
) {
  if (!inherits(object, "MuFiMeshGP")) {
    stop("The object is not of class \"MuFiMeshGP\" \n")
  }
  eps <- .Machine$double.eps
  if (is.null(control$seed)) seed <- sample(1:2^15, 1)
  if (!is.null(control$seed)) ncores <- control$ncores else ncores <- 1
  if (is.null(control)) control <- list(multi.start.n = 20, maxit = 20)
  if (is.null(control$multi.start.n)) control$multi.start.n <- 20
  if (is.null(control$maxit)) control$maxit <- 20
  d <- ncol(object$X)
  if (is.null(Wijs))
    Wijs <- Wij(
      x1 = object$X,
      phi1sq = object$estiP$phi1sq,
      sigma1sq = object$estiP$sigma1sq,
      covtype = object$used_args$covtype
    )
  if (object$used_args$trend.type %in% c("UK", "OK")) {
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = object$X,
        phi1sq = object$estiP$phi1sq,
        sigma1sq = object$estiP$sigma1sq,
        covtype = object$used_args$covtype,
        trend.type = object$used_args$trend.type,
        trend.dim = object$used_args$trend.dim,
        trend.pol = object$used_args$trend.pol,
        interaction = object$used_args$interaction
      )
  }
  if (is.null(DesCand)) {
    if (!is.null(control$DesStart)) {
      DesStart <- control$DesStart
    } else {
      DesStart <- lhs::maximinLHS(control$multi.start.n, d + 1)
      # corners <- matrix(
      #   as.numeric(as.matrix(
      #     do.call(expand.grid, replicate(d + 1, c(0, 1), simplify = F))
      #   )),
      #   ncol = d + 1
      # )
      # DesStart <- rbind(DesStart, corners)
      DesStart[, d + 1] <- min(object$t) +
        DesStart[, d + 1] *
          (max(object$t) - min(object$t))
      if (t.min == t.max) {
        multi.n <- 10
        DesStart <- lhs::maximinLHS(multi.n, d)
      }
    }
    if (t.min != t.max) {
      fn <- function(par)
        return(crit_IMSPEcost.new(
          Des = par,
          object = object,
          cost.func = cost.func,
          cost.new = cost.new,
          Wijs = Wijs,
          Hijs = Hijs
        ))
      if (gr) {
        gr <- function(par) {
          return(grad_crit_IMSPEcost.new(
            Des = par,
            object = object,
            cost.func = cost.func,
            cost.new = cost.new,
            gr_cost.func = gr_cost.func,
            Wijs = Wijs,
            Hijs = Hijs
          ))
        }
      } else gr <- NULL
      local_opt_fun.new <- function(i) {
        out <- optim(
          DesStart[i, ],
          fn = fn,
          gr = gr,
          method = "L-BFGS-B",
          lower = c(rep(eps, d), t.min),
          upper = c(rep(1, d), t.max),
          control = list(maxit = control$maxit, fnscale = -1)
        )
        return(out)
      }
      all_res.new <- parallel::mclapply(
        1:nrow(DesStart),
        local_opt_fun.new,
        mc.cores = ncores
      )
      res_max.new <- which.max(Reduce(
        c,
        lapply(all_res.new, function(x) x$value)
      ))
      par.new <- drop(all_res.new[[res_max.new]]$par)
      res <- list(
        x = matrix(head(par.new, -1), ncol = d),
        t = tail(drop(all_res.new[[res_max.new]]$par), 1),
        value = all_res.new[[res_max.new]]$value,
        new = TRUE,
        id = NULL
      )
    } else if (t.min == t.max) {
      local_opt_fun.new <- function(i) {
        out <- optim(
          DesStart[i, ],
          crit_IMSPEcost.old,
          method = "L-BFGS-B",
          lower = rep(eps, d),
          upper = rep(1, d),
          t = t.max,
          Wijs = Wijs,
          Hijs = Hijs,
          object = object,
          cost.func = cost.func,
          control = list(maxit = control$maxit, fnscale = -1)
        )
        return(out)
      }
      all_res.new <- lapply(1:nrow(DesStart), local_opt_fun.new)
      res_max.new <- which.max(Reduce(
        c,
        lapply(all_res.new, function(x) x$value)
      ))
      par.new <- c(drop(all_res.new[[res_max.new]]$par), t.max)
      res <- list(
        x = matrix(drop(all_res.new[[res_max.new]]$par), ncol = d),
        t = t.max,
        value = all_res.new[[res_max.new]]$value,
        new = TRUE,
        id = NULL
      )
    }
    return(res)
  } else {
    # Needs to be adapted later (case when we use DesCand)
    crit_IMSPEcost_mcl <- function(i, object, Wijs, DesCand) {
      crit_IMSPEcost.new(
        Des = DesCand[i, ],
        object = object,
        cost.func = cost.func,
        cost.new = cost.new,
        Wijs = Wijs
      )
    }
    res <- unlist(parallel::mclapply(
      1:nrow(DesCand),
      crit_IMSPEcost_mcl,
      DesCand = DesCand,
      Wijs = Wijs,
      object = object,
      mc.cores = ncores
    ))
    tmp <- which(duplicated(
      rbind(
        cbind(object$X, object$t),
        DesCand[which.min(res), , drop = FALSE]
      ),
      fromLast = TRUE
    ))
    if (length(tmp) > 0)
      return(list(
        par = DesCand[which.min(res), , drop = FALSE],
        value = min(res),
        new = FALSE,
        id = tmp
      ))
    return(list(
      par = DesCand[which.min(res), , drop = FALSE],
      value = min(res),
      new = TRUE,
      id = NULL
    ))
  }
}
