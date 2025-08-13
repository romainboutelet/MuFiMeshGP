#' Calculates the IMSPE criterion
#'
#' @noRd

crit_IMSPEcost.new <- function(
  Des,
  object,
  cost.func,
  cost.new,
  Wijs = NULL,
  Hijs = NULL
) {
  eps <- sqrt(.Machine$double.eps)
  X <- object$X
  n <- nrow(X)
  d <- ncol(X)
  x <- matrix(Des[1:d], ncol = d)
  t <- tail(Des, 1)
  cost <- cost.func(t) + cost.new
  tt <- object$t
  covtype <- object$used_args$covtype
  iso <- object$used_args$iso
  phi1sq <- object$estiP$phi1sq
  phi2sq <- object$estiP$phi2sq
  sigma1sq <- object$estiP$sigma1sq
  sigma2sq <- object$estiP$sigma2sq
  trend.type <- object$used_args$trend.type
  trend.dim <- object$used_args$trend.dim
  trend.pol <- object$used_args$trend.pol
  interaction <- object$used_args$interaction
  nugget <- object$used_args$nugget
  l <- object$used_args$l
  if (is.null(object$used_args$H.known)) H <- object$estiP$H else
    H <- object$used_args$H.known
  Ki <- object$Ki
  if (is.null(Wijs))
    Wijs <- Wij(x1 = X, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  if (is.null(dim(x))) x <- matrix(x, ncol = d)
  wn <- Wij(
    x1 = X,
    x2 = x,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  )
  w <- Wij(x1 = x, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  k1n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  k2n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  kn <- k1n + k2n
  k <- cov_gen(
    x1 = x,
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
  sigmasq <- max(drop(k - crossprod(kn, crossprod(Ki, kn))), eps)
  a <- -crossprod(Ki, kn) / sigmasq
  R1 <- c(
    crossprod(a, crossprod(Wijs, a)) *
      sigmasq +
      2 * crossprod(wn, a) +
      w / sigmasq
  )
  if (trend.type == "SK") value <- R1 / cost else if (trend.type == "OK") {
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- matrix(1, nrow = nrow(X), ncol = 1)
    newF <- matrix(1, ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      nrow = 1
    )
    SKi <- colSums(Ki)
    v <- sqrt(sigmasq) * sum(a) + newF / sqrt(sigmasq)
    KiF <- matrix(SKi, ncol = 1)
    FKiF <- sum(Ki)
    FKiFi <- 1 / FKiF
    B <- (FKiFi * v^2 * FKiFi) /
      c(1 + v * FKiFi * v)
    S <- sigmasq * a * sum(a) + tcrossprod(a, newF)
    myT <- sum(a) + newF / sigmasq
    R2 <- (B *
      fast_trace_cpp(tcrossprod(KiF, KiF), Wijs) -
      (FKiFi - B) * fast_trace_cpp(tcrossprod(S + 2 * KiF, S), Wijs) -
      2 * myT * (FKiFi - B) %*% crossprod(S + KiF, wn) -
      myT^2 * (FKiFi - B) * w)
    R3 <- sum(diag(B))
    R4 <- -2 *
      (B *
        crossprod(KiF, Hijs) -
        (FKiFi - B) * crossprod(S, Hijs) -
        myT * (FKiFi - B) * newHijs)
    value <- c(R1 + R2 + R3 + R4) / cost
  } else if (trend.type == "UK") {
    if (is.null(object$used_args$Gijs)) {
      Gijs <- Gij(
        d = ncol(X),
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    } else Gijs <- object$used_args$Gijs
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- cbind(rep(1, nrow(X)), regF(X, tt))
    newF <- matrix(c(1, regF(x, t)), ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      ncol = 1
    )
    v <- sqrt(sigmasq) * t(myF) %*% a + newF / sqrt(sigmasq)
    KiF <- Ki %*% myF
    FKiF <- crossprod(myF, KiF)
    FKiFi <- chol2inv(chol(FKiF))
    B <- (FKiFi %*% tcrossprod(v, v) %*% FKiFi) /
      c(1 + crossprod(v, FKiFi) %*% v)
    S <- sigmasq * tcrossprod(a, a) %*% myF + tcrossprod(a, newF)
    myT <- crossprod(a, myF) + t(newF) / sigmasq
    R2 <- c(
      fast_trace_cpp(KiF %*% B, crossprod(KiF, Wijs)) -
        fast_trace_cpp((S + 2 * KiF) %*% tcrossprod(FKiFi - B, S), Wijs) -
        2 * crossprod((KiF + S) %*% tcrossprod((FKiFi - B), myT), wn) -
        myT %*% tcrossprod(FKiFi - B, myT) * w
    )
    R3 <- fast_trace_cpp(B, Gijs)
    R4 <- -2 *
      c(
        fast_trace_cpp(tcrossprod(B, KiF), Hijs) -
          fast_trace_cpp(tcrossprod(FKiFi - B, S), Hijs) -
          crossprod(tcrossprod(FKiFi - B, myT), newHijs)
      )
    value <- c(R1 + R2 + R3 + R4) / cost / sigma1sq
  }
  # bigF <- rbind(myF,c(newF))
  # bigK <- cov_gen(x1 = rbind(X,x), t1 = c(tt,t), phi1sq = phi1sq,
  #                 phi2sq = phi2sq, sigma1sq = sigma1sq,
  #                 sigma2sq = sigma2sq, l = l, covtype = covtype,
  #                 H = H, iso = iso, nugget = nugget)
  # bigKi <- chol2inv(chol(bigK))
  # bigKiF <- bigKi %*% bigF
  # bigFKiFi <- chol2inv(chol(crossprod(bigF, bigKi %*% bigF)))
  # bigWijs <- Wij(x1 = rbind(X,x), phi1sq = phi1sq, sigma1sq = sigma1sq,
  #                 covtype = covtype)
  # M <- tcrossprod(KiF %*% FKiFi, KiF)
  # bigHijs <- Hij(mu = rbind(X,x), phi1sq = phi1sq, sigma1sq = sigma1sq,
  #                covtype = covtype,
  #                trend.type = trend.type,
  #                trend.pol = trend.pol,
  #                interaction = interaction)
  # bigM <-  tcrossprod(bigKiF %*% bigFKiFi, bigKiF)
  # print(bigWijs - rbind(cbind(Wijs,wn),c(wn,w)))
  # print(bigKi - rbind(cbind(Ki + sigmasq*tcrossprod(a,a),a),c(a,1/sigmasq)))
  # print(FKiFi - B - chol2inv(chol(crossprod(bigF, bigKi %*% bigF))))
  # # print(Des)
  # print(fast_trace_cpp(M,Wijs))
  # print(fast_trace_cpp(FKiFi,Gijs))
  # print(-2*fast_trace_cpp(tcrossprod(FKiFi, KiF), Hijs))
  # print(fast_trace_cpp(bigM,bigWijs))
  # print(fast_trace_cpp(bigFKiFi,Gijs))
  # print(-2*fast_trace_cpp(tcrossprod(bigFKiFi, bigKiF), bigHijs))
  # print(R1 - fast_trace_cpp(bigKi,bigWijs) + fast_trace_cpp(Ki,Wijs))
  # print(R2 - fast_trace_cpp(M,Wijs) + fast_trace_cpp(bigM,bigWijs))
  # print(R3 - fast_trace_cpp(FKiFi,Gijs) + fast_trace_cpp(bigFKiFi,Gijs))
  # print(R4 - 2*fast_trace_cpp(tcrossprod(bigFKiFi, bigKiF), bigHijs) +
  #         2*fast_trace_cpp(tcrossprod(FKiFi, KiF), Hijs))
  # print(R1+R2+R3+R4)
  # print(c(R1, R2, R3, R4))
  return(value)
  # return(R1/cost)
}

#' Calculates the IMSPE criterion gradient
#'
#' @noRd

grad_crit_IMSPEcost.new <- function(
  Des,
  object,
  cost.func,
  cost.new,
  gr_cost.func,
  Wijs = NULL,
  Hijs = NULL
) {
  X <- object$X
  d <- ncol(X)
  x <- matrix(Des[1:d], ncol = d)
  t <- tail(Des, 1)
  cost <- cost.func(t) + cost.new
  tt <- object$t
  covtype <- object$used_args$covtype
  iso <- object$used_args$iso
  phi1sq <- object$estiP$phi1sq
  phi2sq <- object$estiP$phi2sq
  sigma1sq <- object$estiP$sigma1sq
  sigma2sq <- object$estiP$sigma2sq
  trend.type <- object$used_args$trend.type
  trend.dim <- object$used_args$trend.dim
  trend.pol <- object$used_args$trend.pol
  interaction <- object$used_args$interaction
  nugget <- object$used_args$nugget
  l <- object$used_args$l
  if (is.null(object$used_args$H.known)) H <- object$estiP$H else
    H <- object$used_args$H.known
  nugget <- nugget
  Ki <- object$Ki
  if (is.null(Wijs))
    Wijs <- Wij(x1 = X, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  if (is.null(dim(x))) x <- matrix(x, ncol = d)
  wn <- c(Wij(
    x1 = X,
    x2 = x,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  ))
  w <- c(Wij(x1 = x, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype))
  dx_wn <- dx_wn(
    x1 = X,
    x2 = x,
    wn = wn,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  )
  dx_w <- c(dx_w(
    x = x,
    w = w,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  ))
  k1n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  k2n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  kn <- k1n + k2n
  k <- cov_gen(
    x1 = x,
    t1 = t,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    sigma1sq = sigma1sq,
    sigma2sq = sigma2sq,
    l = l,
    covtype = covtype,
    H = H,
    iso = iso,
    nugget = 0
  )
  dx_kn <- dx_kn(
    X = X,
    x = x,
    k1n = k1n,
    k2n = k2n,
    phi1sq = phi1sq,
    phi2sq = phi2sq,
    covtype = covtype
  )
  dt_kn <- dt_kn_cpp(t1 = tt, t2 = t, k2n = k2n, H = H, l = l)
  dt_k <- dt_k_cpp(t = t, sigma2sq = sigma2sq, l = l)
  sigmasq <- max(drop(k - crossprod(kn, crossprod(Ki, kn))), nugget)
  dx_sigmasq <- c(dx_sigmasq_cpp(dx_kn = dx_kn, Ki = Ki, kn = kn))
  dt_sigmasq <- c(dt_sigmasq_cpp(dt_kn = dt_kn, Ki = Ki, kn = kn, dt_k = dt_k))
  a <- c(-crossprod(Ki, kn) / sigmasq)
  dx_a <- dx_a_cpp(
    dx_kn = dx_kn,
    Ki = Ki,
    kn = kn,
    sigmasq = sigmasq,
    dx_sigmasq = dx_sigmasq
  )
  dt_a <- c(dt_a_cpp(
    dt_kn = dt_kn,
    Ki = Ki,
    kn = kn,
    sigmasq = sigmasq,
    dt_sigmasq = dt_sigmasq
  ))
  R1 <- crossprod(a, crossprod(Wijs, a)) *
    sigmasq +
    2 * crossprod(wn, a) +
    w / sigmasq
  dx_R1 <- dx_sigmasq *
    c(crossprod(a, c(crossprod(Wijs, a)))) +
    2 * sigmasq * c(crossprod(c(crossprod(Wijs, a)), dx_a)) +
    2 * c(crossprod(dx_wn, a)) +
    2 * c(crossprod(wn, dx_a)) -
    dx_sigmasq / dx_sigmasq^2 * w +
    dx_w / sigmasq
  dt_R1 <- dt_sigmasq *
    c(crossprod(a, c(crossprod(Wijs, a)))) +
    2 * sigmasq * c(crossprod(c(crossprod(Wijs, a)), dt_a)) +
    2 * c(crossprod(wn, dt_a)) -
    dt_sigmasq / dt_sigmasq^2 * w
  dt_cost <- gr_cost.func(t)
  if (trend.type == "SK") {
    gradient <- c(dx_R1 / cost, (dt_R1 * cost - R1 * dt_cost) / cost^2)
  } else if (trend.type == "OK") {
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- matrix(1, nrow = nrow(X), ncol = 1)
    newF <- matrix(1, ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      nrow = 1
    )
    SKi <- colSums(Ki)
    v <- sqrt(sigmasq) * sum(a) + newF / sqrt(sigmasq)
    KiF <- matrix(SKi, ncol = 1)
    FKiF <- sum(Ki)
    FKiFi <- 1 / FKiF
    B <- (FKiFi * v^2 * FKiFi) /
      c(1 + v * FKiFi * v)
    S <- sigmasq * a * sum(a) + tcrossprod(a, newF)
    myT <- sum(a) + newF / sigmasq
    R2 <- -(B *
      fast_trace_cpp(tcrossprod(KiF, KiF), Wijs) -
      (FKiFi - B) * fast_trace_cpp(tcrossprod(S + 2 * KiF, S), Wijs) -
      2 * myT * (FKiFi - B) %*% crossprod(S + KiF, wn) -
      myT^2 * (FKiFi - B) * w)
    R3 <- -sum(diag(B))
    R4 <- -2 *
      (B *
        crossprod(KiF, Hijs) -
        (FKiFi - B) * crossprod(S, Hijs) -
        myT * (FKiFi - B) * newHijs)
    value <- c(R1 + R2 + R3 + R4) / cost
  } else if (trend.type == "UK") {
    if (is.null(object$used_args$Gijs)) {
      Gijs <- Gij(
        d = ncol(X),
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    } else Gijs <- object$used_args$Gijs
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- cbind(rep(1, nrow(X)), regF(X, tt))
    newF <- matrix(c(1, regF(x, t)), ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      ncol = 1
    )
    v <- sqrt(sigmasq) * t(myF) %*% a + newF / sqrt(sigmasq)
    KiF <- Ki %*% myF
    FKiF <- crossprod(myF, KiF)
    FKiFi <- chol2inv(chol(FKiF))
    B <- (FKiFi %*% tcrossprod(v, v) %*% FKiFi) /
      c(1 + crossprod(v, FKiFi) %*% v)
    S <- sigmasq * tcrossprod(a, a) %*% myF + tcrossprod(a, newF)
    myT <- crossprod(a, myF) + t(newF) / sigmasq
    R2 <- -(fast_trace_cpp(tcrossprod(KiF %*% B, KiF), Wijs) -
      fast_trace_cpp((S + 2 * KiF) %*% tcrossprod(FKiFi - B, S), Wijs) -
      2 * myT %*% tcrossprod(FKiFi - B, S + KiF) %*% wn -
      myT %*% tcrossprod(FKiFi - B, myT) * w)
    R3 <- -fast_trace_cpp(B, Gijs)
    R4 <- -2 *
      (fast_trace_cpp(tcrossprod(B, KiF), Hijs) -
        fast_trace_cpp(tcrossprod(FKiFi - B, S), Hijs) -
        c(myT %*% (FKiFi - B) %*% newHijs))
    value <- c(R1 + R2 + R3 + R4) / cost
  }
  # print(gradient)
  return(gradient)
}

#' Calculates the IMSPE criterion
#'
#' @noRd

crit_IMSPEcost.old <- function(
  x,
  t,
  object,
  cost.func,
  Wijs = NULL,
  Hijs = NULL
) {
  eps <- sqrt(.Machine$double.eps)
  X <- object$X
  tt <- object$t
  n <- nrow(X)
  d <- ncol(X)
  cost <- cost.func(t)
  covtype <- object$used_args$covtype
  iso <- object$used_args$iso
  phi1sq <- object$estiP$phi1sq
  phi2sq <- object$estiP$phi2sq
  sigma1sq <- object$estiP$sigma1sq
  sigma2sq <- object$estiP$sigma2sq
  trend.type <- object$used_args$trend.type
  trend.dim <- object$used_args$trend.dim
  trend.pol <- object$used_args$trend.pol
  interaction <- object$used_args$interaction
  nugget <- object$used_args$nugget
  l <- object$used_args$l
  if (is.null(object$used_args$H.known)) H <- object$estiP$H else
    H <- object$used_args$H.known
  Ki <- object$Ki
  if (is.null(Wijs))
    Wijs <- Wij(x1 = X, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  if (is.null(dim(x))) x <- matrix(x, ncol = d)
  wn <- Wij(
    x1 = X,
    x2 = x,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  )
  w <- Wij(x1 = x, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  k1n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  k2n <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  kn <- k1n + k2n
  k <- cov_gen(
    x1 = x,
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
  sigmasq <- max(drop(k - crossprod(kn, crossprod(Ki, kn))), eps)
  a <- -crossprod(Ki, kn) / sigmasq
  R1 <- c(
    crossprod(a, crossprod(Wijs, a)) *
      sigmasq +
      2 * crossprod(wn, a) +
      w / sigmasq
  )
  if (trend.type == "SK") value <- R1 / cost else if (trend.type == "OK") {
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- matrix(1, nrow = nrow(X), ncol = 1)
    newF <- matrix(1, ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      nrow = 1
    )
    SKi <- colSums(Ki)
    v <- sqrt(sigmasq) * sum(a) + newF / sqrt(sigmasq)
    KiF <- matrix(SKi, ncol = 1)
    FKiF <- sum(Ki)
    FKiFi <- 1 / FKiF
    B <- (FKiFi * v^2 * FKiFi) /
      c(1 + v * FKiFi * v)
    S <- sigmasq * a * sum(a) + tcrossprod(a, newF)
    myT <- sum(a) + newF / sigmasq
    R2 <- (B *
      fast_trace_cpp(tcrossprod(KiF, KiF), Wijs) -
      (FKiFi - B) * fast_trace_cpp(tcrossprod(S + 2 * KiF, S), Wijs) -
      2 * myT * (FKiFi - B) %*% crossprod(S + KiF, wn) -
      myT^2 * (FKiFi - B) * w)
    R3 <- sum(diag(B))
    R4 <- -2 *
      (B *
        crossprod(KiF, Hijs) -
        (FKiFi - B) * crossprod(S, Hijs) -
        myT * (FKiFi - B) * newHijs)
    value <- c(R1 + R2 + R3 + R4) / cost
  } else if (trend.type == "UK") {
    if (is.null(object$used_args$Gijs)) {
      Gijs <- Gij(
        d = ncol(X),
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    } else Gijs <- object$used_args$Gijs
    if (is.null(Hijs))
      Hijs <- Hij(
        mu = X,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      )
    regF <- object$regF
    myF <- cbind(rep(1, nrow(X)), regF(X, tt))
    newF <- matrix(c(1, regF(x, t)), ncol = 1)
    newHijs <- matrix(
      Hij(
        mu = x,
        phi1sq = phi1sq,
        sigma1sq = sigma1sq,
        covtype = covtype,
        trend.type = trend.type,
        trend.dim = trend.dim,
        trend.pol = trend.pol,
        interaction = interaction
      ),
      ncol = 1
    )
    v <- sqrt(sigmasq) * t(myF) %*% a + newF / sqrt(sigmasq)
    KiF <- Ki %*% myF
    FKiF <- crossprod(myF, KiF)
    FKiFi <- chol2inv(chol(FKiF))
    B <- (FKiFi %*% tcrossprod(v, v) %*% FKiFi) /
      c(1 + crossprod(v, FKiFi) %*% v)
    S <- sigmasq * tcrossprod(a, a) %*% myF + tcrossprod(a, newF)
    myT <- crossprod(a, myF) + t(newF) / sigmasq
    R2 <- c(
      fast_trace_cpp(KiF %*% B, crossprod(KiF, Wijs)) -
        fast_trace_cpp((S + 2 * KiF) %*% tcrossprod(FKiFi - B, S), Wijs) -
        2 * crossprod((KiF + S) %*% tcrossprod((FKiFi - B), myT), wn) -
        myT %*% tcrossprod(FKiFi - B, myT) * w
    )
    R3 <- fast_trace_cpp(B, Gijs)
    R4 <- -2 *
      c(
        fast_trace_cpp(tcrossprod(B, KiF), Hijs) -
          fast_trace_cpp(tcrossprod(FKiFi - B, S), Hijs) -
          crossprod(tcrossprod(FKiFi - B, myT), newHijs)
      )
    value <- c(R1 + R2 + R3 + R4) / cost / sigma1sq
  }
  # bigF <- rbind(myF,c(newF))
  # bigK <- cov_gen(x1 = rbind(X,x), t1 = c(tt,t), phi1sq = phi1sq,
  #                 phi2sq = phi2sq, sigma1sq = sigma1sq,
  #                 sigma2sq = sigma2sq, l = l, covtype = covtype,
  #                 H = H, iso = iso, nugget = nugget)
  # bigKi <- chol2inv(chol(bigK))
  # bigKiF <- bigKi %*% bigF
  # bigFKiFi <- chol2inv(chol(crossprod(bigF, bigKi %*% bigF)))
  # bigWijs <- Wij(x1 = rbind(X,x), phi1sq = phi1sq, sigma1sq = sigma1sq,
  #                 covtype = covtype)
  # M <- tcrossprod(KiF %*% FKiFi, KiF)
  # bigHijs <- Hij(mu = rbind(X,x), phi1sq = phi1sq, sigma1sq = sigma1sq,
  #                covtype = covtype,
  #                trend.type = trend.type,
  #                trend.pol = trend.pol,
  #                interaction = interaction)
  # bigM <-  tcrossprod(bigKiF %*% bigFKiFi, bigKiF)
  # print(bigWijs - rbind(cbind(Wijs,wn),c(wn,w)))
  # print(bigKi - rbind(cbind(Ki + sigmasq*tcrossprod(a,a),a),c(a,1/sigmasq)))
  # print(FKiFi - B - chol2inv(chol(crossprod(bigF, bigKi %*% bigF))))
  # # print(Des)
  # print(fast_trace_cpp(M,Wijs))
  # print(fast_trace_cpp(FKiFi,Gijs))
  # print(-2*fast_trace_cpp(tcrossprod(FKiFi, KiF), Hijs))
  # print(fast_trace_cpp(bigM,bigWijs))
  # print(fast_trace_cpp(bigFKiFi,Gijs))
  # print(-2*fast_trace_cpp(tcrossprod(bigFKiFi, bigKiF), bigHijs))
  # print(R1 - fast_trace_cpp(bigKi,bigWijs) + fast_trace_cpp(Ki,Wijs))
  # print(R2 - fast_trace_cpp(M,Wijs) + fast_trace_cpp(bigM,bigWijs))
  # print(R3 - fast_trace_cpp(FKiFi,Gijs) + fast_trace_cpp(bigFKiFi,Gijs))
  # print(R4 - 2*fast_trace_cpp(tcrossprod(bigFKiFi, bigKiF), bigHijs) +
  #         2*fast_trace_cpp(tcrossprod(FKiFi, KiF), Hijs))
  # print(R1+R2+R3+R4)
  # print(c(R1, R2, R3, R4))

  return(value)
}

#' Calculates the IMSPE criterion gradient
#'
#' @noRd

deriv_IMSPEcost.old <- function(
  x,
  t,
  object,
  cost.func,
  deriv_cost.func,
  Wijs = NULL,
  out = 2
) {
  X <- object$X
  tt <- object$t
  phi1sq <- object$estiP$phi1sq
  phi2sq <- object$estiP$phi2sq
  sigma1sq <- object$estiP$sigma1sq
  sigma2sq <- object$estiP$sigma2sq
  covtype <- object$used_args$covtype
  iso <- object$used_args$iso
  trend.type <- object$used_args$trend.type
  trend.dim <- object$used_args$trend.dim
  trend.pol <- object$used_args$trend.pol
  interaction <- object$used_args$interaction
  nugget <- object$used_args$nugget
  l <- object$used_args$l
  if (is.null(object$used_args$H.known)) H <- object$estiP$H else
    H <- object$used_args$H.known
  Ki <- object$Ki
  wn <- Wij(
    x1 = X,
    x2 = x,
    phi1sq = phi1sq,
    sigma1sq = sigma1sq,
    covtype = covtype
  )
  w <- Wij(x1 = x, phi1sq = phi1sq, sigma1sq = sigma1sq, covtype = covtype)
  kn1 <- cov_gen(
    x1 = X,
    x2 = x,
    t1 = tt,
    t2 = t,
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
  sigma_sq <- max(
    drop(sigma1sq^2 + sigma2sq^2 * t^l - crossprod(kn1, crossprod(Ki, kn1))),
    .Machine$double.eps
  )
}
