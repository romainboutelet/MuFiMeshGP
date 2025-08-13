#' Generates the covariance matrix
#'
#' @description Generates the covariance matrix for the Gaussian kernel or
#'  Matern kernel
#'
#' @param x1,x2 Design input location matrices. \code{x2} is to be used only to
#'  create cross-covariance matrix.
#' @param t1,t2 Design tunable parameter vectors. \code{t2} is to be used only
#' to create cross-covariance matrix.
#' @param phi1sq,phi2sq,sigma1sq,sigma2sq,H,l hyper-parameters for the
#'  covariance function
#' @param covtype kernel function used: \code{"Gaussian"} is the only available
#'  one at the moment.
#' @param iso If \code{TRUE}, then the covariance function is isotropic.
#' If \code{FALSE}, the covariance function is anisotropic.
#' @param nugget optional
#' @export

cov_gen <- function(
  x1,
  x2 = NULL,
  t1,
  t2 = NULL,
  phi1sq,
  phi2sq,
  sigma1sq,
  sigma2sq,
  H,
  l = 4,
  covtype,
  iso,
  nugget = sqrt(.Machine$double.eps)
) {
  d <- ncol(x1)
  if (iso & length(phi1sq) < d) {
    phi1sq <- rep(phi1sq, d)
    phi2sq <- rep(phi2sq, d)
  }
  if (covtype == "Gaussian") {
    if (!is.null(x2) && is.null(t2)) {
      print("if x2 is not NULL then t2 must be given")
      return(NULL)
    } else if (is.null(x2)) {
      D <- distance.MuFiMeshGP(x1 = x1)
      tmat <- FBMmat1_cpp(t1, H, l)
    } else {
      D <- distance.MuFiMeshGP(x1 = x1, x2 = x2)
      tmat <- FBMmat2_cpp(t1, t2, H, l)
    }
    D1 <- anisoCrossprod_cpp(D, phi1sq)
    D2 <- anisoCrossprod_cpp(D, phi2sq)
    if (is.null(x2)) {
      return(
        sigma1sq *
          (exp(-D1) + diag(nugget, nrow(D1))) +
          sigma2sq * tmat * (exp(-D2) + diag(nugget, nrow(D1)))
      )
    } else {
      if (all(dim(x1) == dim(x2)) && all(x1 == x2))
        print("x2 should not be used in cov_gen")
      return(sigma1sq * exp(-D1) + sigma2sq * exp(-D2) * tmat)
    }
  }
}
