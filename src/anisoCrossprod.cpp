#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat anisoCrossprod_cpp(arma::cube D, arma::vec phisq) {
  int nr = D.n_rows;
  int nc = D.n_cols;
  int ns = D.n_slices;
  arma::mat M(nr,nc);
  
  for (int i = 0; i < ns; i++){
    M += D.slice(i)*phisq(i);
  }
  
  return M;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube anisoCrossprod2_cpp(arma::cube D, arma::vec phisq) {
  int nr = D.n_rows;
  int nc = D.n_cols;
  int ns = D.n_slices;
  arma::cube M(nr,nc,ns);
  
  for (int i = 0; i < ns; i++){
    M.slice(i) = D.slice(i)*phisq(i);
  }
  
  return M;
}
