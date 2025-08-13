#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube anisoDist1_cpp(arma::mat x) {
  int nr = x.n_rows;
  int ns = x.n_cols;
  arma::cube D(nr,nr,ns);
  
  for(int k = 0; k < ns; k++){
    for(int i = 0; i < nr; i++){
      for(int j = 0; j < i; j++){
        D(i,j,k) = pow(x(i,k)-x(j,k),2);
        D(j,i,k) = pow(x(i,k)-x(j,k),2);
      }
    }
  }
  return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat isoDist1_cpp(arma::mat x) {
  int nr = x.n_rows;
  int ns = x.n_cols;
  arma::mat D(nr,nr);
  
  for(int k = 0; k < ns; k++){
    for(int i = 0; i < nr; i++){
      for(int j = 0; j < i; j++){
        D(i,j) += pow(x(i,k)-x(j,k),2);
        D(j,i) += pow(x(i,k)-x(j,k),2);
      }
    }
  }
  return D;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube anisoDist2_cpp(arma::mat x1, arma::mat x2) {
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int ns = x1.n_cols;
  arma::cube D(n1,n2,ns);
  
  for(int k = 0; k < ns; k++){
    for(int i = 0; i < n1; i++){
      for(int j = 0; j < n2; j++){
        D(i,j,k) = pow(x1(i,k)-x2(j,k),2);
      }
    }
  }
  return D;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat isoDist2_cpp(arma::mat x1, arma::mat x2) {
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int ns = x1.n_cols;
  arma::mat D(n1,n2);
  
  for(int k = 0; k < ns; k++){
    for(int i = 0; i < n1; i++){
      for(int j = 0; j < n2; j++){
        D(i,j) += pow(x1(i,k)-x2(j,k),2);
      }
    }
  }
  return D;
}