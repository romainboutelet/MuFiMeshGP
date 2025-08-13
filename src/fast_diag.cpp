#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fast_diag_cpp(NumericMatrix A, NumericMatrix B){
  int nr = A.nrow();
  int nc = A.ncol();
  NumericVector di(nr);
  double mprod;
  
  for(int i = 0; i < nr; i++){
    mprod = 0;
    for(int j = 0; j < nc; j++){
      mprod += A(i, j) * B(j, i);
    }
    di(i) = mprod;
  }
  
  return(di);
}
