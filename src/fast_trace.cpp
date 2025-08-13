#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fast_trace_cpp(NumericMatrix A, NumericMatrix B) {
  double tmp = 0;
  int nrA = A.nrow();
  int ncA = A.ncol();
  const double* ptrA = (const double*) &A(0,0);
  const double* ptrB = (const double*) &B(0,0);
  for(int i = 0; i < nrA; i++){
    ptrA = &A(i,0);
    for(int j = 0; j < ncA; j++, ptrA+=nrA, ptrB++){
      tmp += *ptrA * *ptrB;
    }
  }
  return(tmp);
}
