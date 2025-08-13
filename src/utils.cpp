#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector find_corres(NumericMatrix X0, NumericMatrix X) {
  
  // @param X0 matrix of unique designs
  // @param X matrix of all designs
  // @return vector associating rows of X with those of X0
  
  
  IntegerVector corres(X.nrow());
  
  bool tmp;
  int nX0r = X0.nrow();
  int nX0c = X0.ncol();
  
  for(int i = 0; i < X.nrow(); i++){
    for(int j = 0; j < nX0r; j++){
      tmp = true;
      for(int k = 0; k < nX0c; k++){
        if(X(i,k) != X0(j,k)){
          tmp = false;
          break;
        }
      }
      if(tmp){
        corres(i) = j + 1;
        break;
      }
    }
  }
  return(corres);
}


double c1i_gauss(double x1, double X, double sigma){
  double tmp = -1./2. * sqrt(M_PI/2.) * sigma * exp(-(X - x1) * (X - x1) / (2. * sigma * sigma)) * (erf((X + x1 - 2.) / (sqrt(2.) * sigma)) - erf((X + x1) / (sqrt(2.) * sigma)));
  if(tmp == 0.) return(0.);
  
  return((0.5*exp(-(x1 - X)*(x1 - X) / (2.*sigma*sigma))* (exp(-(x1 + X)*(x1 + X)/(2.*sigma*sigma)) -
         exp(-(2.-(x1 + X))*(2.-(x1 + X))/(2.*sigma * sigma))) - sqrt(2.*M_PI)/4./sigma*(x1 - X) * exp(-(x1 - X)*(x1 - X)/(2.*sigma*sigma)) * 
         (erf((x1 + X)/(sqrt(2.)*sigma)) - erf((x1 + X - 2.)/(sqrt(2.)*sigma))))/tmp);
}

// [[Rcpp::export]]
double c2_gauss_cpp(double x, double t, double w){
  if(w == 0.) return(0.);
  double tmp = -1./2. * sqrt(M_PI/2.) * t * (erf((2. * x  - 2.) / (sqrt(2.) * t)) - erf((2. * x) / (sqrt(2.) * t)));
  if(tmp == 0.) return(0.);
  return((exp(-2. * x * x / (t * t)) - exp(-2.*(1. - x) * (1. - x) / (t * t)))*w/tmp);
}


// [[Rcpp::export]]
NumericVector c1_gauss_cpp(NumericVector X, double x, double sigma, NumericVector W){
  NumericVector cis(X.length());  
  
  for(int i = 0; i < X.length(); i++){
    cis(i) = c1i_gauss(x, X(i), sigma) * W(i);
  }
  
  return(cis);
}