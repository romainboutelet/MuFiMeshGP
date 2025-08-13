#include <Rcpp.h>
using namespace Rcpp;

/// minMat functions

// [[Rcpp::export]]
NumericMatrix minMat1_cpp(NumericVector x) {
  int n = x.length();
  NumericMatrix tmat(n,n);
  
  for(int i = 0; i < n; i++){
    tmat(i,i) = x(i);
    for (int j = 0; j < i; j++){
      if(x(i)<x(j)){
        tmat(i,j) = x(i);
        tmat(j,i) = x(i);
      }else{
        tmat(i,j) = x(j);
        tmat(j,i) = x(j);
      }
    }
  }
  return tmat;
}

// [[Rcpp::export]]
NumericMatrix minMat2_cpp(NumericVector x1, NumericVector x2) {
  int n1 = x1.length();
  int n2 = x2.length();
  NumericMatrix tmat(n1,n2);
  
  for(int i = 0; i < n1; i++){
    for (int j = 0; j < n2; j++){
      if(x1(i)<x2(j)){
        tmat(i,j) = x1(i);
      }else{
        tmat(i,j) = x2(j);
      }
    }
  }
  return tmat;
}

// [[Rcpp::export]]
NumericMatrix FBMmat1_cpp(NumericVector x, double H, double l) {
  int n = x.length();
  NumericMatrix tmat(n,n);
  for(int i = 0; i < n; i++){
    tmat(i,i) = pow(x(i),l);
    for (int j = 0; j < i; j++){
      if (x(i) == x(j)){
        tmat(i,j) = pow(x(i),l);
      } else{
        tmat(i,j) = pow((pow(x(i),2*H) + pow(x(j),2*H)-
          pow(fabs(x(i)-x(j)),2*H))/2, l/(2*H));
      }
      tmat(j,i) = tmat(i,j);
    }
  }
  return tmat;
}

// [[Rcpp::export]]
NumericMatrix FBMmat2_cpp(NumericVector x1, NumericVector x2, double H,
                          double l) {
  int n1 = x1.length();
  int n2 = x2.length();
  NumericMatrix tmat(n1,n2);
  for(int i = 0; i < n1; i++){
    for (int j = 0; j < n2; j++){
      if (x1(i) == x2(j)){
        tmat(i,j) = pow(x1(i),l);
      } else{
        tmat(i,j) = pow((pow(x1(i),2*H) + pow(x2(j),2*H) -
          pow(fabs(x1(i)-x2(j)),2*H))/2, l/(2*H));
      }
    }
  }
  return tmat;
}

// [[Rcpp::export]]
NumericMatrix dK_H_cpp(NumericVector x, double H, double l) {
  int n = x.length();
  NumericMatrix dK_H(n,n);
  for(int i = 0; i < n; i++){
    dK_H(i,i) = 0;
    for (int j = 0; j < i; j++){
      if(x(i) == x(j)){
        dK_H(i,j) = 0;
      }else{
        dK_H(i,j) = l/H * ( log(x(i))*pow(x(i),2*H) + log(x(j))*pow(x(j),2*H) - 
            log(fabs(x(i)-x(j)))*pow(fabs(x(i)-x(j)),2*H) ) / 
          (pow(x(i),2*H) + pow(x(j),2*H) - pow(fabs(x(i)-x(j)),2*H)) - 
          l/(2*pow(H,2)) * ( log((pow(x(i), 2*H) + pow(x(j), 2*H) - 
          pow(fabs(x(i)-x(j)), 2*H)) / 2) );
        dK_H(j,i) = dK_H(i,j);
      }
    }
  }
  return dK_H;
}
