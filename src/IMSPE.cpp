#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double crossprod_cpp(arma::vec v1, arma::vec v2) {
  int n = v1.n_rows;
  double s;
  
  for (int i = 0; i < n; i++){
    s += v1(i)*v2(i);
  }
  
  return s;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec matMulVec_cpp(arma::mat M, arma::vec v) {
  int nr = M.n_rows;
  int nc = M.n_cols;
  arma::vec out(nr);
  
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nc; j++){ 
      out(i) += M(i,j)*v(j);
    }
  }
  
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec vecMulMat_cpp(arma::vec v, arma::mat M) {
  int nr = M.n_rows;
  int nc = M.n_cols;
  arma::vec out(nc);
  
  for (int i = 0; i < nc; i++){
    for (int j = 0; j < nr; j++){ 
      out(i) += v(j)*M(j,i);
    }
  }
  
  return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube sliced_matMulVec_cpp(arma::cube M, arma::vec v) {
  int nr = M.n_rows;
  int nc = M.n_cols;
  int ns = M.n_slices;
  
  arma::cube N(nr,nc,ns);
  
  for (int i = 0; i < ns; i++){
    N.slice(i) = matMulVec_cpp(M.slice(i),v);
  }
  
  return(N);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube sliced_vecMulMat_cpp(arma::cube v, arma::mat M) {
  int nr = v.n_rows;
  int nc = M.n_cols;
  int ns = v.n_slices;
  arma::mat tmp(nr,1);
  arma::cube N(nc,1,ns);
  
  for (int i = 0; i < ns; i++){
    tmp = v.slice(i);
    N.slice(i) = vecMulMat_cpp(tmp.col(0),M);
  }
  
  return(N);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mulMat_cpp(arma::mat M1, arma::mat M2) {
  int nr = M1.n_rows;
  int n = M1.n_cols;
  int nc = M2.n_cols;
  arma::mat M(nr,nc);
  
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nc; j++) {
      for (int k = 0; k < n; k++) {
        M(i,j) += M1(i,k) * M2(k,j);
      }
    }
  }
  return M;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube sliced_mulMat_cpp(arma::cube D, arma::mat M) {
  int nr = D.n_rows;
  int nc = D.n_cols;
  int ns = D.n_slices;
  
  arma::cube N(nr,nc,ns);
  
  for (int i = 0; i < ns; i++){
    N.slice(i) = mulMat_cpp(D.slice(i),M);
  }
  
  return N;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec hadam_vec_cpp(arma::vec v1, arma::vec v2){
  int n = v1.n_rows;
  
  arma::vec v(n);
  
  for (int i = 0; i < n; i++){
    v(i) = v1(i)*v2(i);
  }
  
  return(v);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat hadam_mat_cpp(arma::mat M1, arma::mat M2){
  int nr = M1.n_rows;
  int nc = M1.n_cols;
  
  arma::mat M(nr,nc);
  
  for (int i = 0; i < nc; i++){
    M.col(i) = hadam_vec_cpp(M1.col(i),M2.col(i));
  }
  
  return(M);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat dx_kn_gauss_cpp(arma::mat X, arma::vec x, arma::vec k1n,
                          arma::vec k2n, arma::vec phi1sq, arma::vec phi2sq) {
  int nr = X.n_rows;
  int ns = X.n_cols;
  arma::mat dx_kn(nr,ns);
  arma::vec tmp1(nr);
  arma::vec tmp2(nr);
  
  for (int i = 0; i < ns; i++){
    for (int j = 0; j < nr; j++){
      tmp1(j) = 2*phi1sq(i)*(X(j,i) - x(i));
      tmp2(j) = 2*phi2sq(i)*(X(j,i) - x(i));
    }
    dx_kn.col(i) = hadam_vec_cpp(tmp1,k1n) + hadam_vec_cpp(tmp2,k2n);
  }
  
  return dx_kn;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dt_kn_cpp(arma::vec t1, double t2, arma::vec k2n, double H,
                       double l) {
  int nr = t1.n_rows;
  arma::vec dt_kn(nr);
  double tmp;
  
  for (int i = 0; i < nr; i++){
    tmp = (pow(t2,2*H-1) + (t1(i)-t2) * pow(fabs(t1(i)-t2),2*H-2)) / 
      (pow(t1(i),2*H) + pow(t2,2*H) - pow(fabs(t1(i)-t2),2*H));
    dt_kn(i) = l*tmp*k2n(i);
  }
  
  return dt_kn;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double dt_k_cpp(double t, double sigma2sq, double l) {
  double dt_k;

  dt_k = sigma2sq*l*pow(t,l-1);
  
  return dt_k;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Wijs1_gauss_cpp(arma::mat x, arma::vec phi1sq, double sigma1sq)
{
  int nr = x.n_rows;
  int nc = x.n_cols;
  arma::mat Wijs(nr,nr);
  double tmp;
  
  for (int i = 0; i < nr; i++){
    for (int j = 0; j < nr; j++){
      tmp = 1;
      for (int k = 0; k < nc; k++){
        tmp *= sqrt(2*M_PI/phi1sq(k))/4*exp(-pow(x(i,k)-x(j,k),2)*phi1sq(k)/2) *
          (erf(sqrt(phi1sq(k)/2)*(2-x(i,k)-x(j,k))) +
          erf(sqrt(phi1sq(k)/2)*(x(i,k)+x(j,k))));
      }
      Wijs(i,j) = sigma1sq*sigma1sq*tmp;
    }
  }
  
  return Wijs;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat Wijs2_gauss_cpp(arma::mat x1, arma::mat x2, arma::vec phi1sq,
                          double sigma1sq) {
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int nc = x1.n_cols;
  arma::mat Wijs(n1,n2);
  double tmp;
  
  for (int i = 0; i < n1; i++){
    for (int j = 0; j < n2; j++){
      tmp = 1;
      for (int k = 0; k < nc; k++){
        tmp *= sqrt(2*M_PI/phi1sq(k))/4*exp(-pow(x1(i,k)-x2(j,k),2)*phi1sq(k)/2) *
          (erf(sqrt(phi1sq(k)/2)*(2-x1(i,k)-x2(j,k))) +
          erf(sqrt(phi1sq(k)/2)*(x1(i,k)+x2(j,k))));
      }
      Wijs(i,j) = sigma1sq*sigma1sq*tmp;
    }
  }
  
  return Wijs;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat dx_wn_gauss_cpp(arma::mat x1, arma::vec x2, arma::vec wn,
                          arma::vec phi1sq, double sigma1sq) {
  int n = x1.n_rows;
  int m = x1.n_cols;
  arma::mat dx_wn(n,m);
  double tmp;
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      if (m > 1){
        tmp = wn(i)/sqrt(2*M_PI/phi1sq(j))/4*exp(-pow(x1(i,j)-x2(j),2)*phi1sq(j)/2) *
          (erf(sqrt(phi1sq(j)/2)*(2-x1(i,j)-x2(j))) +
          erf(sqrt(phi1sq(j)/2)*(x1(i,j)+x2(j))));
      }
      else {
        tmp = sigma1sq*sigma1sq;
      }
      dx_wn(i,j) = tmp*1/2*exp(-phi1sq(j)*pow(x1(i,j)-x2(j),2)/2) * (
          exp(-phi1sq(j)*pow(x1(i,j)+x2(j),2)/2) - 
          exp(-phi1sq(j)*pow(x1(i,j)+x2(j)-2,2)/2) -
          sqrt(phi1sq(j)*M_PI/8) * (erf(sqrt(phi1sq(j)/2)*(2-x1(i,j)-x2(j))) +
          erf(sqrt(phi1sq(j)/2)*(x1(i,j)+x2(j))))
      );
    }
  }
  return dx_wn;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dx_w_gauss_cpp(arma::mat x, arma::vec phi1sq,double sigma1sq) {
  int m = x.n_cols;
  arma::vec dx_w(m);
  double tmp;
  
  for (int j = 0; j < m; j++){
    tmp = 1;
    for (int k = 0; k < j; k++){
      tmp *= sqrt(2*M_PI/phi1sq(k))/4 *
        (erf(sqrt(phi1sq(k)/2)*(2-2*x(0,k)))+erf(sqrt(phi1sq(k)/2)*(2*x(0,k))));
    }
    tmp *= exp(-2*phi1sq(j)*pow(x(0,j),2)) - exp(-2*phi1sq(j)*pow(1-x(0,j),2));
    for (int k = j+1; k < m; k++){
      tmp *= sqrt(2*M_PI/phi1sq(k))/4 *
        (erf(sqrt(phi1sq(k)/2)*(2-2*x(0,k)))+erf(sqrt(phi1sq(k)/2)*(2*x(0,k))));
    }
    dx_w(j) = sigma1sq*sigma1sq*tmp;
  }
  return dx_w;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sigmasq_cpp(arma::mat Ki, arma::vec kn, double k) {
  int n = Ki.n_cols;
  double sigmasq;
  arma::vec tmp(n);
  
  for (int i = 0; i < n; i++){
    sigmasq = k - crossprod_cpp(vecMulMat_cpp(kn,Ki),kn);
  }
  
  return sigmasq;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dx_sigmasq_cpp(arma::mat dx_kn, arma::mat Ki, arma::vec kn) {
  int n = dx_kn.n_cols;
  int m = dx_kn.n_rows;
  arma::vec dx_sigmasq(n);
  arma::vec tmp(m);
  
  for (int i = 0; i < n; i++){
    tmp = vecMulMat_cpp(dx_kn.col(i),Ki);
    dx_sigmasq(i) = -2*crossprod_cpp(tmp,kn);
  }
  
  return dx_sigmasq;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double dt_sigmasq_cpp(arma::vec dt_kn, arma::mat Ki, arma::vec kn,double dt_k) {
  int n = dt_kn.n_rows;
  
  double dt_sigmasq;
  arma::vec tmp(n);

  tmp = vecMulMat_cpp(dt_kn,Ki);
  dt_sigmasq = dt_k - 2*crossprod_cpp(tmp,kn);
  
  return dt_sigmasq;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat dx_a_cpp(arma::mat dx_kn, arma::mat Ki, arma::vec kn, double sigmasq,
                   arma::vec dx_sigmasq) {
  int n = Ki.n_cols;
  int m = dx_kn.n_cols;
  arma::mat dx_a(n,m);

  for (int i = 0; i < m; i++){
    dx_a.col(i) = dx_sigmasq(i)/pow(dx_sigmasq(i),2)*matMulVec_cpp(Ki,kn) - 
      matMulVec_cpp(Ki,dx_kn.col(i))/sigmasq;
  }
  
  return dx_a;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec dt_a_cpp(arma::vec dt_kn, arma::mat Ki, arma::vec kn, double sigmasq,
                double dt_sigmasq) {
  int n = dt_kn.n_rows;
  arma::vec dt_a(n);

  dt_a = dt_sigmasq/pow(dt_sigmasq,2)*matMulVec_cpp(Ki,kn) -
    matMulVec_cpp(Ki,dt_kn)/sigmasq;
  
  return dt_a;
}
