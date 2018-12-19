// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace Rcpp;     

// [[Rcpp::export]]
double vSB(arma::vec x, double p){
  double n = x.n_elem;
  double ct = 1;
  double m = arma::mean(x);
  arma::vec z = x - m;
  double v = arma::sum(z % z) / n;

  for (int i = 1; i < n; i++) {
    arma::vec aux1 = z.subvec(i, n-1);
    arma::vec aux2 = z.subvec(0, n-i-1);
    double aux3 = arma::sum(aux1 % aux2)/n;
    double aux4 = (1-ct/n)*pow(1-p,i) + (ct/n)*pow(1-p,n-i);
    v += 2*aux4*aux3;
    ct += 1;
  }
  
  return(v);
  
}

// [[Rcpp::export]]
arma::vec vSB_mat(arma::mat x, double p){
  int n = x.n_cols;
  arma::vec out = arma::zeros(n);
  for (int i = 1; i < (n+1); i++) {
    out(i-1) = vSB(x.col(i-1), p);
  }
  return(out);
}

// [[Rcpp::export]]
double vHACC(arma::vec v, int k){
  int nv = v.n_rows;
  arma::mat vcv = (v.t() * v) / nv;
  
  if (k > 0){
    arma::vec w = arma::ones(k) - arma::linspace<arma::vec>(1, k, k)/(k+1);
    for (int i = 1; i < (k+1); i++) {
      arma::uvec aux1 = arma::linspace<arma::uvec>(i, nv-1, nv-i);
      arma::uvec aux2 = arma::linspace<arma::uvec>(0, nv-i-1, nv-i);
      arma::mat cov = arma::trans(v.elem(aux1)) * v.elem(aux2) / nv;
      vcv += w(i-1)*(cov + cov.t());
    }
  }
  
  double out = as_scalar(vcv);
  return(out);
  
}

// [[Rcpp::export]]
double t_testC(arma::vec x, int k, double eps){
  double m = arma::mean(x);
  double s = sqrt(vHACC(x-m, k));
  if (s > eps){
    double out = (sqrt(x.n_rows)*m)/s;
    return(out);
  } else {
    double out_na = NA_REAL;
    return(out_na);
  }
}

// [[Rcpp::export]]
arma::vec t_matC(arma::mat x, int k, double eps){
  int n = x.n_cols;
  arma::vec out = arma::zeros(n);
  for (int i = 1; i < (n+1); i++) {
    out(i-1) = t_testC(x.col(i-1), k, eps);
  }
  return(out);
}

// [[Rcpp::export]]
arma::vec vHACC_mat(arma::mat x, int k){
  int n = x.n_cols;
  arma::vec out = arma::zeros(n);
  for (int i = 1; i < (n+1); i++) {
    out(i-1) = vHACC(x.col(i-1), k);
  }
  return(out);
}
