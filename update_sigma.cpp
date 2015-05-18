// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector update_sigma(arma::mat Theta, arma::mat Y, arma::mat M, arma::vec lambda, int nu){
  int n = Y.n_cols;
  int K = Y.n_rows;
  int R = Theta.n_cols;
  double N = n;
  arma::vec Den(n);
  arma::mat Thetasq(n, R);
  arma::mat Residual(K, n);
  arma::mat Res_sq(K, n);
  arma::mat Tmp(n, K);
  //rowvec scale;
  //rowvec sigma(K);
  arma::rowvec scale;
  arma::rowvec sigma(K);
  
  for (int j = 0; j < R; j++){
    Thetasq.col(j) = pow(Theta.col(j), 2);
  }
  for (int i = 0; i < n; i++){
    Den(i) = sum(Thetasq.row(i));
  }
  Residual = Y - M*arma::trans(Theta);
  for(int i = 0; i < K; i++){
    Res_sq.row(i) = pow(Residual.row(i),2); 
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < K; j++){
      Tmp(i,j) = Res_sq(j,i)/Den(i);
    }  
  }
  scale = (nu*lambda)/(N + nu) + (1/(N + nu))*(sum(Tmp));
  for (int k = 0; k < K; k++){
    sigma(k) = ((N + nu)*scale(k))/(R::rchisq(N + nu));
  }

  return (wrap(sigma));
}