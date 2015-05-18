#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]


double update_lambda(NumericVector sigma, int nu){
    int d = sigma.size();
    double shape = d*nu/2;
    double invscale;
    NumericVector sigmainv(d);
    for(int i = 0; i < d; i++){
        sigmainv(i) = 1/sigma(i); 
    }
    invscale = (nu*0.5)*sum(sigmainv);
    double lambda;
    lambda = R::rgamma(shape, 1/invscale);
    return lambda;
}
