#include <RcppArmadillo.h>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;
//using namespace std;

/* Random shuffle function to use instead of sample() */
inline int randWrapper(const int n) { return floor(unif_rand()*n); }
Rcpp::IntegerVector randomShuffle(Rcpp::IntegerVector a) {
    /* clone a into b to leave a alone */
    Rcpp::IntegerVector b = Rcpp::clone(a);
    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}

/*Function to square and sum each element of a vector*/
double sumSquared(arma::vec myvec){
  int nelem = myvec.size();
  double sumsq = 0;
  vec squared(nelem);
  
  for (int i = 0; i < nelem; i++){
    squared(i) = pow(myvec(i), 2);
  }
  sumsq = sum(squared);
  return sumsq;
}


// [[Rcpp::export]]

List update_ThetaStar(vec ThetaStar, vec y, mat M, vec sigma, vec Tx, vec sig, int iter){
  /* Variable declaration */
  /*Integers:   K        = number of markers
                R        = number of plants
                j        = index of current value
                nestim   = number of proportions being estimated */
                
  int K             = M.n_rows;
  int R             = M.n_cols;
  int j             = 0;
  int nestim        = 0;
  
  /* Double precision:  sumThetak    = sum of proportions not being estimated
                        theta        = proposed value of Theta
                        thetaNew     = new value of Theta given that it fulfils requirements
                        ThetaSq      = Sum(theta_ij^2) for current value
                        newThetaSq   = Sum(theta_ij^2) for new value
                        Diff         = Difference in log posteriors for proposed vs current value of theta */
                        
  double sumThetak  = 0;
  double theta      = 0;
  double thetaNew   = 0;
  double ThetaSq    = 0;
  double newThetaSq = 0;
  double Diff       = 0;
  
  /* Integer vectors:   plantj   = sequence 1 to R
                        tmpj     = shuffled elements of plantj
                        lastj    = proportion not to be estimated 
                        inTheta  = proportions to be estimated
                        tmpe     = temporal to assign values
                        ignore   = index containing elements of theta not to be estimated 
                        estim    = unsigned vector index of elements of theta to be estimated */  
                        
  IntegerVector plantj; 
  IntegerVector tmpj;
  IntegerVector lastj;
  IntegerVector inTheta;
  IntegerVector tmpe;
  IntegerVector ignore(2);
  uvec estim(R - 2);
  
  /* Double Vectors:    rho          = 1 or 0 if proportion was accepted
                        newThetaStar = vector of new values of Theta
                        Thetak       = vector of values not estimated of Theta
                        noise        = 1-element vector, draw from Normal(0, 1)
                        muStar       = mu vector for current value of Theta
                        muStarNew    = mu vector for new value of Theta
                        Res          = residual for current value
                        newRes       = residual for new value
                        alpha        = 1-element vector, draw from Uniform(0, 1) */
                        
  vec rho          = zeros(R);
  vec newThetaStar = zeros(R);
  vec Thetak;
  vec noise(1);
  vec muStar;
  vec muStarNew;
  vec Res;
  vec newRes;
  vec alpha(1);
  
  /* List: Out = Contains 3 elements:   ThetaStar
                                        rho
                                        sig */
  List Out;
  
  /* Kernel of the sampling function */
  plantj    = seq_len(R);
  tmpj      = randomShuffle(plantj);
  lastj     = tmpj(0);
  inTheta   = Rcpp::setdiff(tmpj, lastj);
  
  for(IntegerVector::iterator jit = inTheta.begin(); jit != inTheta.end(); jit++){
    j = (int)*jit - 1;

    if(R == 2){
      Thetak = zeros(1);
    }
    else{
      nestim = R - 2;
      ignore(0) = (int)*jit;
      ignore(1) = lastj(0);
      tmpe = Rcpp::setdiff(plantj, ignore);
      
      for(int i = 0; i < nestim; i++){
        estim(i) = tmpe(i) - 1; 
      }      
      Thetak = ThetaStar.elem(estim);
    }
    sumThetak = sum(Thetak);
    
    /* Adjusting variance for proposal*/
    if((Tx(j) > 0.4) && ((iter - 1)%100) == 0){
      sig(j) = sig(j)*5;
    }
    else if((Tx(j) < 0.3) && ((iter - 1)%100) == 0){
      sig(j) = sig(j)/5;
    }
    
    /* proposed value for theta = old value plus random noise from normal(0, 1) */
    noise = rnorm(1);
    theta = ThetaStar(j) + std::sqrt(sig(j))*noise(0);
    
    if(theta > 0 && theta < (1 - sumThetak)){
      thetaNew = theta;
      
      /* (re)Assembling ThetaStar into a new vector newThetaStar 
      in the 2 plant case the one of the thetas is 1 - thetaNew */
      if (R == 2){
       newThetaStar(j) = thetaNew;
       newThetaStar(lastj(0) - 1) = 1 - thetaNew;
      }
      else{
        newThetaStar(j) = thetaNew;
        newThetaStar(lastj(0) - 1) = 1 - sumThetak - thetaNew; 
        newThetaStar.elem(estim) = ThetaStar.elem(estim);
      }  
      
      /* Squaring and summing each element of ThetaStar and newThetaStar */ 
      ThetaSq = sumSquared(ThetaStar);
      newThetaSq = sumSquared(newThetaStar);
      
      /* Calculating mu for current and new value of ThetaStar */
      muStar = M*ThetaStar;
      muStarNew = M*newThetaStar;
      
      /* Residuals for current and new value of Theta */
      Res = y - muStar;
      newRes = y - muStarNew;
      
      /* log of the ratio between new value and current value of Theta*/
      Diff = as_scalar(0.5*((trans(newRes/sigma)*newRes)/newThetaSq - (trans(Res/sigma)*Res)/ThetaSq) + (K/2)*log(newThetaSq/ThetaSq));
      }
    else{
      Diff = 2.5e25;
      }
    /* Accept / Reject */
    alpha = runif(1);
    if(Diff < 0.0 || exp(-Diff) > alpha(0)){
      ThetaStar = newThetaStar;
      rho(j) = 1;
    }
    else{
      rho(j) = 0;
    }
  }
  
  Out["ThetaStar"] = trans(ThetaStar);
  Out["rho"] = rho; 
  Out["sig"] = sig;
  return (wrap(Out));
}