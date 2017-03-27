#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double StaticLLK_Univ(arma::vec vY, arma::vec vTheta, int iT, std::string Dist){
  double dLLK = 0.0;
  for(int i=0;i<iT;i++){
    dLLK += ddist_univ(vY(i), vTheta, Dist, true);
  }
  return dLLK;
}

//[[Rcpp::export]]
double StaticLLK_Multi(arma::mat mY, arma::vec vTheta, int iT, int iN, std::string Dist){
  double dLLK = 0.0;
  for(int i=0;i<iT;i++){
    dLLK += ddist_multi(mY.col(i), vTheta, iN, Dist, true);
  }
  return dLLK;
}
