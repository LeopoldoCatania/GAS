#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec EvaluateLogScore_Univ(arma::mat mTheta, arma::vec vY, std::string Dist, int iT){

  arma::vec vLS(iT);
  int t;

  for(t=0;t<iT;t++){
    vLS(t) = ddist_univ(vY(t), mTheta.col(t), Dist, true);
  }

  return vLS;

}

//[[Rcpp::export]]
arma::vec EvaluateLogScore_Multi(arma::mat mTheta, arma::mat mY, std::string Dist, int iT, int iN){

  arma::vec vLS(iT);
  int t;

  for(t=0;t<iT;t++){
    vLS(t) = ddist_multi(mY.col(t), mTheta.col(t), iN, Dist, true);
  }

  return vLS;

}
