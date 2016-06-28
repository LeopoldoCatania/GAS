#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec EvaluatePit_Univ(arma::mat mTheta, arma::vec vY, std::string Dist, int iT){

  arma::vec vU(iT);
  int t;

  for(t=0;t<iT;t++){
    vU(t) = pdist_univ(vY(t), mTheta.col(t), Dist);
  }

  return vU;

}

