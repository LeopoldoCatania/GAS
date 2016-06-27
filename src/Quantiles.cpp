#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat Quantiles(arma::mat mTheta, std::string Dist, arma::vec vProbs){
  int iT = mTheta.n_cols;
  int iP = vProbs.size();

  arma::mat mQuantile(iP,iT);

  int t,p;
  for(t=0;t<iT;t++){
    for(p=0;p<iP;p++){
      mQuantile(p,t) = qdist_univ(vProbs(p), mTheta.col(t), Dist);
    }
  }
  return mQuantile.t();
}

