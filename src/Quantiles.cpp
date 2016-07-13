#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat Quantiles(arma::mat mTheta, std::string Dist, arma::vec vProbs){
  int iT = mTheta.n_cols;
  int iP = vProbs.size();

  int t,p;

  arma::mat mQuantile;

  if(iP == iT){ // guess you have a time series

    mQuantile = zeros(1,iT);

    for(t=0;t<iT;t++){
        mQuantile(0,t) = qdist_univ(vProbs(t), mTheta.col(t), Dist);
    }

  }else{

    mQuantile = zeros(iP,iT);

    for(t=0;t<iT;t++){
      for(p=0;p<iP;p++){
        mQuantile(p,t) = qdist_univ(vProbs(p), mTheta.col(t), Dist);
      }
    }
  }
  return mQuantile.t();
}

