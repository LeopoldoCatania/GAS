#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::mat EvalMoments_univ(arma::mat mTheta, std::string Dist){
  int iT = mTheta.n_cols;

  arma::mat mMoments(4,iT);

  int t;
  for(t=0;t<iT;t++){
    mMoments.col(t) = mdist_univ(mTheta.col(t), Dist);
  }
  return mMoments.t();
}

//[[Rcpp::export]]
List EvalMoments_multi(arma::mat mTheta, std::string Dist, int iN){
  int iT = mTheta.n_cols;

  arma::mat  mMu(iN,iT);
  arma::cube mCov(iN,iN,iT);

  int t;
  for(t=0;t<iT;t++){
    mMu.col(t)    = mdist_multi_mean(mTheta.col(t), Dist, iN);
    mCov.slice(t) = mdist_multi_cov(mTheta.col(t), Dist, iN);
  }

  List lMoments;

  lMoments["mean"] = mMu.t();
  lMoments["cov"]  = mCov;

  return lMoments;
}
