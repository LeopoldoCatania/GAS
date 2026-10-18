#include "MomentWarnings.h"
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
    mMoments.col(t) = mdist_univ_quiet(mTheta.col(t), Dist);
  }
  GAS_warn_univ_moments(mMoments, Dist);
  return mMoments.t();
}

//[[Rcpp::export]]
List EvalMoments_multi(arma::mat mTheta, std::string Dist, int iN){
  int iT = mTheta.n_cols;

  arma::mat  mMu(iN,iT);
  arma::cube mCov(iN,iN,iT);

  int t;
  for(t=0;t<iT;t++){
    mMu.col(t)    = mdist_multi_mean_quiet(mTheta.col(t), Dist, iN);
    mCov.slice(t) = mdist_multi_cov_quiet(mTheta.col(t), Dist, iN);
  }

  std::string missing;
  if (!mMu.is_finite()) missing = "mean";
  if (!mCov.is_finite()) missing += missing.empty() ? "covariance" : ", covariance";
  GAS_warn_moments(Dist, missing);

  List lMoments;

  lMoments["mean"] = mMu.t();
  lMoments["cov"]  = mCov;

  return lMoments;
}
