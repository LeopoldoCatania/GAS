#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::mat EvalMoments(arma::mat mTheta, std::string Dist){
  int iT = mTheta.n_cols;

  arma::mat mMoments(4,iT);

  int t;
  for(t=0;t<iT;t++){
    mMoments.col(t) = mdist_univ(mTheta.col(t), Dist);
  }
  return mMoments.t();
}



