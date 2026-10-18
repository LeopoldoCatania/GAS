#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat Quantiles(arma::mat mTheta, std::string Dist, arma::vec vProbs) {
  if (vProbs.n_elem == 0 || !vProbs.is_finite() ||
      arma::any(vProbs < 0) || arma::any(vProbs > 1))
    Rcpp::stop("probs must be a nonempty vector of finite probabilities in [0,1].");
  arma::mat out(mTheta.n_cols,vProbs.n_elem);
  for (arma::uword t=0;t<mTheta.n_cols;++t)
    for (arma::uword p=0;p<vProbs.n_elem;++p)
      out(t,p)=qdist_univ(vProbs(p),mTheta.col(t),Dist);
  return out;
}

