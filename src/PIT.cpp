#include <RcppArmadillo.h>
#include "DistWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec EvaluatePit_Univ(arma::mat mTheta, arma::vec vY, std::string Dist, int iT) {
  if (iT < 0 || vY.n_elem < static_cast<arma::uword>(iT) ||
      mTheta.n_cols < static_cast<arma::uword>(iT))
    Rcpp::stop("Incompatible PIT data and parameter dimensions.");
  const bool discrete=Dist=="ber" || Dist=="poi" || Dist=="negbin" || Dist=="skellam";
  arma::vec u(iT);
  for (int t=0;t<iT;++t) {
    const double hi=pdist_univ(vY(t),mTheta.col(t),Dist);
    if (discrete) {
      if (!R_FINITE(vY(t)) || vY(t) != floor(vY(t)))
        Rcpp::stop("Discrete PIT requires finite integer observations.");
      const double lo=pdist_univ(vY(t)-1.0,mTheta.col(t),Dist);
      if (!R_FINITE(lo) || !R_FINITE(hi) || lo > hi || lo < 0 || hi > 1)
        Rcpp::stop("Invalid CDF values in randomized PIT.");
      u(t)=lo+Rf_runif(0,1)*(hi-lo);
    } else u(t)=hi;
  }
  return u;
}

