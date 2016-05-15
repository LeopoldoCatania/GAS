#include <RcppArmadillo.h>
#include "GASInnovation.h"
#include "DistWrap.h"
#include "Mapping.h"

using namespace Rcpp;
using namespace arma;

List GASFilter_univ(arma::vec vY, arma::vec vKappa, arma::vec mA, arma::vec mB, int iT, int iK, std::string Dist, std::string ScalingType){

  int i;
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;
  mTheta.col(0) = vKappa;

  //initialise Likelihood
  vLLK(0) = ddist_univ(vY(0), mTheta.col(0), Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_univ(vY(i), mTheta.col(i-1), mTheta_tilde.col(i-1), iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vIntercept + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters(mTheta_tilde.col(i),Dist, iK);
    if(i<iT){
      vLLK(i) = ddist_univ(vY(i), mTheta.col(i), Dist, true);
      dLLK   += vLLK(i);
    }
  }

  List out;

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return out;
}
