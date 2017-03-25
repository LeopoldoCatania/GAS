#include <RcppArmadillo.h>
#include "GASInnovation.h"
#include "DistWrap.h"
#include "Mapping.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List GASFilter_univ(arma::vec vY, arma::vec vKappa, arma::mat mA, arma::mat mB, int iT, int iK, std::string Dist, std::string ScalingType){

  int i;
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  mTheta_tilde.col(0) = arma::inv( eye(iK,iK) - mB) * vKappa;;

  mTheta.col(0) = MapParameters_univ(mTheta_tilde.col(0),Dist, iK);

  //initialise Likelihood
  vLLK(0) = ddist_univ(vY(0), mTheta.col(0), Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_univ(vY(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1), iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vKappa + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_univ(mTheta_tilde.col(i),Dist, iK);
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

//[[Rcpp::export]]
List GASFilter_multi(arma::mat mY, arma::vec vKappa, arma::mat mA, arma::mat mB, int iT, int iN, int iK, std::string Dist, std::string ScalingType){

  int i;
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  mTheta_tilde.col(0) = arma::inv( eye(iK,iK) - mB) * vKappa;

  mTheta.col(0) = MapParameters_multi(mTheta_tilde.col(0),Dist, iN, iK);

  //initialise Likelihood
  vLLK(0) = ddist_multi(mY.col(0), mTheta.col(0),iN, Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_multi(mY.col(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1),iN, iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vKappa + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_multi(mTheta_tilde.col(i),Dist,iN, iK);
    if(i<iT){
      vLLK(i) = ddist_multi(mY.col(i), mTheta.col(i), iN,Dist, true);
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
