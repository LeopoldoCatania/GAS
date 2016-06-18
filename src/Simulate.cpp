#include <RcppArmadillo.h>
#include "Utils.h"
#include "DistWrap.h"
#include "Mapping.h"
#include "GASInnovation.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List SimulateGAS_univ(int iT, arma::vec vKappa, arma::mat mA, arma::mat mB, std::string Dist, std::string ScalingType){

  int i,iK = NumberParameters(Dist);

  arma::vec vY(iT);
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;

  mTheta_tilde.col(0) = vKappa;
  mTheta.col(0) = MapParameters_univ(mTheta_tilde.col(0),Dist, iK);

  //initialise observables
  vY(0) = rdist_univ(mTheta.col(0),Dist);

  //initialise Likelihood
  vLLK(0) = ddist_univ(vY(0), mTheta.col(0), Dist, true);
  dLLK   += vLLK(0);

  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_univ(vY(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1), iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vIntercept + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_univ(mTheta_tilde.col(i),Dist, iK);
    if(i<iT){
      vY(i)   = rdist_univ(mTheta.col(i),Dist);
      vLLK(i) = ddist_univ(vY(i), mTheta.col(i), Dist, true);
      dLLK   += vLLK(i);
    }
  }

  List out;

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["vY"]   = vY;
  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return out;
}
//[[Rcpp::export]]
List SimulateGAS_multi(int iT, int iN, arma::vec vKappa, arma::mat mA, arma::mat mB, std::string Dist, std::string ScalingType){

  int i,iK = NumberParameters(Dist,iN);

  arma::mat mY(iN,iT);
  arma::vec vLLK(iT);
  double dLLK = 0;

  //initialise parameter
  arma::mat mTheta_tilde(iK,iT+1);
  arma::mat mTheta(iK,iT+1);
  arma::mat mInnovations(iK,iT);

  //initialise Dynamics
  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;

  mTheta_tilde.col(0) = vKappa;
  mTheta.col(0) = MapParameters_multi(mTheta_tilde.col(0),Dist, iN, iK);

  //initialise observables
  mY.col(0) = rdist_multi(mTheta.col(0),iN,Dist);

  //initialise Likelihood
  vLLK(0) = ddist_multi(mY.col(0), mTheta.col(0), iN,Dist, true);
  dLLK   += vLLK(0);
  // Dynamics
  for(i=1;i<iT+1;i++){
    mInnovations.col(i-1) = GASInnovation_multi(mY.col(i-1), mTheta.col(i-1), mTheta_tilde.col(i-1),iN, iK, Dist, ScalingType);
    mTheta_tilde.col(i)   = vIntercept + mA * mInnovations.col(i-1) + mB *  mTheta_tilde.col(i-1);
    mTheta.col(i)         = MapParameters_multi(mTheta_tilde.col(i),Dist, iN, iK);
    if(i<iT){
      mY.col(i) = rdist_multi(mTheta.col(i),iN,Dist);
      vLLK(i)   = ddist_multi(mY.col(i), mTheta.col(i), iN, Dist, true);
      dLLK     += vLLK(i);
    }
  }

  List out;

  out["mTheta"]       = mTheta;
  out["mInnovations"] = mInnovations;
  out["mTheta_tilde"] = mTheta_tilde;

  out["mY"]   = mY;
  out["vLLK"] = vLLK;
  out["dLLK"] = dLLK;

  return out;
}
