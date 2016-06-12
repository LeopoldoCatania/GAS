#include <RcppArmadillo.h>
#include "Mapping.h"
#include "IMWrap.h"
#include "ScoreWrap.h"

using namespace Rcpp;
using namespace arma;

arma::vec GASInnovation_univ(double dY, arma::vec vTheta, arma::vec vTheta_tilde, int iK, std::string Dist, std::string ScalingType){

  arma::vec vS_tilde(iK);

  arma::vec vScore       = Score_univ(dY,vTheta,Dist);
  arma::mat mJ           = MapParametersJacobian(vTheta_tilde, Dist, iK);
  arma::vec Score_tilde  = mJ.t() * vScore ;

  if(ScalingType=="Identity") {

    vS_tilde = Score_tilde;
  }
  if(ScalingType=="Inv") {

    arma::mat mIM = IM_univ(vTheta, Dist);

    arma::mat mIM_Tilde_inv = mJ.t() * pinv(mIM) * mJ;

    vS_tilde = mIM_Tilde_inv * Score_tilde;

  }
  if(ScalingType=="InvSqrt") {

    arma::mat mIM= IM_univ(vTheta,Dist);

    arma::mat K = chol(mIM);

    vS_tilde = pinv(K) * Score_tilde;

  }
  //
  return vS_tilde;
}