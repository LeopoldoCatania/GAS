#include <RcppArmadillo.h>
#include "Mapping.h"
#include "IMWrap.h"
#include "Utils.h"
#include "ScoreWrap.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec GASInnovation_univ(double dY, arma::vec vTheta, arma::vec vTheta_tilde, int iK, std::string Dist, std::string ScalingType){

  arma::vec vS_tilde(iK);

  arma::vec vScore       = Score_univ(dY,vTheta,Dist);
  arma::mat mJ           = MapParametersJacobian_univ(vTheta_tilde, Dist, iK);
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
  vS_tilde = Thresholding_vec(vS_tilde, 1e5);
  vS_tilde = NaN2Zero(vS_tilde);
  vS_tilde = InfRemover_vec(vS_tilde, 1e5);

  return vS_tilde;
}

//[[Rcpp::export]]
arma::vec GASInnovation_multi(arma::vec vY, arma::vec vTheta, arma::vec vTheta_tilde, int iN, int iK,  std::string Dist, std::string ScalingType){

  arma::vec vScore   = Score_multi(vY, vTheta, iN, Dist);
  arma::mat mJ       = MapParametersJacobian_multi(vTheta_tilde, Dist, iN, iK);
  arma::vec vS_tilde = mJ.t() * vScore ;

  vS_tilde = Thresholding_vec(vS_tilde, 1e5);
  vS_tilde = NaN2Zero(vS_tilde);
  vS_tilde = InfRemover_vec(vS_tilde, 1e5);

  return vS_tilde;
}
