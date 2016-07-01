#include <RcppArmadillo.h>
#include "DistWrap.h"
#include "Mapping.h"
#include "GASInnovation.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List uGASMultiForcast(arma::vec vTheta_tp1, arma::vec vKappa, arma::mat mA, arma::mat mB,
                      int iH, int iB, int iK, std::string Dist, std::string ScalingType, bool bReturnsDraws){

  arma::cube cTheta = zeros(iK,iH,iB);
  arma::mat  mY     = zeros(iB,iH);

  int h,b;
  for(b=0;b<iB;b++){
    cTheta.slice(b).col(0) = vTheta_tp1;
    mY(b,0)                = rdist_univ(vTheta_tp1, Dist);
  }

  arma::vec vInnovation(iK);
  arma::vec vTheta_tilde_t(iK);
  arma::vec vTheta_tilde_tp1(iK);

  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;

  for(h=1;h<iH;h++){
     for(b=0;b<iB;b++){
       vTheta_tilde_t         = UnmapParameters_univ(cTheta.slice(b).col(h-1), Dist, iK);
       vInnovation            = GASInnovation_univ(mY(b,h-1), cTheta.slice(b).col(h-1),
                                                   vTheta_tilde_t, iK, Dist, ScalingType);
       vTheta_tilde_tp1       = vIntercept + mA * vInnovation + mB * vTheta_tilde_t;
       cTheta.slice(b).col(h) = MapParameters_univ(vTheta_tilde_tp1, Dist, iK);
       mY(b,h)                = rdist_univ(cTheta.slice(b).col(h), Dist);
    }
  }

  List lOut;

  if(bReturnsDraws) lOut["mY"]     = mY;
  lOut["cTheta"] = cTheta;

  return lOut;

}

//[[Rcpp::export]]
List mGASMultiForcast(arma::vec vTheta_tp1, arma::vec vKappa, arma::mat mA, arma::mat mB,
                      int iH, int iB, int iK, int iN, std::string Dist, std::string ScalingType,
                      bool bReturnsDraws){

  arma::cube cTheta = zeros(iK,iH,iB);
  arma::cube cY     = zeros(iN,iH,iB);

  int h,b;
  for(b=0;b<iB;b++){
    cTheta.slice(b).col(0) = vTheta_tp1;
    cY.slice(b).col(0)     = rdist_multi(vTheta_tp1, iN, Dist);
  }

  arma::vec vInnovation(iK);
  arma::vec vTheta_tilde_t(iK);
  arma::vec vTheta_tilde_tp1(iK);

  arma::vec vIntercept = ( eye(iK,iK) - mB) * vKappa;

  for(h=1;h<iH;h++){
    for(b=0;b<iB;b++){
      vTheta_tilde_t         = UnmapParameters_multi(cTheta.slice(b).col(h-1), Dist, iN, iK);
      vInnovation            = GASInnovation_multi(cY.slice(b).col(h-1), cTheta.slice(b).col(h-1),
                                                  vTheta_tilde_t, iN, iK, Dist, ScalingType);
      vTheta_tilde_tp1       = vIntercept + mA * vInnovation + mB * vTheta_tilde_t;
      cTheta.slice(b).col(h) = MapParameters_multi(vTheta_tilde_tp1, Dist, iN, iK);
      cY.slice(b).col(h)     = rdist_multi(cTheta.slice(b).col(h), iN, Dist);
    }
  }

  List lOut;

  if(bReturnsDraws) lOut["cY"] = cY;
  lOut["cTheta"] = cTheta;

  return lOut;

}
