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

