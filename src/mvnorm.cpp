#include <RcppArmadillo.h>
#include "Utils.h"
#include "SafeFun.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat rmvnorm_mat(int iN, arma::vec vMu, arma::mat mSigma) {
  int incols = mSigma.n_cols;
  arma::mat mY = arma::randn(iN, incols);
  return arma::repmat(vMu, 1, iN).t() + mY * chol_safe(mSigma);
}
arma::mat rmvnorm_ThetaParam(arma::vec vTheta,int iN, int iJ) { //iJ = # of draws

  int iK = NumberParameters("mvnorm", iN);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iK-1);

  arma::mat mD = diagmat(vSigma);
  arma::mat mR = build_mR(vRho, iN);

  arma::mat mSigma = mD * mR * mD;

  arma::mat mY = rmvnorm_mat(iJ, vMu, mSigma);

  return mY;

}

const double log2pi = std::log(2.0 * M_PI);
double dmvnorm(arma::vec vY,
               arma::vec vMu,
               arma::mat mSigma,
               bool bLog = false) {
  int iN = vY.size();
  double dOut=0;
  arma::mat mrooti = arma::trans(arma::inv(trimatu(chol_safe(mSigma))));
  double drootisum = arma::sum(log(mrooti.diag()));
  double dconstants = -(static_cast<double>(iN)/2.0) * log2pi;

  arma::vec vZ = mrooti *  (vY - vMu) ;
  dOut = dconstants - 0.5 * arma::sum(vZ%vZ) + drootisum;

  if (bLog == false) {
    dOut = exp(dOut);
  }
  return dOut;
}
double dmvnorm_ThetaParam(arma::vec vY,
                          arma::vec vTheta,
                          int iN,
                          bool bLog = false) {
  int iL = 2*iN + iN*(iN-1)/2;

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iL-1);

  arma::mat mD = diagmat(vSigma);
  arma::mat mR = build_mR(vRho, iN);

  arma::mat mSigma = mD * mR * mD;

  double dPDF = dmvnorm(vY, vMu, mSigma, bLog);

  return dPDF;

}
arma::vec RhoScore(arma::vec vR, arma::mat mD, arma::vec vY, arma::vec vMu, int iN){

  arma::mat mRho_S = zeros(iN,iN);
  arma::mat mU     = zeros(iN,iN);

  arma::mat mR = build_mR(vR,iN);

  int i, j;
  arma::vec vV   = mD.i() * (vY - vMu);
  arma::mat mR_i = mR.i();

  for(i = 0;i<iN;i++){
    for(j = 0;j<iN;j++){
      if(i!=j){
        mU(i,j) = 1.0;
        mRho_S(i,j) = as_scalar(vV.t() * mR_i * mU * mR_i * vV);
        mU(i,j) = 0.0;
      }
    }
  }

  mRho_S = mRho_S - mR_i.t();

  arma::vec vR_s = build_vR(mRho_S, iN);

  return vR_s;

}
arma::vec DScore(arma::mat mD, arma::mat mR, arma::vec vY, arma::vec vMu, int iN){

  arma::mat mU = zeros(iN,iN);

  arma::vec vD_s(iN);

  int i;
  arma::mat mR_i = mR.i();
  arma::mat mD_i = mD.i();

  arma::vec vX = vY-vMu;

  for(i = 0;i<iN;i++){
    mU(i,i)   = 1.0;
    vD_s(i) =  - mD_i(i,i) - 0.5 * as_scalar(vX.t() * (-mD_i * mU * mD_i * mR_i * mD_i - mD_i * mR_i * mD_i * mU * mD_i) * vX);
    mU(i,i)   = 0.0;
  }

  return vD_s;

}
arma::vec MuScore(arma::vec vMu, arma::mat mD, arma::mat mR, arma::vec vY, int iN){

  arma::vec vU(iN);
  vU.fill(0);

  arma::mat mSigma   = mD*mR*mD;
  arma::mat mSigma_i = mSigma.i();
  arma::vec vX       = vY - vMu;

  arma::vec vMu_s(iN);
  for(int i = 0;i<iN;i++){
    vU(i) = 1.0;
    vMu_s(i) = -0.5 * as_scalar(-vU.t()*mSigma_i*vX - vX.t() * mSigma_i * vU);
    vU(i) = 0.0;
  }

  return vMu_s;

}
arma::vec mvnorm_Score(arma::vec vY, arma::vec vTheta, int iN){

  int iL = 2*iN + iN*(iN-1)/2;

  arma::vec vScore(iL);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iL-1);

  arma::mat mD = zeros(iN,iN);
  mD.diag() = vSigma;

  arma::mat mR = build_mR(vRho,iN);

  arma::vec vMu_s    = MuScore(vMu, mD, mR, vY, iN);
  arma::vec vSigma_s = DScore(mD, mR, vY, vMu, iN);
  arma::vec vRho_s   = RhoScore(vRho, mD, vY, vMu, iN) ;


  vScore.subvec(0,iN-1)    = vMu_s;
  vScore.subvec(iN,2*iN-1) = vSigma_s;
  vScore.subvec(2*iN,iL-1) = vRho_s;

  return vScore;

}
arma::vec mMVNORM_mean(arma::vec vTheta, int iN){
  arma::vec vMu    = vTheta.subvec(0,iN-1);
  return vMu;
}

arma::mat mMVNORM_cov(arma::vec vTheta, int iN){

  int iK = NumberParameters("mvnorm", iN);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iK-1);

  arma::mat mD = diagmat(vSigma);
  arma::mat mR = build_mR(vRho, iN);

  arma::mat mSigma = mD * mR * mD;

  return mSigma;
}
