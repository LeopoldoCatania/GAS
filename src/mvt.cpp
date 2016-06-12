#include <RcppArmadillo.h>
#include "mvnorm.h"
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double dmvt(arma::vec vY,
            arma::vec vMu,
            arma::mat mSigma,
            double dNu,
            bool dLog= false) {

  double dDetSigma = det(mSigma);
  if(dDetSigma<1e-50){
    dDetSigma = 1e-50;
  }
  double dN = vMu.size()*1.0;

  arma::vec vZ = arma::inv(arma::trans(chol(mSigma)))*(vY-vMu);

  double dFoo = as_scalar(vZ.t()*vZ);

  double dLLK = Rf_lgammafn((dNu + dN)*0.5) - Rf_lgammafn(dNu*0.5) - 0.5*dN * log(dNu) - 0.5*dN*log(M_PI*1.0) - 0.5*log(dDetSigma) -
    0.5*(dN+dNu)*log(1.0 + dFoo/dNu);

  if (!dLog) {
    dLLK = exp(dLLK);
  }
  return dLLK;

}

double rigamma_d(double dA, double dB){
  double dOut = (1/Rf_rgamma(dA,1.0/dB));
  return dOut;
}

arma::vec rigamma_vec(int iN, double dA, double dB){
  arma::vec vOut(iN);
  for(int i=0;i<iN;i++){
    vOut(i)=rigamma_d(dA, dB);
  }
  return vOut;
}

arma::mat rmvt_mat(int iN, arma::vec vMu, arma::mat mSigma, double dNu) {
  int dD=vMu.size();
  arma::mat mZ = rmvnorm_mat(iN, vMu, mSigma);
  arma::vec mW = rigamma_vec(iN, 0.5*dNu, 0.5*dNu);
  arma::mat mY(iN,dD);

  for(int i=0;i<iN;i++){
    mY.row(i)=vMu.t() + sqrt(mW(i))*mZ.row(i);
  }
  return mY;
}

arma::vec RhoScore_mvt(arma::vec vR, arma::mat mD, arma::vec vY, arma::vec vMu, double dNu, int iN){

  arma::mat mRho_S = zeros(iN,iN);
  arma::mat mU     = zeros(iN,iN);

  arma::mat mR = build_mR(vR,iN);

  int i, j;
  arma::vec vV   = mD.i() * (vY - vMu);
  arma::mat mR_i = mR.i();

  arma::mat mSigma = mD*mR*mD;

  arma::vec vZ = arma::inv(arma::trans(chol(mSigma)))*(vY-vMu);

  double dFoo = as_scalar(vZ.t()*vZ);

  double dConst = (dNu + iN*1.0)/(2.0*(1.0 + dFoo/dNu)*dNu);

  for(i = 0;i<iN;i++){
    for(j = 0;j<iN;j++){
      if(i!=j){
        mU(i,j) = 1.0;
        mRho_S(i,j) = as_scalar(vV.t() * mR_i * mU * mR_i * vV);
        mU(i,j) = 0.0;
      }
    }
  }

  mRho_S = mRho_S*dConst*2.0 - mR_i.t();

  arma::vec vR_s = build_vR(mRho_S, iN);

  return vR_s;

}
arma::vec MuScore_mvt(arma::vec vMu, arma::mat mD, arma::mat mR, arma::vec vY, double dNu,int iN){

  arma::vec vU(iN);
  vU.fill(0);

  arma::mat mSigma   = mD*mR*mD;
  arma::mat mSigma_i = mSigma.i();
  arma::vec vX       = vY - vMu;

  arma::vec vZ = arma::inv(arma::trans(chol(mSigma)))*(vY-vMu);

  double dFoo = as_scalar(vZ.t()*vZ);

  double dConst = (dNu + iN*1.0)/(2.0*(1.0 + dFoo/dNu)*dNu);

  arma::vec vMu_s(iN);
  for(int i = 0;i<iN;i++){
    vU(i) = 1.0;
    vMu_s(i) = as_scalar(vU.t()*mSigma_i*vX);
    vU(i) = 0.0;
  }

  return vMu_s*dConst*2.0;

}
arma::vec DScore_mvt(arma::mat mD, arma::mat mR, arma::vec vY, arma::vec vMu, double dNu, int iN){

  arma::mat mU = zeros(iN,iN);

  arma::vec vD_s(iN);

  int i;
  arma::mat mR_i = mR.i();
  arma::mat mD_i = mD.i();
  arma::mat mSigma   = mD*mR*mD;

  arma::vec vX = vY-vMu;

  arma::vec vZ = arma::inv(arma::trans(chol(mSigma)))*(vY-vMu);

  double dFoo = as_scalar(vZ.t()*vZ);

  double dConst = (dNu + iN*1.0)/(2.0*(1.0 + dFoo/dNu)*dNu);

  for(i = 0;i<iN;i++){
    mU(i,i)   = 1.0;
    vD_s(i) =  - mD_i(i,i) - dConst* as_scalar(vX.t() * (-mD_i * mU * mD_i * mR_i * mD_i - mD_i * mR_i * mD_i * mU * mD_i) * vX);
    mU(i,i)   = 0.0;
  }

  return vD_s;

}
double NuScore_mvt(arma::mat mD, arma::mat mR, arma::vec vY, arma::vec vMu, double dNu, int iN){

  arma::mat mSigma   = mD*mR*mD;
  arma::vec vZ = arma::inv(arma::trans(chol(mSigma)))*(vY-vMu);

  double dK = as_scalar(vZ.t()*vZ);

  double dNu_s = 0.5*Rf_digamma(0.5*(dNu+iN)) -0.5*Rf_digamma(0.5*dNu) - 0.5*iN/dNu -
    0.5*( log(1.0 + dK/dNu) - dK*(dNu+iN)/(dNu*dNu*(1+dK/dNu)) );

  return dNu_s;

}
arma::vec mvt_Score(arma::vec vTheta, arma::vec vY, int iN){

  int iL = 2*iN + iN*(iN-1)/2 + 1;

  arma::vec vScore(iL);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iL-2);
  double dNu       = vTheta(iL-1);

  arma::mat mD = zeros(iN,iN);
  mD.diag() = vSigma;

  arma::mat mR = build_mR(vRho,iN);

  arma::vec vMu_s    = MuScore_mvt(vMu, mD, mR, vY,dNu, iN);
  arma::vec vSigma_s = DScore_mvt(mD, mR, vY, vMu,dNu, iN);
  arma::vec vRho_s   = RhoScore_mvt(vRho, mD, vY, vMu,dNu, iN) ;
  double dNu_s       = NuScore_mvt(mD, mR, vY, vMu, dNu, iN) ;


  vScore.subvec(0,iN-1)    = vMu_s;
  vScore.subvec(iN,2*iN-1) = vSigma_s;
  vScore.subvec(2*iN,iL-2) = vRho_s;
  vScore(iL-1)             = dNu_s;

  return vScore;

}


