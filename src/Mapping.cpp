#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double Map(double dX, double dL,double dU) {
  double dMap =  dL + ( dU - dL ) / (1.0 + exp( - dX ));
  return dMap;
}
double Unmap(double dG, double dL,double dU) {
  double dUnmap = log((dG-dL)/(dU-dG));
  return dUnmap;
}
//[[Rcpp::export]]
arma::vec Map_Vec(arma::vec vX, double dL ,double dU) {
  arma::vec vMap =  dL + ( dU - dL ) / (1.0 + exp( - vX ));
  return vMap;
}
//[[Rcpp::export]]
arma::vec unmapVec_C(arma::vec vG, double dL, double dU) {
  arma::vec vUnmap = log((vG-dL)/(dU-vG));
  return vUnmap;
}
double MapDeriv(double dX, double dL, double dU){
  // double dDeriv=exp(-dX)*(dU-dL)/( pow(1.0+exp(-dX),2.0) );

  double dDeriv = -dX + log(dU-dL) - 2.0*log(1.0 + exp(-dX));

  return exp(dDeriv);
}

//######################### UNIVARIATE #####################################
//[[Rcpp::export]]
arma::vec MapParameters_univ(arma::vec vTheta_tilde, std::string Dist, int iK){

  arma::vec vTheta(iK);

  if(Dist=="ast"){
    double dMu_tilde    = vTheta_tilde(0);
    double dSigma_tilde = vTheta_tilde(1);
    double dAlpha_tilde = vTheta_tilde(2);
    double dNu1_tilde   = vTheta_tilde(3);
    double dNu2_tilde   = vTheta_tilde(4);

    double dMu    = dMu_tilde;
    double dSigma = exp(dSigma_tilde);
    double dAlpha = Map(dAlpha_tilde,0.01,0.99);
    double dNu1   = exp(dNu1_tilde) + 4.01;
    double dNu2   = exp(dNu2_tilde) + 4.01;

    vTheta(0) = dMu;
    vTheta(1) = dSigma;
    vTheta(2) = dAlpha;
    vTheta(3) = dNu1;
    vTheta(4) = dNu2;

  }
  if(Dist=="ast1"){
    double dMu_tilde    = vTheta_tilde(0);
    double dSigma_tilde = vTheta_tilde(1);
    double dAlpha_tilde = vTheta_tilde(2);
    double dNu1_tilde   = vTheta_tilde(3);

    double dMu    = dMu_tilde;
    double dSigma = exp(dSigma_tilde);
    double dAlpha = Map(dAlpha_tilde,0.01,0.99);
    double dNu1   = exp(dNu1_tilde) + 4.01;

    vTheta(0) = dMu;
    vTheta(1) = dSigma;
    vTheta(2) = dAlpha;
    vTheta(3) = dNu1;
  }
  if(Dist=="std"){

    double dMu_tilde  = vTheta_tilde(0);
    double dPhi_tilde = vTheta_tilde(1);
    double dNu_tilde  = vTheta_tilde(2);

    double dMu  = dMu_tilde;
    double dPhi = exp(dPhi_tilde);
    double dNu  = exp(dNu_tilde) + 2.01;

    vTheta(0) = dMu;
    vTheta(1) = dPhi;
    vTheta(2) = dNu;
  }
  if(Dist=="norm"){

    double dMu_tilde     = vTheta_tilde(0);
    double dSigma2_tilde = vTheta_tilde(1);

    double dMu     = dMu_tilde;
    double dSigma2 = exp(dSigma2_tilde);

    vTheta(0) = dMu;
    vTheta(1) = dSigma2;
  }

  return InfRemover_vec(vTheta);
}

//[[Rcpp::export]]
arma::vec UnmapParameters_univ(arma::vec vTheta, std::string Dist, int iK){

  arma::vec vTheta_tilde(iK);

  if(Dist=="ast"){
    double dMu    = vTheta(0);
    double dSigma = vTheta(1);
    double dAlpha = vTheta(2);
    double dNu1   = vTheta(3);
    double dNu2   = vTheta(4);

    double dMu_tilde    = dMu;
    double dSigma_tilde = log(dSigma);
    double dAlpha_tilde = Unmap(dAlpha,0.01,0.99);
    double dNu1_tilde   = log(dNu1 - 4.01);
    double dNu2_tilde   = log(dNu2 - 4.01);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma_tilde;
    vTheta_tilde(2) = dAlpha_tilde;
    vTheta_tilde(3) = dNu1_tilde;
    vTheta_tilde(4) = dNu2_tilde;

  }
  if(Dist=="ast1"){
    double dMu    = vTheta(0);
    double dSigma = vTheta(1);
    double dAlpha = vTheta(2);
    double dNu1   = vTheta(3);

    double dMu_tilde    = dMu;
    double dSigma_tilde = log(dSigma);
    double dAlpha_tilde = Unmap(dAlpha,0.01,0.99);
    double dNu1_tilde   = log(dNu1 - 4.01);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma_tilde;
    vTheta_tilde(2) = dAlpha_tilde;
    vTheta_tilde(3) = dNu1_tilde;

  }
  if(Dist=="std"){

    double dMu  = vTheta(0);
    double dPhi = vTheta(1);
    double dNu  = vTheta(2);

    double dMu_tilde  = dMu;
    double dPhi_tilde = log(dPhi);
    double dNu_tilde  = log(dNu - 2.01);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dPhi_tilde;
    vTheta_tilde(2) = dNu_tilde;
  }

  if(Dist=="norm"){

    double dMu     = vTheta(0);
    double dSigma2 = vTheta(1);

    double dMu_tilde     = dMu;
    double dSigma2_tilde = log(dSigma2);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma2_tilde;
  }

  return vTheta_tilde;
}
arma::mat MapParametersJacobian_univ(arma::vec vTheta_tilde, std::string Dist, int iK){

  arma::mat mJ=zeros(iK,iK);

  if(Dist=="ast"){
    double dSigma_tilde = vTheta_tilde(1);
    double dAlpha_tilde = vTheta_tilde(2);
    double dNu1_tilde   = vTheta_tilde(3);
    double dNu2_tilde   = vTheta_tilde(4);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma_tilde);;
    mJ(2,2) = MapDeriv(dAlpha_tilde,0.01,0.99);;
    mJ(3,3) = exp(dNu1_tilde);
    mJ(4,4) = exp(dNu2_tilde);

  }
  if(Dist=="ast1"){
    double dSigma_tilde = vTheta_tilde(1);
    double dAlpha_tilde = vTheta_tilde(2);
    double dNu1_tilde   = vTheta_tilde(3);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma_tilde);;
    mJ(2,2) = MapDeriv(dAlpha_tilde,0.01,0.99);;
    mJ(3,3) = exp(dNu1_tilde);
  }
  if(Dist=="std"){

    double dPhi_tilde = vTheta_tilde(1);
    double dNu_tilde  = vTheta_tilde(2);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dPhi_tilde);
    mJ(2,2) = exp(dNu_tilde);
  }

  if(Dist=="norm"){

    double dSigma2_tilde = vTheta_tilde(1);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma2_tilde);
  }
  arma::vec vJ_safe =  InfRemover_vec(mJ.diag());
  mJ.diag() = vJ_safe;
  return mJ;
}

//######################### MULTIVARIATE #####################################

arma::mat HalfR(arma::vec vPhi){
  int k = vPhi.size();

  int n=(1+sqrt(1+8*k))/2;

  arma::mat phi=FillUpperTriangular(vPhi,n);

  arma::mat c=cos(phi);
  arma::mat s=sin(phi);

  arma::mat foo=cumprodMat_removeLastRow(s);
  arma::vec baz(n);
  baz.fill(1);

  arma::mat foo2=Up_rbind_C(foo,baz);

  arma::mat X = c % foo2;
  return X;
}

arma::mat MapR_C(arma::vec vPhi, int iN){

  arma::mat X = HalfR(vPhi);

  arma::mat R = X.t() * X;

  return(R);
}
//[[Rcpp::export]]
arma::vec UnMapR_C(arma::vec vRho, int iN){

  arma::vec vPhi(iN*(iN-1)/2);

  if(iN==2){
    vPhi(0) = acos(vRho(0));
  }
  if(iN==3){
    vPhi(0) = acos(vRho(0));
    vPhi(1) = acos(vRho(1));
    vPhi(2) = acos((vRho(2)-vRho(0)*vRho(1))/(sin(vPhi(0))*sin(vPhi(1))));
  }
  if(iN==4){
    vPhi(0) = acos(vRho(0));
    vPhi(1) = acos(vRho(1));
    vPhi(2) = acos((vRho(2)-vRho(0)*vRho(1))/(sin(vPhi(0))*sin(vPhi(1))));
    vPhi(3) = acos(vRho(3));
    vPhi(4) = acos(  (vRho(4)-cos(vPhi(0))*cos(vPhi(3)))/(sin(vPhi(3))*sin(vPhi(0))));
    vPhi(5) = acos( (vRho(5)-cos(vPhi(1))*cos(vPhi(3)) - cos(vPhi(2))*cos(vPhi(4))*sin(vPhi(1))*sin(vPhi(3)))/(sin(vPhi(1))*sin(vPhi(3))*sin(vPhi(2))*sin(vPhi(4))));
  }

  return vPhi;
}
arma::vec mvnormMap(arma::vec vTheta_tilde, int iN, int iK){

  arma::vec vTheta(iK);

  arma::vec vMu_tilde    = vTheta_tilde.subvec(0,iN-1);
  arma::vec vSigma_tilde = vTheta_tilde.subvec(iN,2*iN-1);
  arma::vec vRho_tilde   = vTheta_tilde.subvec(2*iN,iK-1);

  arma::vec vSigma = exp(vSigma_tilde);

  arma::mat mR = MapR_C(vRho_tilde, iN);

  arma::vec vR = build_vR(mR,iN);

  vTheta.subvec(0,iN-1) = vMu_tilde;
  vTheta.subvec(iN,2*iN-1) = vSigma;
  vTheta.subvec(2*iN,iK-1) = vR;

  return vTheta;
}
arma::vec mvtMap(arma::vec vTheta_tilde, int iN, int iK){

  arma::vec vTheta(iK);

  arma::vec vMu_tilde    = vTheta_tilde.subvec(0,iN-1);
  arma::vec vSigma_tilde = vTheta_tilde.subvec(iN,2*iN-1);
  arma::vec vRho_tilde   = vTheta_tilde.subvec(2*iN,iK-2);
  double dNu_tilde       = vTheta_tilde(iK-1);

  arma::vec vSigma = exp(vSigma_tilde);

  vSigma = InfRemover_vec(vSigma);
  vSigma = ZeroRemover_v(vSigma);

  double dNu       = exp(dNu_tilde) + 2.01;

  arma::mat mR = MapR_C(vRho_tilde, iN);

  arma::vec vR = build_vR(mR,iN);

  vTheta.subvec(0,iN-1)    = vMu_tilde;
  vTheta.subvec(iN,2*iN-1) = vSigma;
  vTheta.subvec(2*iN,iK-2) = vR;
  vTheta(iK-1)             = dNu;

  return vTheta;
}

//[[Rcpp::export]]
arma::vec MapParameters_multi(arma::vec vTheta_tilde, std::string Dist,int iN, int iK){

  arma::vec vTheta(iK);

  if(Dist=="mvnorm"){
    vTheta = mvnormMap(vTheta_tilde, iN, iK);
  }
  if(Dist=="std"){
    vTheta = mvtMap(vTheta_tilde, iN, iK);
  }
  return InfRemover_vec(vTheta);
}

arma::mat Jacobian_MapR(arma::vec vPhi, int iN){

  int i,j,iC=0,iK = iN*(iN-1)/2;

  arma::mat mPhi=zeros(iN,iN);
  arma::mat mC=zeros(iN,iN);
  arma::mat mS=zeros(iN,iN);

  for(i = 0;i<iN;i++){
    for(j = i;j<iN;j++){
      if(i!=j){
        mPhi(i,j) = vPhi(iC);
        mC(i,j) = cos(vPhi(iC));
        mS(i,j) = sin(vPhi(iC));
        iC+=1;
      }
    }
  }


  arma::mat mJ = zeros(iK,iK);

  if(iN==2){
    mJ(0,0) = -sin(mPhi(0,1));
  }

  if(iN==3){
    mJ(0,0) = -sin(mPhi(0,1));
    //
    mJ(1,1) = -sin(mPhi(0,2));
    mJ(2,0) = -mC(0,2) * sin(mPhi(0,1)) + mS(0,2)*mC(1,2)*cos(mPhi(0,1));
    mJ(2,1) = -mC(0,1) * sin(mPhi(0,2)) + mS(0,1)*mC(1,2)*cos(mPhi(0,2));
    mJ(2,2) = -mS(0,1)*mS(0,2)*sin(mPhi(1,2));
  }
  if(iN==4){
    mJ(0,0) = -sin(mPhi(0,1));
    mJ(1,1) = -sin(mPhi(0,2));
    mJ(2,0) = -mC(0,2) * sin(mPhi(0,1)) + mS(0,2)*mC(1,2)*cos(mPhi(0,1));
    mJ(2,1) = -mC(0,1) * sin(mPhi(0,2)) + mS(0,1)*mC(1,2)*cos(mPhi(0,2));
    mJ(2,2) = -mS(0,1)*mS(0,2)*sin(mPhi(1,2));
    //
    mJ(3,3) = -sin(mPhi(0,3));
    mJ(4,0) = -mC(0,3)*sin(mPhi(0,1)) + mC(1,3)*mS(0,1)*cos(mPhi(0,1));
    mJ(4,3) = -sin(mPhi(0,3))*mC(0,1) + mC(1,3)*cos(mPhi(0,3))*mS(0,1);
    mJ(4,4) = sin(mPhi(1,3))*mS(0,3)*mS(0,1);
    mJ(5,1) = -sin(mPhi(0,2))*mC(0,3) + mC(1,2)*mC(1,3)*cos(mPhi(0,2))*mS(0,3) + mS(1,2)*cos(mPhi(0,2))*mS(0,3)*mS(1,3)*mC(2,3);
    mJ(5,2) = -sin(mPhi(1,2))*mC(1,3)*mS(0,2)*mS(0,3) + cos(mPhi(1,2))*mS(0,2)*mS(0,3)*mS(1,3)*mC(2,3);
    mJ(5,3) = -sin(mPhi(0,3))*mC(0,2) + cos(mPhi(0,3))*mC(1,2)*mC(1,3)*mS(0,2) + cos(mPhi(0,3))*mS(1,2)*mS(0,2)*mS(1,3)*mC(2,3);
    mJ(5,4) = -sin(mPhi(1,3))*mC(1,2)*mS(0,2)*mS(0,3) + cos(mPhi(1,3))*mS(1,2)*mS(0,2)*mS(0,3)*mC(2,3);
    mJ(5,5) = -sin(mPhi(2,3))*mS(1,2)*mS(0,2)*mS(0,3)*mS(1,3);
  }

  if(iN>4) mJ.diag().ones();

  return mJ;
}
arma::mat Jacobian_MapD(arma::vec vSigma_tilde, int iN){

  arma::mat mJ = zeros(iN,iN);

  for(int i=0;i<iN;i++) mJ(i,i) = exp(vSigma_tilde(i));

  return mJ;
}
arma::mat Jacobian_mvnormMap(arma::vec vTheta_tilde, int iN, int iK){

  arma::vec vSigma_tilde = vTheta_tilde.subvec(iN,2*iN-1);
  arma::vec vRho_tilde   = vTheta_tilde.subvec(2*iN,iK - 1);

  arma::mat mJ=eye(iK,iK);

  mJ.submat(iN,iN,2*iN-1,2*iN-1) = Jacobian_MapD(vSigma_tilde, iN);
  mJ.submat(2*iN,2*iN,iK-1,iK-1)     = Jacobian_MapR(vRho_tilde, iN);

  return mJ;

}
arma::mat Jacobian_mvtMap(arma::vec vTheta_tilde, int iN, int iK){

  arma::vec vSigma_tilde = vTheta_tilde.subvec(iN,2*iN-1);
  arma::vec vRho_tilde   = vTheta_tilde.subvec(2*iN,iK - 2);
  double dNu             = vTheta_tilde(iK-1);

  arma::mat mJ=eye(iK,iK);

  mJ.submat(iN,iN,2*iN-1,2*iN-1) = Jacobian_MapD(vSigma_tilde, iN);
  mJ.submat(2*iN,2*iN,iK-2,iK-2)     = Jacobian_MapR(vRho_tilde, iN);
  mJ(iK-1,iK-1) = dNu;

  return mJ;

}

arma::mat MapParametersJacobian_multi(arma::vec vTheta_tilde, std::string Dist, int iN, int iK){

  arma::mat mJ(iK,iK);

  if(Dist == "mvnorm"){
    mJ = Jacobian_mvnormMap(vTheta_tilde, iN,iK);
  }
  if(Dist == "mvt"){
    mJ = Jacobian_mvtMap(vTheta_tilde, iN,iK);
  }

  return mJ;

}



