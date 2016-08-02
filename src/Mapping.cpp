#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

const double dLowerShape = 4.00;
const double dUpperShape = 50.0;

const double dLowerSkewFS = 0.10;
const double dUpperSkewFS = 2.00;

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
    double dNu1   = Map(dNu1_tilde,dLowerShape,dUpperShape);//exp(dNu1_tilde) + dLowerShape;
    double dNu2   = Map(dNu2_tilde,dLowerShape,dUpperShape);//exp(dNu2_tilde) + dLowerShape;

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
    double dNu1   = Map(dNu1_tilde,dLowerShape,dUpperShape);//exp(dNu1_tilde) + dLowerShape;

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
    double dNu  = exp(dNu_tilde) + dLowerShape;

    vTheta(0) = dMu;
    vTheta(1) = dPhi;
    vTheta(2) = dNu;
  }
  if(Dist=="sstd"){

    double dMu_tilde    = vTheta_tilde(0);
    double dSigma_tilde = vTheta_tilde(1);
    double dXi_tilde    = vTheta_tilde(2);
    double dNu_tilde    = vTheta_tilde(3);


    double dMu    = dMu_tilde;
    double dSigma = exp(dSigma_tilde);
    double dXi    = Map(dXi_tilde, dLowerSkewFS, dUpperSkewFS);
    double dNu    = exp(dNu_tilde) + dLowerShape;

    vTheta(0) = dMu;
    vTheta(1) = dSigma;
    vTheta(2) = dXi;
    vTheta(3) = dNu;
  }
  if(Dist=="norm"){

    double dMu_tilde     = vTheta_tilde(0);
    double dSigma2_tilde = vTheta_tilde(1);

    double dMu     = dMu_tilde;
    double dSigma2 = exp(dSigma2_tilde);

    vTheta(0) = dMu;
    vTheta(1) = dSigma2;
  }
  if(Dist=="snorm"){

    double dMu_tilde     = vTheta_tilde(0);
    double dSigma2_tilde = vTheta_tilde(1);
    double dDelta_tilde  = vTheta_tilde(2);

    double dMu     = dMu_tilde;
    double dSigma2 = exp(dSigma2_tilde);
    double dDelta  = Map(dDelta_tilde, dLowerSkewFS, dUpperSkewFS);

    vTheta(0) = dMu;
    vTheta(1) = dSigma2;
    vTheta(2) = dDelta;

  }
  if(Dist == "poi"){
    double dMu_tilde = vTheta_tilde(0);
    double dMu       = exp(dMu_tilde);

    vTheta(0) = dMu;
  }
  if(Dist == "ber"){
    double dPi_tilde = vTheta_tilde(0);
    double dPi       = 1.0/(1.0 + exp(-dPi_tilde));

    vTheta(0) = dPi;
  }
  if(Dist == "gamma"){
    double dAlpha_tilde = vTheta_tilde(0);
    double dBeta_tilde  = vTheta_tilde(1);

    double dAlpha       = exp(dAlpha_tilde);
    double dBeta        = exp(dBeta_tilde);

    vTheta(0) = dAlpha;
    vTheta(1) = dBeta;
  }
  if(Dist == "exp"){
    double dMu_tilde = vTheta_tilde(0);
    double dMu       = exp(dMu_tilde);

    vTheta(0) = dMu;
  }
  if(Dist == "beta"){
    double dAlpha_tilde = vTheta_tilde(0);
    double dBeta_tilde  = vTheta_tilde(1);

    double dAlpha       = exp(dAlpha_tilde);
    double dBeta        = exp(dBeta_tilde);

    vTheta(0) = dAlpha;
    vTheta(1) = dBeta;
  }
  if(Dist=="ald"){

    double dTheta_tilde  = vTheta_tilde(0);
    double dSigma_tilde  = vTheta_tilde(1);
    double dKappa_tilde  = vTheta_tilde(2);

    double dTheta = dTheta_tilde;
    double dSigma = exp(dSigma_tilde);
    double dKappa = exp(dKappa_tilde);

    vTheta(0) = dTheta;
    vTheta(1) = dSigma;
    vTheta(2) = dKappa;
  }
  return InfRemover_vec(vTheta);
}
//[[Rcpp::export]]
arma::vec UnmapParameters_univ(arma::vec vTheta, std::string Dist, int iK = -9999){

  if(iK == -9999) iK = vTheta.size();

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
    double dNu1_tilde   = Unmap(dNu1,dLowerShape,dUpperShape);// log(dNu1 - dLowerShape);
    double dNu2_tilde   = Unmap(dNu2,dLowerShape,dUpperShape);// log(dNu2 - dLowerShape);

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
    double dNu1_tilde   = Unmap(dNu1,dLowerShape,dUpperShape);// log(dNu1 - dLowerShape);

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
    double dNu_tilde  = log(dNu - dLowerShape);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dPhi_tilde;
    vTheta_tilde(2) = dNu_tilde;
  }
  if(Dist=="sstd"){

    double dMu    = vTheta(0);
    double dSigma = vTheta(1);
    double dXi    = vTheta(2);
    double dNu    = vTheta(3);


    double dMu_tilde    = dMu;
    double dSigma_tilde = log(dSigma);
    double dXi_tilde    = Unmap(dXi, dLowerSkewFS, dUpperSkewFS);
    double dNu_tilde    = log(dNu - dLowerShape);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma_tilde;
    vTheta_tilde(2) = dXi_tilde;
    vTheta_tilde(3) = dNu_tilde;
  }

  if(Dist=="norm"){

    double dMu     = vTheta(0);
    double dSigma2 = vTheta(1);

    double dMu_tilde     = dMu;
    double dSigma2_tilde = log(dSigma2);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma2_tilde;
  }

  if(Dist=="snorm"){

    double dMu     = vTheta(0);
    double dSigma2 = vTheta(1);
    double dDelta  = vTheta(2);

    double dMu_tilde     = dMu;
    double dSigma2_tilde = exp(dSigma2);
    double dDelta_tilde  = Unmap(dDelta, dLowerSkewFS, dUpperSkewFS);

    vTheta_tilde(0) = dMu_tilde;
    vTheta_tilde(1) = dSigma2_tilde;
    vTheta_tilde(2) = dDelta_tilde;

  }

  if(Dist == "poi"){
    double dMu       = vTheta(0);
    double dMu_tilde = log(dMu);

    vTheta_tilde(0) = dMu_tilde;
  }
  if(Dist == "ber"){
    double dPi       = vTheta(0);
    double dPi_tilde = log(dPi / (1.0 - dPi));

    vTheta_tilde(0) = dPi_tilde;
  }
  if(Dist == "gamma"){
    double dAlpha = vTheta(0);
    double dBeta  = vTheta(1);

    double dAlpha_tilde = log(dAlpha);
    double dBeta_tilde  = log(dBeta);

    vTheta_tilde(0) = dAlpha_tilde;
    vTheta_tilde(1) = dBeta_tilde;
  }
  if(Dist == "exp"){
    double dMu       = vTheta(0);
    double dMu_tilde = log(dMu);

    vTheta_tilde(0) = dMu_tilde;
  }
  if(Dist == "beta"){
    double dAlpha = vTheta(0);
    double dBeta  = vTheta(1);

    double dAlpha_tilde = log(dAlpha);
    double dBeta_tilde  = log(dBeta);

    vTheta_tilde(0) = dAlpha_tilde;
    vTheta_tilde(1) = dBeta_tilde;
  }

  if(Dist=="ald"){

    double dTheta  = vTheta(0);
    double dSigma  = vTheta(1);
    double dKappa  = vTheta(2);

    double dTheta_tilde = dTheta;
    double dSigma_tilde = log(dSigma);
    double dKappa_tilde = log(dKappa);

    vTheta_tilde(0) = dTheta_tilde;
    vTheta_tilde(1) = dSigma_tilde;
    vTheta_tilde(2) = dKappa_tilde;
  }

  return vTheta_tilde;
}

//[[Rcpp::export]]
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
    mJ(3,3) = MapDeriv(dNu1_tilde,dLowerShape,dUpperShape);; //exp(dNu1_tilde);
    mJ(4,4) = MapDeriv(dNu2_tilde,dLowerShape,dUpperShape);; //exp(dNu2_tilde);

  }
  if(Dist=="ast1"){
    double dSigma_tilde = vTheta_tilde(1);
    double dAlpha_tilde = vTheta_tilde(2);
    double dNu1_tilde   = vTheta_tilde(3);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma_tilde);;
    mJ(2,2) = MapDeriv(dAlpha_tilde,0.01,0.99);;
    mJ(3,3) = MapDeriv(dNu1_tilde,dLowerShape,dUpperShape);; //exp(dNu1_tilde);
  }
  if(Dist=="std"){

    double dPhi_tilde = vTheta_tilde(1);
    double dNu_tilde  = vTheta_tilde(2);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dPhi_tilde);
    mJ(2,2) = exp(dNu_tilde);
  }
  if(Dist=="sstd"){

    double dSigma_tilde = vTheta_tilde(1);
    double dXi_tilde    = vTheta_tilde(2);
    double dNu_tilde    = vTheta_tilde(3);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma_tilde);
    mJ(2,2) = MapDeriv(dXi_tilde,dLowerSkewFS, dUpperSkewFS);
    mJ(3,3) = exp(dNu_tilde);
  }

  if(Dist=="norm"){

    double dSigma2_tilde = vTheta_tilde(1);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma2_tilde);
  }
  if(Dist=="snorm"){

    double dSigma2_tilde = vTheta_tilde(1);
    double dDelta_tilde  = vTheta_tilde(2);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma2_tilde);
    mJ(2,2) = MapDeriv(dDelta_tilde, dLowerSkewFS, dUpperSkewFS);;

  }
  if(Dist=="poi"){

    double dMu_tilde = vTheta_tilde(0);

    mJ(0,0) = exp(dMu_tilde);
  }
  if(Dist=="ber"){

    double dPi_tilde = vTheta_tilde(0);

    mJ(0,0) = exp(-dPi_tilde)/pow(1.0 + exp(-dPi_tilde),2.0);
  }
  if(Dist=="gamma"){

    double dAlpha_tilde = vTheta_tilde(0);
    double dBeta_tilde  = vTheta_tilde(1);

    mJ(0,0) = exp(dAlpha_tilde);
    mJ(1,1) = exp(dBeta_tilde);

  }
  if(Dist=="exp"){

    double dMu_tilde = vTheta_tilde(0);

    mJ(0,0) = exp(dMu_tilde);
  }
  if(Dist=="beta"){

    double dAlpha_tilde = vTheta_tilde(0);
    double dBeta_tilde  = vTheta_tilde(1);

    mJ(0,0) = exp(dAlpha_tilde);
    mJ(1,1) = exp(dBeta_tilde);

  }
  if(Dist=="ald"){

    double dSigma_tilde = vTheta_tilde(1);
    double dKappa_tilde  = vTheta_tilde(2);

    mJ(0,0) = 1;
    mJ(1,1) = exp(dSigma_tilde);
    mJ(2,2) = exp(dKappa_tilde);
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

//[[Rcpp::export]]
arma::mat MapR_C(arma::vec vPhi, int iN){

  arma::mat X = HalfR(vPhi);

  arma::mat R = X.t() * X;

  return(R);
}

//[[Rcpp::export]]
arma::vec UnMapR_C(arma::vec vRho, int iN){

  arma::mat mPhi = zeros(iN,iN);
  arma::mat mX   = zeros(iN,iN);
  arma::mat mR   = build_mR(vRho, iN);

  mX(0,0) = 1.0;

  int i,j,k,l;

  double dFoo1 = 0.0;
  double dFoo2 = 0.0;
  double dFoo3 = 0.0;

  for(i=1;i<iN;i++){
    mPhi(i,0) = acos(mR(i,0));
    mX(i,0)   = cos(mPhi(i,0));
  }

  if(iN>2){
    for(j=1;j<iN;j++){
      for(i=1;i<=j;i++){
        dFoo1 = 1.0;
        for(k=0;k<=i-1;k++){
          dFoo1 *= sin(mPhi(j,k));
        }
        if(i==j){
          mX(j,j) = dFoo1;
        }else{
          dFoo2 = 0.0;
          if(i>1){
            for(k=1;k<=i-1;k++){
              dFoo3 = 1.0;
              for(l=0;l<=(k-1);l++){
                dFoo3 *= sin(mPhi(i,l))*sin(mPhi(j,l));
              }
              dFoo2 += cos(mPhi(i,k))*cos(mPhi(j,k))*dFoo3;
            }
          }

          dFoo3 = 1.0;
          for(l=0;l<=(i-1);l++){
            dFoo3 *= sin(mPhi(i,l))*sin(mPhi(j,l));
          }
          mPhi(j,i) = acos((mR(i,j) - cos(mPhi(i,0))*cos(mPhi(j,0)) - dFoo2)/dFoo3);
          mX(i,j)   = cos(mPhi(j,i)) * dFoo1;
        }
      }
    }
  }

  arma::vec vPhi(iN*(iN-1.0)/2.0);
  int iC=0;

  for(j=0;j<iN-1;j++){
    for(i=j+1;i<iN;i++){
      vPhi(iC) = mPhi(i,j);
      iC++;
    }
  }

  return vPhi;
}

arma::vec mvnormMap(arma::vec vTheta_tilde, int iN, int iK){

  arma::vec vTheta(iK);

  arma::vec vMu_tilde    = vTheta_tilde.subvec(0,iN-1);
  arma::vec vSigma_tilde = vTheta_tilde.subvec(iN,2*iN-1);
  arma::vec vRho_tilde   = vTheta_tilde.subvec(2*iN,iK-1);

  arma::vec vSigma = exp(vSigma_tilde);
  vSigma           = Thresholding_vec(vSigma, 1e10);
  vSigma           = ZeroRemover_v(vSigma);

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

  vSigma = Thresholding_vec(vSigma, 1e10);
  vSigma = ZeroRemover_v(vSigma);

  double dNu       = exp(dNu_tilde) + dLowerShape;

  arma::mat mR = MapR_C(vRho_tilde, iN);

  arma::vec vR = build_vR(mR,iN);

  vTheta.subvec(0,iN-1)    = vMu_tilde;
  vTheta.subvec(iN,2*iN-1) = vSigma;
  vTheta.subvec(2*iN,iK-2) = vR;
  vTheta(iK-1)             = dNu;

  return vTheta;
}

arma::vec mvnormUnmap(arma::vec vTheta, int iN, int iK){

  arma::vec vTheta_tilde(iK);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iK-1);

  arma::vec vSigma_tilde = log(vSigma);

  arma::vec vRho_tilde = UnMapR_C(vRho, iN);

  vTheta_tilde.subvec(0,iN-1) = vMu;
  vTheta_tilde.subvec(iN,2*iN-1) = vSigma_tilde;
  vTheta_tilde.subvec(2*iN,iK-1) = vRho_tilde;

  return vTheta_tilde;
}

arma::vec mvtUnmap(arma::vec vTheta, int iN, int iK){

  arma::vec vTheta_tilde(iK);

  arma::vec vMu    = vTheta.subvec(0,iN-1);
  arma::vec vSigma = vTheta.subvec(iN,2*iN-1);
  arma::vec vRho   = vTheta.subvec(2*iN,iK-2);
  double dNu       = vTheta(iK-1);

  arma::vec vSigma_tilde = log(vSigma);

  double dNu_tilde       = log(dNu - dLowerShape);

  arma::vec vRho_tilde = UnMapR_C(vRho, iN);

  vTheta_tilde.subvec(0,iN-1)    = vMu;
  vTheta_tilde.subvec(iN,2*iN-1) = vSigma_tilde;
  vTheta_tilde.subvec(2*iN,iK-2) = vRho_tilde;
  vTheta_tilde(iK-1)             = dNu_tilde;

  return vTheta_tilde;
}

//[[Rcpp::export]]
arma::vec MapParameters_multi(arma::vec vTheta_tilde, std::string Dist,int iN, int iK){

  arma::vec vTheta(iK);

  if(Dist=="mvnorm"){
    vTheta = mvnormMap(vTheta_tilde, iN, iK);
  }
  if(Dist=="mvt"){
    vTheta = mvtMap(vTheta_tilde, iN, iK);
  }
  return InfRemover_vec(vTheta);
}
//[[Rcpp::export]]
arma::vec UnmapParameters_multi(arma::vec vTheta, std::string Dist,int iN, int iK){

  arma::vec vTheta_tilde(iK);

  if(Dist=="mvnorm"){
    vTheta = mvnormUnmap(vTheta, iN, iK);
  }
  if(Dist=="mvt"){
    vTheta = mvtUnmap(vTheta, iN, iK);
  }
  return InfRemover_vec(vTheta);
}

//[[Rcpp::export]]
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

  //
  // int iK = iN*(iN-1)/2;
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
    mJ(3,0) = -mC(0,2) * sin(mPhi(0,1)) + mS(0,2)*mC(1,2)*cos(mPhi(0,1));
    mJ(3,1) = -mC(0,1) * sin(mPhi(0,2)) + mS(0,1)*mC(1,2)*cos(mPhi(0,2));
    mJ(2,2) = -mS(0,1)*mS(0,2)*sin(mPhi(1,2));
    //
    mJ(3,3) = -sin(mPhi(0,3));
    mJ(4,0) = -mC(0,3)*sin(mPhi(0,1)) + mC(1,3)*mS(0,1)*cos(mPhi(0,1));
    mJ(4,2) = -sin(mPhi(0,3))*mC(0,1) + mC(1,3)*cos(mPhi(0,3))*mS(0,1);
    mJ(4,4) = sin(mPhi(1,3))*mS(0,3)*mS(0,1);
    mJ(5,1) = -sin(mPhi(0,2))*mC(0,3) + mC(1,2)*mC(1,3)*cos(mPhi(0,2))*mS(0,3) + mS(1,2)*cos(mPhi(0,2))*mS(0,3)*mS(1,3)*mC(2,3);
    mJ(5,3) = -sin(mPhi(1,2))*mC(1,3)*mS(0,2)*mS(0,3) + cos(mPhi(1,2))*mS(0,2)*mS(0,3)*mS(1,3)*mC(2,3);
    mJ(5,2) = -sin(mPhi(0,3))*mC(0,2) + cos(mPhi(0,3))*mC(1,2)*mC(1,3)*mS(0,2) + cos(mPhi(0,3))*mS(1,2)*mS(0,2)*mS(1,3)*mC(2,3);
    mJ(5,4) = -sin(mPhi(1,3))*mC(1,2)*mS(0,2)*mS(0,3) + cos(mPhi(1,3))*mS(1,2)*mS(0,2)*mS(0,3)*mC(2,3);
    mJ(5,5) = -sin(mPhi(2,3))*mS(1,2)*mS(0,2)*mS(0,3)*mS(1,3);
  }

  if(iN>4) mJ.diag().ones();
  // mJ.diag().ones();
  return mJ;
}


//[[Rcpp::export]]
arma::vec IndexesFinder(int iC, int iN){
  int l,m,iC_c = 0;

  arma::vec vIndexes(2);

  for(l = 0; l<iN; l++){
    for(m = 0; m<=l; m++){
      if(l!=m){
        if(iC_c == iC){
          vIndexes(0) = l;
          vIndexes(1) = m;
        }
        iC_c++;
      }
    }
  }

  return vIndexes;
}


// //[[Rcpp::export]]
// arma::mat Jacobian_MapR2(arma::vec vPhi, int iN){
//
//   int i,j,iC=0,iK = iN*(iN-1)/2;
//
//   arma::mat mPhi = zeros(iN,iN);
//   arma::mat mS   = zeros(iN,iN);
//   arma::mat mJ   = zeros(iK,iK);
//
//   for(i = 0;i<iN;i++){
//     for(j=0;j<=i;j++){
//       if(i!=j){
//         mPhi(i,j) = vPhi(iC);
//         mS(i,j)   = sin(vPhi(iC));
//         iC+=1;
//       }
//     }
//   }
//
//   double dFoo = 1.0;
//
//   arma::vec vIndexes_num(2);
//   arma::vec vIndexes_den(2);
//
//   //diagonal elements before
//   for(i = 0;i<iN-1;i++){
//     mJ(i,i) = -sin(vPhi(i));
//   }
//
//   int l,m,h,k;
//
//   for(i = iN-1;i<iK;i++){
//     for(j=0;j<=i;j++){
//       vIndexes_num = IndexesFinder(i, iN);
//       vIndexes_den = IndexesFinder(j, iN);
//
//       l = vIndexes_num(0);
//       m = vIndexes_num(1);
//       h = vIndexes_den(0);
//       k = vIndexes_num(0);
//
//
//     }
//   }
//
//   return mJ;
//
// }

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
  double dNu_tilde       = vTheta_tilde(iK-1);

  arma::mat mJ=eye(iK,iK);

  mJ.submat(iN,iN,2*iN-1,2*iN-1) = Jacobian_MapD(vSigma_tilde, iN);
  mJ.submat(2*iN,2*iN,iK-2,iK-2)     = Jacobian_MapR(vRho_tilde, iN);
  mJ(iK-1,iK-1) = exp(dNu_tilde);

  return mJ;

}

//[[Rcpp::export]]
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



