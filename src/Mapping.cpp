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
//[[Rcpp::export]]
arma::vec MapParameters(arma::vec vTheta_tilde, std::string Dist, int iK){

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
arma::vec UnmapParameters(arma::vec vTheta, std::string Dist, int iK){

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
arma::mat MapParametersJacobian(arma::vec vTheta_tilde, std::string Dist, int iK){

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
