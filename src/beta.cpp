#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dBETA(double dY, double dAlpha, double dBeta, bool bLog=false) {

  double dLPDF = (dAlpha - 1.0)*log(dY) + (dBeta - 1.0)*log(1.0 - dY) + Rf_lgammafn(dAlpha + dBeta) -
                  Rf_lgammafn(dAlpha) - Rf_lgammafn(dBeta);

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}
double pBETA(double dY, double dAlpha, double dBeta) {

  double dP = Rf_pbeta(dY, dAlpha, dBeta, 1, 0);

  return dP;

}
double qBETA(double dP, double dAlpha, double dBeta) {

  double dQ = Rf_qbeta(dP,dAlpha,dBeta,1,0);

  return dQ;

}
//
double rBETA(double dAlpha, double dBeta){
  double dY = Rf_rbeta(dAlpha, dBeta);

  return dY;
}
//
arma::vec mBETA(double dAlpha, double dBeta){
  arma::vec vMoments(4);
  vMoments(0) = dAlpha/(dAlpha + dBeta);
  vMoments(1) = dAlpha*dBeta/( pow(dAlpha + dBeta, 2.0)*(dAlpha + dBeta + 1.0));
  vMoments(2) = (2.0*(dBeta - dAlpha)*pow(dAlpha + dBeta + 1.0, 0.5))/((dAlpha + dBeta + 2.0)*pow(dAlpha*dBeta,0.5));
  vMoments(3) = 6.0*( pow(dAlpha - dBeta, 2.0) * (dAlpha + dBeta + 1.0) - dAlpha*dBeta*(dAlpha + dBeta + 2.0) )/(dAlpha*dBeta*(dAlpha + dBeta + 2.0)*(dAlpha + dBeta + 3.0)) + 3.0;
  return vMoments;
}
//
arma::vec beta_Score(double dY, arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::vec vScore(2);

  double dAlpha_s = log(dY) + Rf_digamma(dAlpha + dBeta) - Rf_digamma(dAlpha);
  double dBeta_s  = log(1.0 - dY) + Rf_digamma(dAlpha + dBeta) - Rf_digamma(dBeta);

  vScore(0) = dAlpha_s;
  vScore(1) = dBeta_s;

  return vScore;

}
arma::mat beta_IM( arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::mat mIM = zeros(2,2);

  double dTrig =  Rf_trigamma(dAlpha + dBeta);

  mIM(0,0) = Rf_trigamma(dAlpha) - dTrig;
  mIM(1,1) = Rf_trigamma(dBeta)  - dTrig;
  mIM(1,0) = - dTrig;

  return mIM;
}
//
//
