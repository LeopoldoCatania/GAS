#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dGAMMA(double dY, double dAlpha, double dBeta, bool bLog=false) {

  double dLPDF = Rf_dgamma(dY, dAlpha, 1.0/dBeta, 1);

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}
//
double pGAMMA(double dY, double dAlpha, double dBeta) {

  double dP = Rf_pgamma(dY, dAlpha, 1.0/dBeta, 1, 0);

  return dP;

}
double qGAMMA(double dP, double dAlpha, double dBeta){

  double dQ = Rf_qgamma(dP, dAlpha, 1.0/dBeta, 1 ,0 );

  return dQ;

}
double rGAMMA(double dAlpha, double dBeta){
  double dY = Rf_rgamma(dAlpha, 1.0/dBeta);

  return dY;
}

arma::vec mGAMMA(double dAlpha, double dBeta){
  arma::vec vMoments(4);
  vMoments(0) = dAlpha/dBeta;
  vMoments(1) = dAlpha/pow(dBeta,2.0);
  vMoments(2) = 2.0*pow(dAlpha,0.5);
  vMoments(3) = 6.0/dAlpha + 3.0;
  return vMoments;
}

arma::vec gamma_Score(double dY, arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::vec vScore(2);

  double dAlpha_s = log(dBeta) + log(dY) - Rf_digamma(dAlpha);
  double dBeta_s  = dAlpha/dBeta - dY;

  vScore(0) = dAlpha_s;
  vScore(1) = dBeta_s;

  return vScore;

}
arma::mat gamma_IM(arma::vec vTheta){

  double dAlpha = vTheta(0);
  double dBeta  = vTheta(1);

  arma::mat mIM(2,2);

  mIM(0,0) = - Rf_trigamma(dAlpha);
  mIM(0,1) = 1.0/dBeta;
  mIM(1,0) = mIM(0,1);
  mIM(1,1) = -dAlpha/pow(dBeta,2.0);

  return mIM;
}


