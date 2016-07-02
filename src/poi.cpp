#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dPOI(double dY, double dMu, bool bLog=false) {

  double dLPDF = Rf_dpois(dY,dMu,1);

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}

double pPOI(double dY, double dMu) {

  double dP = Rf_ppois(dY, dMu, 1, 0);

  return dP;

}
double qPOI(double dP, double dMu){

  double dQ = Rf_qpois(dP, dMu, 1, 0);

  return dQ;

}
double rPOI(double dMu){
  double dY = Rf_rpois(dMu);

  return dY;
}

arma::vec mPOI(double dMu){
  arma::vec vMoments(4);
  vMoments(0) = dMu;
  vMoments(1) = dMu;
  vMoments(2) = pow(dMu, - 0.5);
  vMoments(3) = 1.0/dMu;
  return vMoments;
}

arma::vec poi_Score(double dY, double dMu){

  arma::vec vScore(1);

  double dMu_s = -1.0 + dY / dMu;

  vScore(0) = dMu_s;

  return vScore;

}
arma::mat poi_IM(double dMu){

  arma::mat mIM=zeros(1,1);

  mIM(0,0) = 1.0 / dMu;

  return mIM;
}


