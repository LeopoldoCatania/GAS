#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dEXP(double dY, double dMu, bool bLog=false) {

  double dLPDF = log(dMu) - dMu*dY;

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}
double pEXP(double dY, double dMu) {

  double dP = Rf_pexp(dY,dMu,1,0);

  return dP;

}
double qEXP(double dP, double dMu) {

  double dQ = Rf_qexp(dP,dMu, 1,0);

  return dQ;

}
//
double rEXP(double dMu){
  double dY = Rf_rexp(dMu);

  return dY;
}
//
arma::vec mEXP(double dMu){
  arma::vec vMoments(4);
  vMoments(0) = 1.0/dMu;
  vMoments(1) = 1.0/pow(dMu, 2.0);
  vMoments(2) = 2.0;
  vMoments(3) = 9;
  return vMoments;
}
//
arma::vec exp_Score(double dY, double dMu){

  arma::vec vScore(1);

  vScore(0) = 1.0/dMu - dY;

  return vScore;

}
arma::mat exp_IM(double dMu){

  arma::mat mIM = zeros(1,1);

  mIM(0,0) = 1.0/pow(dMu,2.0);

  return mIM;
}
