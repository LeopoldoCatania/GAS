#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dEXP(double dY, double dMu, bool bLog=false) {
  if (!R_FINITE(dMu) || dMu <= 0) return NA_REAL;
  return Rf_dexp(dY, 1.0/dMu, bLog ? 1 : 0);
}
double pEXP(double dY, double dMu) {

  // dMu is the rate; Rmath expects the scale.
  return Rf_pexp(dY, 1.0 / dMu, 1, 0);

}
double qEXP(double dP, double dMu) {

  return Rf_qexp(dP, 1.0 / dMu, 1, 0);

}
double rEXP(double dMu) {

  return Rf_rexp(1.0 / dMu);

}

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
