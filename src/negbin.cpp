#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dNEGBIN(double dY, double dPi, double dNu, bool bLog = false) {

  double dLPDF = Rf_dnbinom(dY, dNu, dPi, 1);

  if(!bLog) dLPDF = exp(dLPDF);

  return dLPDF;

}

double pNEGBIN(double dY, double dPi, double dNu) {

  double dP = Rf_pnbinom(dY, dNu, dPi, 1, 0);

  return dP;

}
double qNEGBIN(double dP, double dPi, double dNu){

  double dQ = Rf_qnbinom(dP, dNu, dPi, 1, 0);

  return dQ;

}
double rNEGBIN(double dPi, double dNu){
  double dY = Rf_rnbinom(dNu, dPi);

  return dY;
}

arma::vec mNEGBIN(double dPi, double dNu){
  arma::vec vMoments(4);
  vMoments(0) = dNu * (1.0 - dPi)/dPi;
  vMoments(1) = vMoments(0) + pow(vMoments(0), 2.0)/dNu;
  vMoments(2) = (1.0 + dPi)/pow(dPi * dNu, 0.5);
  vMoments(3) = 6.0/dNu + pow(1.0 - dPi, 2.0)/(dPi * dNu) + 3.0;
  return vMoments;
}

arma::vec negbin_Score(double dY, arma::vec vTheta) {

  double dPi = vTheta(0);
  double dNu = vTheta(1);

  arma::vec vScore(2);

  double dPi_s = dNu/dPi - dY/(1.0 - dPi);
  double dNu_s = Rf_digamma(dNu + dY) - Rf_digamma(dNu) + log(dPi);

  vScore(0) = dPi_s;
  vScore(1) = dNu_s;

  return vScore;

}
arma::mat negbin_IM(arma::vec vTheta) {

  double dPi = vTheta(0);
  double dNu = vTheta(1);

  arma::mat mIM = zeros(2, 2);

  mIM(0, 0) = dNu/(pow(dPi, 2.0) * (1.0 - dPi));
  mIM(1, 1) = 1.0;
  //NOTE: WE ARE DOING THIS BY PURPOSE. SCALING IS ALLOWED ONLY FOR \pi. 19/05/2017.

  return mIM;
}


