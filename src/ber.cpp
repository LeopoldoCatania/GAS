#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dBER(double dY, double dPi, bool bLog=false) {
  return Rf_dbinom(dY,1.0,dPi,bLog ? 1 : 0);
}

double pBER(double dY, double dPi) {
  return Rf_pbinom(dY,1.0,dPi,1,0);
}
double qBER(double dP, double dPi) {
  return Rf_qbinom(dP,1.0,dPi,1,0);
}

double rBER(double dPi){

  double dU = Rf_runif(0,1);

  double dY = 1.0;

  if(dU > dPi) dY = 0.0;

  return dY;
}

arma::vec mBER(double dPi){
  arma::vec vMoments(4);
  vMoments(0) = dPi;
  vMoments(1) = dPi*(1.0 - dPi);
  vMoments(2) = (1.0 - 2.0 * dPi)/pow(dPi*(1.0 - dPi), 0.5);
  vMoments(3) =  3.0 + (1.0 - 6.0 * dPi * (1.0 - dPi)) / (dPi * (1.0 - dPi));
  return vMoments;
}

arma::vec ber_Score(double dY, double dPi){

  arma::vec vScore(1);

  double dPi_s = (dY - dPi)/(dPi * (1.0 - dPi));

  vScore(0) = dPi_s;

  return vScore;

}
arma::mat ber_IM(double dPi){

  arma::mat mIM=zeros(1,1);

  mIM(0,0) = 1.0 / (dPi * (1.0 - dPi));

  return mIM;
}

