#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dBER(double dY, double dPi, bool bLog=false) {

  double dLPDF = 0.0;

  if(dY == 1.0) dLPDF = log(dPi);
  if(dY == 0.0) dLPDF = log(1.0 - dPi);

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}

double pBER(double dY, double dPi) {

  double dP = 1.0 - dPi;
  if(dY == 1.0) dP = 1.0;

  return dP;

}
double qBER(double dP, double dPi){

  double dQ = 0.0;
  if(dP > 1.0 - dPi) dQ = 1.0;

  return dQ;

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
  vMoments(3) = (1.0 - 6.0 * dPi*(1.0 - dPi))/dPi*(1.0 - dPi) + 3.0;
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

