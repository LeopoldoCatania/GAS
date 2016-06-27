#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

arma::vec norm_Score(double dY, arma::vec vTheta){

  double dMu=vTheta(0);
  double dSigma2=vTheta(1);

  double dMu_s     =(dY - dMu)/dSigma2;

  double dSigma2_s =-0.5*(1.0 - pow(dY-dMu,2)/dSigma2 )/dSigma2;

  arma::vec vScore(2);

  vScore(0)=dMu_s;
  vScore(1)=dSigma2_s;

  return vScore;

}
arma::mat norm_IM(arma::vec vTheta){

  double dSigma2 = vTheta(1);

  arma::mat mIM=zeros(2,2);

  mIM(0,0) = 1.0/dSigma2;
  mIM(1,1) = 1.0/(2.0*pow(dSigma2,2.0));

  return mIM;

}

double dNORM(double dY, double dMu, double dSigma2, bool bLog = false){

  double dLPDF = -0.5 * log(2.0 * M_PI) - 0.5 * log(dSigma2) - 0.5*pow(dY-dMu,2.0)/dSigma2;

  if(!bLog){
    dLPDF = exp(dLPDF);
  }

  return dLPDF;

}

arma::vec mNORM(double dMu, double dSigma2){
  arma::vec vMoments(4);
  vMoments(0) = dMu;
  vMoments(1) = dSigma2;
  vMoments(2) = 0.0;
  vMoments(3) = 3.0*pow(dSigma2,2.0);
  return vMoments;
}
