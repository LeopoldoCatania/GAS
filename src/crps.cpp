#include <RcppArmadillo.h>
#include "DistWrap.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

double WeightFun(double dZ, std::string sType, double dA = 0.0, double dB = 1.0){
  double dW = 0.0;

  if(sType=="uniform")  dW = 1.0;
  if(sType=="center")   dW = Rf_dnorm4(dZ, dA, dB, 0.0);
  if(sType=="tails")    dW = (1.0 - Rf_dnorm4(dZ, dA, dB, 0.0)/Rf_dnorm4(0.0, dA, dB, 0.0));
  if(sType=="tail_r")  dW = Rf_pnorm5(dZ, dA, dB, 1.0, 0.0);
  if(sType=="tail_l")  dW = (1.0 - Rf_pnorm5(dZ, dA, dB, 1.0, 0.0));

  return dW;
}

double wCRPS_internal(double dY, arma::vec vTheta, std::string Dist, double dLower, double dUpper,
             std::string sType = "uniform", int iB = 1000,  double dA = 0.0, double dB = 1.0){

  double dY_i = 0.0;

  double dWCRPS = 0.0;
  double dF     = 0.0;
  double dInd   = 0.0;
  double dW     = 0.0;

  double dBrierPS = 0.0;

  for(int b = 0; b<iB; b++){
    dY_i = dLower + b * (dUpper - dLower)/ iB;

    dF = pdist_univ(dY_i, vTheta, Dist);
    dInd = IndicatorLess(dY, dY_i);

    dBrierPS = pow(dF - dInd, 2.0);
    dW       = WeightFun(dY_i, sType, dA, dB);

    dWCRPS += dW * dBrierPS;
  }

  dWCRPS *= (dUpper - dLower)/(iB - 1.0);

  return dWCRPS;
}

arma::vec wCRPS_series(arma::vec vY, arma::mat mTheta, std::string Dist, double dLower, double dUpper,
                            std::string sType = "uniform", int iB = 1000,  double dA = 0.0, double dB = 1.0){

  int iT = vY.size();
  arma::vec vWCRPS(iT);

  for(int i=0;i<iT;i++){
    vWCRPS(i) = wCRPS_internal(vY(i), mTheta.col(i), Dist, dLower, dUpper, sType, iB, dA, dB);
  }

  return vWCRPS;

}

//[[Rcpp::export]]
arma::mat mWCRPS_backtest(arma::vec vY, arma::mat mTheta, std::string Dist, double dLower, double dUpper,
                          int iB = 1000,  double dA = 0.0, double dB = 1.0){
  int iT = vY.size();
  arma::mat mWCRPS(iT, 5);

  mWCRPS.col(0) = wCRPS_series(vY, mTheta, Dist, dLower, dUpper, "uniform", iB,  dA, dB);
  mWCRPS.col(1) = wCRPS_series(vY, mTheta, Dist, dLower, dUpper, "center", iB,  dA, dB);
  mWCRPS.col(2) = wCRPS_series(vY, mTheta, Dist, dLower, dUpper, "tails", iB,  dA, dB);
  mWCRPS.col(3) = wCRPS_series(vY, mTheta, Dist, dLower, dUpper, "tail_r", iB,  dA, dB);
  mWCRPS.col(4) = wCRPS_series(vY, mTheta, Dist, dLower, dUpper, "tail_l", iB,  dA, dB);

  return mWCRPS;

}


