#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double dALD(double dY, double dTheta, double dSigma, double dKappa, bool bLog = false){

  double dLPDF = 0.5 * log(2.0) - log(dSigma) + log(dKappa) - log(1.0 + pow(dKappa, 2.0));

  if(dY<dTheta){
    dLPDF += pow(2.0,0.5)*(dY - dTheta)/(dSigma * dKappa);
  }else{
    dLPDF += - dKappa * pow(2.0,0.5)*(dY - dTheta)/(dSigma);
  }

  if(!bLog) dLPDF = exp(dLPDF);

  return dLPDF;
}
double pALD(double dY, double dTheta, double dSigma, double dKappa){

  double dP = 0.0;

  if(dY<dTheta){

    dP = pow(dKappa, 2.0)/(1.0 + pow(dKappa,2.0)) * exp(-pow(2.0 , 0.5)* abs3(dY - dTheta)/(dSigma * dKappa));

  }else{

    dP = 1.0 - 1.0/(1.0 + pow(dKappa,2.0)) * exp(-pow(2.0 , 0.5)* dKappa * abs3(dY - dTheta)/(dSigma));

  }

  return dP;

}
double rALD(double dTheta, double dSigma, double dKappa){

  double dU1 = Rf_runif(0.0, 1.0);
  double dU2 = Rf_runif(0.0, 1.0);

  double dY = dTheta + dSigma * ( dKappa*log(dU1) - log(dU2)/dKappa );

  return dY;

}
double qALD(double dP, double dTheta, double dSigma, double dKappa,
            double lower=-150, double upper=150, int maxiter=1e4, double eps=1e-7) {
  double a=lower;
  double b=upper;

  double x=lower;
  double x1=upper;
  int iter = 1;
  double fa,fx;
  //check
  fa = pALD(a,  dTheta, dSigma, dKappa) - dP;
  fx = pALD(x1, dTheta, dSigma, dKappa) - dP;

  if(fa*fx>0){
    Rprintf("Bisection Error: upper and lower function evaluations have same sign");
    return NA_LOGICAL;
  }

  do
  {
    fa=pALD(a, dTheta, dSigma, dKappa) - dP;
    fx=pALD(x, dTheta, dSigma, dKappa) - dP;

    if (fa*fx < 0){
      b=x;
    }else{
      a=x;
    }

    x1=(a+b)/2.0;
    iter++;

    if (abs3(x1-x) < eps)
    {
      return x1;
    }
    x=x1;
  }while(iter<maxiter);

  Rprintf("Bisection Warning: Maximum numeber of iteration reached");
  return NA_LOGICAL;
}
arma::vec mALD(double dTheta, double dSigma, double dKappa){

  arma::vec vMoments(4);

  vMoments(0) = dTheta + dSigma*(1.0/dKappa - dKappa)/pow(2.0, 0.5);
  vMoments(1) = pow(dSigma,2.0)*(1.0/pow(dKappa,2.0) + pow(dKappa,2.0))/2.0;
  vMoments(2) = 2.0*(1.0/pow(dKappa, 3.0) - pow(dKappa, 3.0))/pow(1.0/pow(dKappa,2.0) + pow(dKappa, 2.0), 1.5);
  vMoments(3) = 6.0 - 12.0/pow( 1.0/pow(dKappa, 2.0) + pow(dKappa, 2.0)  , 2.0);

  return vMoments;
}
arma::vec ald_Score(double dY, arma::vec vTheta){

  arma::vec vScore(3);

  double dTheta = vTheta(0);
  double dSigma = vTheta(1);
  double dKappa = vTheta(2);

  double dS2 = pow(2.0, 0.5);

  if(dY < dTheta){
    vScore(0) = -dS2/(dSigma*dKappa);
    vScore(1) = -1.0/dSigma - dS2*(dY - dTheta)/(pow(dSigma, 2.0)*dKappa);
    vScore(2) = 1.0/dKappa - 2.0*dKappa/(1.0 + pow(dKappa, 2.0)) - dS2*(dY-dTheta)/(dSigma*pow(dKappa, 2.0));
  }else{
    vScore(0) = dKappa*dS2/(dSigma);
    vScore(1) = -1.0/dSigma + dS2*dKappa*(dY - dTheta)/(pow(dSigma, 2.0));
    vScore(2) = 1.0/dKappa - 2.0*dKappa/(1.0 + pow(dKappa, 2.0)) - dS2*(dY-dTheta)/(dSigma);
  }
  return vScore;
}
arma::mat ald_IM(arma::vec vTheta){

  arma::mat mIM=zeros(3,3);

  double dSigma = vTheta(1);
  double dKappa = vTheta(2);

  double uu = 2.0/pow(dSigma,2.0);
  double du = -2.0*pow(2.0, 0.5)/(dSigma*(1.0 + pow(dKappa, 2.0)));
  double dd = 1.0/pow(dKappa, 2.0) + 4.0/pow(1.0 + pow(dKappa,2.0),2.0);
  double td = -(1.0 - pow(dKappa, 2.0))/(dSigma*dKappa*(1.0 + pow(dKappa, 2.0)));
  double tt = 1.0/pow(dSigma, 2.0);

  mIM(0,0) = uu;

  mIM(2,0) = du;
  mIM(0,2) = du;

  mIM(1,1) = tt;

  mIM(2,1) = td;
  mIM(1,2) = td;

  mIM(2,2) = dd;

  return mIM;

}


