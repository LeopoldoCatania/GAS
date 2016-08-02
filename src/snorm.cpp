#include <RcppArmadillo.h>
#include "ber.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

double dSNORM(double dY, double dMu, double dSigma2, double dDelta, bool bLog = false){

  double dLPDF = log(2.0*dDelta) - log(1.0 + pow(dDelta, 2.0)) - 0.5*log(dSigma2);

  double dZ = 0.0;

  if (dY < dMu){
    dZ = (dDelta*(dY - dMu))/pow(dSigma2, 0.5);
  }else{
    dZ = (dY - dMu)/(dDelta * pow(dSigma2, 0.5));
  }

  dLPDF += Rf_dnorm4(dZ, 0.0, 1.0, 1);

  if(!bLog) dLPDF = exp(dLPDF);

  return dLPDF;

}


double pSNORM(double dY, double dMu, double dSigma2, double dDelta){

  double dP = 0.0;

  double dZ = 0.0;

  if (dY < dMu){

    dZ = (dDelta*(dY - dMu))/pow(dSigma2, 0.5);

    dP = 2.0/(1.0 + pow(dDelta, 2.0)) * Rf_pnorm5(dZ, 0.0,1.0,1,0);

  }else{

    dZ = (dY - dMu)/(dDelta * pow(dSigma2, 0.5));

    dP = (1.0 - pow(dDelta, 2.0))/ (1.0 + pow(dDelta, 2.0)) +  2.0 * pow(dDelta, 2.0)/(1.0 + pow(dDelta, 2.0)) * Rf_pnorm5(dZ, 0.0,1.0,1,0);

  }

  return dP;

}


double rSNORM(double dMu, double dSigma2, double dDelta){

  double dZ0 = Rf_rnorm(0.0, 1.0);
  double dZ1 = Rf_rnorm(0.0, 1.0);

  double dPi = 1.0/(1.0 + pow(dDelta, 2.0));

  double dU = rBER(dPi);

  double dY = dMu + pow(dSigma2, 0.5) * (dDelta*(1.0 - dU)*abs3(dZ0) - dU*abs3(dZ1)/dDelta);

  return dY;

}


double qSNORM(double dP, double dMu, double dSigma2, double dDelta,
            double lower=-150, double upper=150, int maxiter=1e4, double eps=1e-7) {
  double a=lower;
  double b=upper;

  double x=lower;
  double x1=upper;
  int iter = 1;
  double fa,fx;
  //check
  fa = pSNORM(a,  dMu, dSigma2, dDelta) - dP;
  fx = pSNORM(x1, dMu, dSigma2, dDelta) - dP;

  if(fa*fx>0){
    Rprintf("Bisection Error: upper and lower function evaluations have same sign");
    return NA_LOGICAL;
  }

  do
  {
    fa=pSNORM(a, dMu, dSigma2, dDelta) - dP;
    fx=pSNORM(x, dMu, dSigma2, dDelta) - dP;

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


arma::vec mSNORM(double dMu, double dSigma2, double dDelta){
  arma::vec vMoments(4);

  double dSigma = pow(dSigma2, 0.5);
  double dDelta2 = pow(dDelta, 2.0);

  vMoments(0) = dMu + pow(2.0, 0.5)*dSigma*(dDelta2 - 1.0)/(pow(M_PI, 0.5)*dDelta);
  vMoments(1) = (dSigma2*((M_PI - 2.0)*pow(dDelta, 6.0) + 2.0*dDelta2*(dDelta2 + 1.0) + M_PI - 2.0))/(M_PI * dDelta2*(1.0 + dDelta2));
  vMoments(2) = pow(2.0, 0.5)*(1.0 - dDelta2)*pow(1.0 + dDelta2, 0.5)*((M_PI - 4.0)*pow(dDelta, 6.0) + 2.0*(2.0 - M_PI)*dDelta2*(1.0 + dDelta2) + M_PI - 4.0)/pow((M_PI - 2.0)*pow(dDelta, 6.0) + 2.0*dDelta2*(1.0 + dDelta2) + (M_PI - 2.0), 3.0/2.0);
  vMoments(3) = ( (1.0 + dDelta2)*((3.0*pow(M_PI, 2.0)- 4.0*M_PI - 12.0)*(pow(dDelta, 10.0) + 1.0) + 4.0*(9.0 - 2.0*M_PI)*dDelta2*(pow(dDelta, 6.0) + 1.0)) )/pow( (M_PI - 2.0)*pow(dDelta, 6.0) + 2.0*dDelta2*(1.0 + dDelta2) + (M_PI - 2.0), 2.0) +
                (1.0 + dDelta2)*(12.0*(M_PI - 2.0)*pow(dDelta, 4.0)*(dDelta2 + 1.0))/pow((M_PI - 2.0)*pow(dDelta, 6.0) + 2.0*dDelta2*(1.0 + dDelta2) + (M_PI - 2.0), 2.0);

  return vMoments;

}

arma::vec snorm_Score(double dY, arma::vec vTheta){

  double dMu     = vTheta(0);
  double dSigma2 = vTheta(1);
  double dDelta  = vTheta(2);

  double dDelta2 = pow(dDelta, 2.0);

  arma::vec vScore(3);

  double dI = 1.0;
  if(dY >= dMu) dI = 0.0;

  double dI_c = 1.0 - dI;

  double ddMu     = dDelta2/dSigma2 * (dY - dMu)*dI + dI_c*(dY - dMu)/(dDelta2*dSigma2);
  double ddSigma2 = -1.0/(2.0*dSigma2) + dDelta2*pow(dY - dMu, 2.0)/(2.0*pow(dSigma2, 2.0)) * dI + dI_c*pow(dY - dMu, 2.0)/(2.0*dDelta2*pow(dSigma2, 2.0));
  double ddDelta  = 1.0/dDelta - 2.0*dDelta/(1.0 + dDelta2) - dDelta/dSigma2 * pow(dY - dMu, 2.0)*dI + dI_c * pow(dY - dMu, 2.0)/(pow(dDelta, 3.0)*dSigma2);

  vScore(0) = ddMu;
  vScore(1) = ddSigma2;
  vScore(2) = ddDelta;

  return vScore;
}

arma::mat snorm_IM(arma::vec vTheta){

  double dMu     = vTheta(0);
  double dSigma2 = vTheta(1);
  double dDelta  = vTheta(2);
  double dDelta2 = pow(dDelta, 2.0);

  arma::mat mIM = zeros(3,3);

  double uu = 1.0/dSigma2;
  double tu = 8.0/(pow(2.0*M_PI*dSigma2, 0.5)*(1.0 + dDelta2));
  double dd = 1.0/(2.0*pow(dSigma2, 2.0));
  double td = (dDelta2 - 1.0)/(dDelta*(1.0 + dDelta2)*dSigma2);
  double tt = 2.0/dDelta2 + 4.0/pow(1.0 + dDelta2 , 2.0);

  mIM(0,0) = uu;
  mIM(2,0) = tu;
  mIM(0,2) = tu;
  mIM(1,1) = dd;
  mIM(2,1) = td;
  mIM(1,2) = td;
  mIM(2,2) = tt;

  return mIM;

}
