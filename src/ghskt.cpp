//// [[Rcpp::depends(RcppEigen)]]
//// [[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include "Utils.h"
// #include <RcppNumerical.h>

//using namespace Numer;
using namespace arma;
using namespace Rcpp;

//// These functions are taken principally from the rugarch package of Ghalanos (2016) and have been slightly
//// modified to fit the GAS package

arma::vec paramghskt(double betabar, double nu)
{
  arma::vec param(4);
  double delta = sqrt(1/( ((2 * betabar*betabar)/((nu-2)*(nu-2)*(nu-4))) + (1/(nu-2)) ));
  double beta = betabar/delta;
  double mu = -( (beta * (delta*delta))/(nu-2));

  if (beta == 0) {
    beta = beta + 1e-12;
  }

  param(0) = nu;
  param(1) = beta;
  param(2) = delta;
  param(3) = mu;

  return param;
}

arma::mat Jacobian_paramghskt(double dBetaBar, double dNu){

  arma::vec param = paramghskt(dBetaBar, dNu);

  double dDelta = param(2);

  double dDelta_dBetaBar = -(2.0 * dBetaBar/(pow(dNu - 2.0, 2.0) * (dNu - 4.0))) * pow( 2.0 *
                                 pow(dBetaBar, 2.0)/(pow(dNu - 2.0, 2.0) * (dNu - 4.0)) + 1.0/(dNu - 2.0), -1.5);

  double dDelta_dNu  = -0.5 * pow( 2.0 * pow(dBetaBar, 2.0)/(pow(dNu - 2.0, 2.0) * (dNu - 4.0)) + 1.0/(dNu - 2.0), -1.5) *
                        ( -2.0 * pow(dBetaBar, 2.0) * (2.0*(dNu - 2.0)*(dNu - 4.0) + pow(dNu - 2.0, 2.0))/(pow(dNu - 2.0, 4.0)*pow(dNu - 4.0, 2.0))
                            - 1.0/pow(dNu - 2.0, 2.0) );

  double dBeta_dBetaBar = (dDelta - dBetaBar * dDelta_dBetaBar)/pow(dDelta, 2.0);

  double dBeta_dNu = -dBetaBar/pow(dDelta, 2.0) * dDelta_dNu;

  double dMu_dBetaBar = -(dDelta + dBetaBar * dDelta_dBetaBar)/(dNu - 2.0);

  double dMu_dNu = -dBetaBar * (dDelta_dNu * (dNu - 2.0) - dDelta)/pow(dNu - 2.0, 2.0);

  double dNu_dBetaBar = 0.0;

  double dNu_dNu = 1.0;

  arma::mat mJ(4, 2);
  mJ.zeros();

  mJ(0, 0) = dNu_dBetaBar;
  mJ(0, 1) = dNu_dNu;

  mJ(1, 0) = dBeta_dBetaBar;
  mJ(1, 1) = dBeta_dNu;

  mJ(2, 0) = dDelta_dBetaBar;
  mJ(2, 1) = dDelta_dNu;

  mJ(3, 0) = dMu_dBetaBar;
  mJ(3, 1) = dMu_dNu;

  return mJ;

}

double dghsktstd(double x, double betabar, double nu)
{
  arma::vec param = paramghskt(betabar, nu);

  double beta  = param(1);
  double delta = param(2);
  double mu    = param(3);

  double pdf = ((1.0 - nu)/2.0) * log(2.0) + nu * log(delta) + ((nu + 1.0)/2.0) * log(abs3(beta))
    + log(Rf_bessel_k(sqrt(beta*beta * (delta*delta + (x - mu)*(x - mu))), (nu + 1.0)/2.0, 2.0))
    - sqrt(beta*beta * (delta*delta + (x - mu)*(x - mu))) + beta * (x - mu) -
    Rf_lgammafn(nu/2.0) - log(M_PI)/2.0 - ((nu + 1.0)/2.0) * log(delta*delta + (x - mu)*(x - mu))/2.0;

  pdf = exp(pdf);

  return pdf;
}

double rsghst(double betabar, double nu)
{
  // Note: for values of nu<5 it is likely that sd(r) will have a very large variance
  // Existence of moment conditions (vs lower parameter bounds) are defined in the paper.
  arma::vec param = paramghskt(betabar, nu);

  double beta  = param(1);
  double delta = param(2);
  double mu    = param(3);

  double y = 1.0/Rf_rgamma(nu/2.0, 2.0/(delta*delta));
  double sigma = sqrt(y);
  double z = Rf_rnorm(0.0, 1.0);
  double ans =  mu + beta * sigma*sigma + sigma * z;

  return ans;
}


double dGHSKT(double dY, double dMuBar, double dSigma, double dBetaBar, double dNu, bool bLog = false ){
  double dZ = (dY - dMuBar)/dSigma;
  double dPDF = dghsktstd(dZ, dBetaBar, dNu)/dSigma;

  if(bLog) dPDF = log(dPDF);

  return dPDF;
}

//// [[Rcpp::depends(RcppEigen)]]
//// [[Rcpp::depends(RcppNumerical)]]
// class GHSKTPDF: public Func
// {
// private:
//   double dMu;
//   double dSigma;
//   double dBetaBar;
//   double dNu;
// public:
//   GHSKTPDF(double dMu_, double dSigma_, double dBetaBar_, double dNu_) : dMu(dMu_), dSigma(dSigma_), dBetaBar(dBetaBar_), dNu(dNu_) {}
//
//   double operator()(const double& dY) const
//   {
//     return dGHSKT(dY, dMu, dSigma, dBetaBar, dNu, false);
//   }
// };

double pGHSKT(double dY, double dMuBar, double dSigma, double dBetaBar, double dNu){

  double lower = dMuBar - 50 * dSigma, upper = dY;

  if (lower >= upper) {
    lower = upper - 50 * dSigma;
  }

  // double err_est;
  // int err_code;
  //
  // GHSKTPDF f(dMuBar, dSigma, dBetaBar, dNu);
  //
  // const double res = integrate(f, lower, upper, err_est, err_code);

  // double out = res;

  double out = NA_REAL;

  return out;

}

double rGHSKT(double dMuBar, double dSigma, double dBetaBar, double dNu){

  double dY = dMuBar + rsghst(dBetaBar, dNu)*dSigma;
  return dY;
}

double qGHSKT(double dP, double dMuBar, double dSigma, double dBetaBar, double dNu,
              int maxiter = 1e4, double eps = 1e-7) {

  double lower = dMuBar - 150 * dSigma;
  double upper = dMuBar + 150 * dSigma;

  double a=lower;
  double b=upper;

  double x=lower;
  double x1=upper;
  int iter = 1;
  double fa,fx;
  //check
  fa = pGHSKT(a, dMuBar, dSigma, dBetaBar, dNu) - dP;
  fx = pGHSKT(x1, dMuBar, dSigma, dBetaBar, dNu) - dP;

  if(fa * fx > 0){
    Rprintf("Bisection Error: upper and lower function evaluations have same sign");
    return NA_LOGICAL;
  }

  do
  {
    fa = pGHSKT(a, dMuBar, dSigma, dBetaBar, dNu) - dP;
    fx = pGHSKT(x, dMuBar, dSigma, dBetaBar, dNu) - dP;

    if (fa * fx < 0){
      b = x;
    }else{
      a = x;
    }

    x1 = (a + b)/2.0;
    iter++;

    if (abs3(x1 - x) < eps)
    {
      return x1;
    }
    x = x1;
  }while(iter<maxiter);

  Rprintf("Bisection Warning: Maximum numeber of iteration reached");
  return NA_LOGICAL;
}

double ghstexkurt(double betabar, double nu){

  double ans = 0.0;

  if(nu < 8.0){
    ans = NA_REAL;
  } else{
   arma::vec param = paramghskt(betabar, nu);

    double beta  = param(1);
    double delta = param(2);

    double beta2 = beta*beta;
    double delta2 = delta*delta;
    double k1 = 6.0/( pow(2.0 * beta2 * delta2 + (nu - 2.0) * (nu - 4.0), 2.0));
    double k21 = (nu - 2.0) * (nu - 2.0) * (nu - 4.0);
    double k22 = (16 * beta2 * delta2 * (nu - 2.0) * (nu - 4.0))/(nu - 6.0);
    double k23 = (8 * pow(beta2, 2.0) * pow(delta2, 2.0) * (5.0 * nu - 22.0))/((nu - 6.0) * (nu - 8.0));
    ans = k1 * (k21 + k22 + k23);
  }
  return ans;
}

double ghstskew(double betabar, double nu){

  double ans = 0.0;

  if(nu < 6.0){
    ans = NA_REAL;
  } else{
    arma::vec param = paramghskt(betabar, nu);
    double beta  = param(1);
    double delta = param(2);

    double beta2 = beta*beta;
    double delta2 = delta*delta;

    ans = ((2.0 * sqrt(nu - 4.0) * beta * delta)/( pow(2 * beta2 * delta2 + (nu - 2.0) * (nu - 4.0), 3/2) )) *
      (3.0 * (nu - 2.0) + ((8.0 * beta2 * delta2)/(nu - 6.0)));
  }
  return ans;
}

arma::vec mGHSKT(double dMuBar, double dSigma, double dBetaBar, double dNu){

  arma::vec vMoments(4);
  vMoments(0) = dMuBar;
  vMoments(1) = pow(dSigma, 2.0);
  vMoments(2) = ghstskew(dBetaBar, dNu);
  vMoments(3) = ghstexkurt(dBetaBar, dNu) + 3.0;

  return vMoments;
}

arma::vec ghsk_Score_OriginalParametrization(double dY, double dNu, double dBeta, double dDelta, double dMu) {

  double dK = pow(pow(dBeta, 2.0) * (pow(dDelta, 2.0) + pow(dY - dMu, 2.0)), 0.5);
  double dBessel = exp(log(Rf_bessel_k(dK, (dNu + 1.0)/2.0, 2.0)) - dK);
  double dBessel_Deriv_X = ModBesselThird_Deriv_X(dK, (dNu + 1.0)/2.0);
  double dBessel_Deriv_Nu = ModBesselThird_Deriv_Nu(dK, (dNu + 1.0)/2.0);

  double dFullDerivBessel = dBessel_Deriv_X/(dBessel * dK);

  double ddMu = -dBeta * ( dFullDerivBessel * (dBeta*(dY - dMu)) + 1.0) +
    (dNu + 1.0) * (dY - dMu)/(2.0*(pow(dDelta, 2.0) + pow(dY - dMu, 2.0)));

  double ddDelta = dNu/dDelta + pow(dBeta, 2.0) * dDelta * dFullDerivBessel -
    dDelta*(dNu + 1.0)/(2.0 * (pow(dDelta, 2.0) + pow(dY - dMu, 2.0)));

  double dSign = 1.0;
  if (dBeta < 0.0) {
    dSign = -1.0;
  }

  double ddBeta = (dNu + 1.0)/(2.0 * abs3(dBeta)) * dSign + dFullDerivBessel * dBeta * (pow(dDelta, 2.0) + pow(dY - dMu, 2.0)) + (dY - dMu);

  double ddNu = -log(2.0) * 0.5 + log(dDelta) + log(abs3(dBeta))*0.5 + 0.5 * dBessel_Deriv_Nu/dBessel -
    0.5 * Rf_digamma(0.5 * dNu) - log(pow(dDelta, 2.0) + pow(dY - dMu, 2.0))/4.0 ;

  double ddY = dBeta * (1.0 + dBeta * (dY - dMu) * dFullDerivBessel) - 0.5*(dNu + 1.0) * (dY - dMu)/(pow(dDelta, 2.0) + pow(dY - dMu, 2.0));

  arma::vec vScore_Original(5);

  vScore_Original.zeros();

  vScore_Original(0) = ddNu;
  vScore_Original(1) = ddBeta;
  vScore_Original(2) = ddDelta;
  vScore_Original(3) = ddMu;
  vScore_Original(4) = ddY;

  return vScore_Original;
}

arma::vec ghskt_Score(double dY, arma::vec vTheta) {

  double dMuBar = vTheta(0);
  double dSigma = vTheta(1);
  double dBetaBar = vTheta(2);
  double dNu      = vTheta(3);

  arma::vec param = paramghskt(dBetaBar, dNu);

  double dBeta  = param(1);
  double dDelta = param(2);
  double dMu    = param(3);

  arma::mat mJ = Jacobian_paramghskt(dBetaBar, dNu);

  double dZ = (dY - dMuBar)/dSigma;

  arma::vec vScore_Full = ghsk_Score_OriginalParametrization(dZ, dNu, dBeta, dDelta, dMu);

  arma::vec vScore_Original = vScore_Full.subvec(0, 3);

  arma::vec vScore_BetaBar_and_Nu = mJ.t() * vScore_Original;

  double ddZ = vScore_Full(4);

  double ddMuBar = -ddZ/dSigma;
  double ddSigma = -1.0/(dSigma) * (1.0 + ddZ * dZ);

  arma::vec vScore(4);

  vScore(0) = ddMuBar;
  vScore(1) = ddSigma;
  vScore(2) = vScore_BetaBar_and_Nu(0);
  vScore(3) = vScore_BetaBar_and_Nu(1);

  return vScore;

}

arma::mat ghskt_IM(arma::vec vTheta){

  arma::mat mIM = eye(4, 4);

  return mIM;

}

