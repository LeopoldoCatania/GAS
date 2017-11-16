#include <RcppArmadillo.h>
#include "ber.h"
#include "Utils.h"

using namespace arma;
using namespace Rcpp;

//// These functions are taken principally from the rugarch package of Ghalanos (2016) and have been slightly
//// modified to fit the GAS package

double dnormstd(const double x)
{
  double pdf;
  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * PI );
  if(pdf == 0.0) pdf = 0.0 + 2.22507e-24;
  return pdf;
}

double rsnorm(const double xi)
{
  double weight, z, rr, m1, mu, sigma, xx, ans;
  weight = xi / (xi + 1.0/xi);
  z =  Rf_runif(-weight, 1.0 - weight);
  xx = (z < 0)? 1.0/xi : xi;
  rr = -1.0 * abs3(Rf_rnorm(0.0, 1.0))/xx * sign_C(z);
  m1 = 2.0/sqrt(2.0 * PI);
  mu = m1 * (xi - 1.0/xi);
  sigma = sqrt((1 - (m1 * m1)) * ( (xi * xi) + 1.0/(xi* xi) ) + 2 * (m1 * m1) - 1.0);
  ans = (rr - mu ) / sigma;
  return(ans);
}

double dsnormstd(const double x, const double xi)
{
  double pdf;
  double mu, sigma,z, xxi, g;
  double m1 = 2.0/sqrt(2.0*PI);
  double m12 = m1*m1;
  double xi2 = xi*xi;
  mu = m1*(xi-1.0/xi);
  sigma = sqrt((1.0-m12)*(xi2+1.0/xi2)+2.0*m12-1.0);
  z = x*sigma+mu;
  xxi = (z<0)? 1.0/xi : xi;
  g = 2.0/(xi + 1.0/xi);
  pdf = g * dnormstd(z/xxi)*sigma;
  return pdf;
}

double psnorm(const double q, const double mu, const double sigma, const double xi)
{
  double qx = (q-mu)/sigma;
  double m1 = 2.0/sqrt(2*PI);
  double mux = m1 * (xi - 1.0/xi);
  double sig = sqrt((1.0-m1*m1)*(xi*xi+1.0/(xi*xi)) + 2.0*m1*m1 - 1.0);
  double z = qx*sig + mux;
  double Xi = (z<0)?1.0/xi:xi;
  double g = 2.0/(xi + 1.0/xi);
  double p = Heaviside(z, 0) - signum(z) * g * Xi * Rf_pnorm5(-abs3(z)/Xi, 0, 1, 1, 0);
  return( p );
}

double qsnorm(const double p, const double xi)
{
  double m1 = 2.0/sqrt(2*PI);
  double mu = m1 * (xi - 1.0/xi);
  double sigma = sqrt((1.0-m1*m1)*(xi*xi+1.0/(xi*xi)) + 2.0*m1*m1 - 1.0);
  double g = 2.0/(xi + 1.0/xi);
  double z = p-0.5;
  double Xi = (z<0)?1.0/xi:xi;
  double tmp = (Heaviside(z, 0) - signum(z) * p)/ (g* Xi);
  double q = (-1.0*signum(z)* Rf_qnorm5(tmp, 0, Xi, 1, 0) - mu)/sigma;
  return( q );
}


double snormskew( double dXi )
{
  double m1 = 2.0/sqrt(2.0 * M_PI);
  double m2 = 1.0;
  double m3 = 4/sqrt(2.0 * M_PI);

  double dSkew = (dXi - 1.0/dXi) * ( ( m3 + 2.0 * pow(m1, 3.0) - 3.0 * m1 * m2 ) * ( pow(dXi,2.0) + (1.0/pow(dXi, 2.0) ) ) + 3.0 * m1 * m2 - 4.0 * pow(m1,3.0) )/
    ( pow(( (m2 - pow(m1, 2.0)) * ( pow(dXi,2.0) + 1.0/pow(dXi, 2.0) ) + 2.0 * pow(m1, 2.0) - m2), 1.5)  );

  return dSkew;
}

double dSNORM(double dY, double dMu, double dSigma, double dXi, bool bLog = false){

  double dZ = (dY - dMu)/dSigma;
  double dPDF = dsnormstd(dZ, dXi)/dSigma;

  if(bLog) dPDF = log(dPDF);

  return dPDF;

}

double pSNORM(double dY, double dMu, double dSigma, double dXi){

  double dP = psnorm(dY, dMu, dSigma, dXi);
  return dP;

}

double rSNORM(double dMu, double dSigma, double dXi){

  double dY = dMu + rsnorm(dXi)*dSigma;
  return dY;

}

double qSNORM(double dP, double dMu, double dSigma, double dXi) {

  double dQ = qsnorm(dP, dXi)*dSigma + dMu;

  return dQ;
}

arma::vec mSNORM(double dMu, double dSigma, double dXi){

  arma::vec vMoments(4);

  vMoments(0) = dMu;
  vMoments(1) = pow(dSigma, 2.0);
  vMoments(2) = snormskew(dXi);
  vMoments(3) = 0.0;

  return vMoments;

}

arma::vec snorm_Score(double dY, arma::vec vTheta){

  double dMu    = vTheta(0);
  double dSigma = vTheta(1);
  double dXi     = vTheta(2);

  arma::vec vScore(3);

  double dK = (dY - dMu)/dSigma;

  double ddK_mu = -1.0/dSigma;
  double ddK_sigma = -dK/dSigma;

  double m1 = 2.0/sqrt(2.0*PI);
  double m12 = m1*m1;
  double xi2 = dXi * dXi;

  double dMu_tilde = m1 * (dXi - 1.0/dXi);
  double dSigma_tilde = sqrt((1.0-m12)*(xi2+1.0/xi2)+2.0*m12-1.0);

  double dZ = dK * dSigma_tilde + dMu_tilde;

  double dXi_star = (dZ<0)? 1.0/dXi : dXi;

  //derivative log dsnormstd wrt dK
  double dA = - dZ * dSigma_tilde / pow(dXi_star, 2.0);

  double ddMu = ddK_mu * dA;
  double ddSigma = ddK_sigma * dA - 1.0/dSigma;

  // derivative wrt dXi

  double dG = 2.0/(dXi + 1.0/dXi);
  double dH = dZ/dXi_star;

  double ddG_xi = -2.0*(1.0 - 1.0/pow(dXi, 2.0))/pow(dXi + 1.0/dXi, 2.0);
  double ddXi_Star_xi = (dZ<0)? -1.0/pow(dXi, 2.0) : 1.0;
  double ddSigma_tilde_xi = -(1.0 - m12)*(dXi - 1.0/pow(dXi, 3.0))/dSigma_tilde;
  double ddMu_tilde_xi = m1 * (1.0 + 1.0/pow(dXi, 2.0));

  double ddZ_xi = dK * ddSigma_tilde_xi + ddMu_tilde_xi;

  double ddH_xi = (ddZ_xi * dXi_star - dZ * ddXi_Star_xi )/pow(dXi_star, 2.0);

  double ddXi = ddG_xi/dG - dH*ddH_xi +  ddSigma_tilde_xi/dSigma_tilde;

  vScore(0) = ddMu;
  vScore(1) = ddSigma;
  vScore(2) = ddXi;

  return vScore;
}

arma::mat snorm_IM(arma::vec vTheta){

  arma::mat mIM = eye(3,3);

  return mIM;

}
