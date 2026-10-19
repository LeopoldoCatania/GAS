#include <RcppArmadillo.h>
#include "ber.h"
#include "Utils.h"
#include "TwoPiece.h"
#include "TwoPieceInformation.h"

using namespace arma;
using namespace Rcpp;

//// These functions are taken principally from the rugarch package of Ghalanos (2016) and have been slightly
//// modified to fit the GAS package

double dnormstd(const double x) { return Rf_dnorm4(x, 0, 1, 0); }

double rsnorm(const double xi)
{
  double weight, z, rr, m1, mu, sigma, xx, ans;
  weight = xi / (xi + 1.0/xi);
  z =  Rf_runif(-weight, 1.0 - weight);
  xx = (z < 0)? 1.0/xi : xi;
  rr = -1.0 * abs3(Rf_rnorm(0.0, 1.0))/xx * sign_C(z);
  m1 = 2.0/sqrt(2.0 * M_PI);
  mu = m1 * (xi - 1.0/xi);
  sigma = sqrt((1 - (m1 * m1)) * ( (xi * xi) + 1.0/(xi* xi) ) + 2 * (m1 * m1) - 1.0);
  ans = (rr - mu ) / sigma;
  return(ans);
}


double GAS_log_snorm(double x, double xi) {
  if (ISNAN(x) || !R_FINITE(xi) || xi <= 0) return NA_REAL;
  const GAS_TwoPiece t(xi, sqrt(2.0/M_PI));
  const double z = x*t.sd+t.mean;
  const double u = z < 0 ? z*xi : z/xi;
  return log(2.0)-log(xi+1.0/xi)+log(t.sd)+Rf_dnorm4(u,0,1,1);
}
double dsnormstd(const double x, const double xi) { return exp(GAS_log_snorm(x, xi)); }

double psnorm(const double q, const double mu, const double sigma, const double xi) {
  if (ISNAN(q) || !R_FINITE(mu) || !R_FINITE(sigma) || sigma <= 0 ||
      !R_FINITE(xi) || xi <= 0) return NA_REAL;
  const GAS_TwoPiece t(xi, sqrt(2.0/M_PI));
  const double z = (q-mu)/sigma*t.sd + t.mean;
  if (z <= 0) return exp(log(2.0)+t.logleft+Rf_pnorm5(z*xi,0,1,1,1));
  return -expm1(log(2.0)+t.logright+Rf_pnorm5(-z/xi,0,1,1,1));
}

double qsnorm(const double p, const double xi) {
  if (!R_FINITE(p) || p < 0 || p > 1 || !R_FINITE(xi) || xi <= 0)
    return NA_REAL;
  if (p == 0) return R_NegInf;
  if (p == 1) return R_PosInf;
  const GAS_TwoPiece t(xi, sqrt(2.0/M_PI));
  double z;
  if (p == exp(t.logleft)) z = 0.0;
  else if (log(p) < t.logleft) {
    const double lp = log(p)-log(2.0)-t.logleft;
    z = (Rf_qnorm5(lp,0,1,1,1))/xi;
  } else {
    const double lp = log1p(-p)-log(2.0)-t.logright;
    z = -xi*(Rf_qnorm5(lp,0,1,1,1));
  }
  return (z-t.mean)/t.sd;
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

double dSNORM(double dY, double dMu, double dSigma, double dXi, bool bLog=false) {
  if (!R_FINITE(dMu) || !R_FINITE(dSigma) || dSigma <= 0) return NA_REAL;
  const double ld = GAS_log_snorm((dY-dMu)/dSigma, dXi)-log(dSigma);
  return bLog ? ld : exp(ld);
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
  if (!R_FINITE(dMu) || !R_FINITE(dSigma) || dSigma <= 0) return NA_REAL;
  return dMu+dSigma*qsnorm(dP,dXi);
}

arma::vec mSNORM(double dMu, double dSigma, double dXi) {

  const double c = sqrt(2.0 / M_PI);
  const double c2 = c * c;

  const double delta = dXi - 1.0 / dXi;
  const double delta2 = delta * delta;

  // Variance and central moments of the unstandardized
  // two-piece normal variable.
  const double variance =
    1.0 + (1.0 - c2) * delta2;

  const double central3 =
    c * delta * (1.0 + (2.0 * c2 - 1.0) * delta2);

  const double central4 =
    3.0 +
    (9.0 - 10.0 * c2) * delta2 +
    (3.0 - 2.0 * c2 - 3.0 * c2 * c2) * delta2 * delta2;

  arma::vec vMoments(4);

  // The package standardizes the two-piece variable before
  // applying dMu and dSigma.
  vMoments(0) = dMu;
  vMoments(1) = dSigma * dSigma;
  vMoments(2) = central3 / pow(variance, 1.5);
  vMoments(3) = central4 / (variance * variance);

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

  double m1 = 2.0/sqrt(2.0*M_PI);
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
  double ddSigma_tilde_xi = (1.0 - m12) * (dXi - pow(dXi, -3.0)) / dSigma_tilde;
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
  return GAS_two_piece_information(vTheta, false);
}
