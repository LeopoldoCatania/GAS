#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double dSKELLAM(double dY, double dMu, double dSigma2, bool bLog = false) {

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double dB = ModBesselFirst(2.0 * pow(dMu1 * dMu2, 0.5), abs3(dY));

  double dPDF = -(dMu1 + dMu2) + 0.5 * dY * log(dMu1/dMu2) + log(dB);

  if (!bLog) {
    dPDF = exp(dPDF);
  }

  return dPDF;
}


// This code has been taken from the R code of the skellam package of
// Jerry W Lewis, Patrick E Brown and Michail Tsagris (see
// see http://r-forge.r-project.org/projects/healthqueues) and slightly modified to fit the GAS package.
double pSKELLAM(double dY, double dMu, double dSigma2) {

  // double xm = -floor(dY) - 0.5;
  // double s  = log(0.5 * (xm + pow(pow(xm, 2.0) + 4.0 * dMu1 * dMu2, 0.5))/dMu2);
  // double K  = dMu2 * (exp(s) - 1.0) + dMu1 * (exp(-s) - 1.0);
  // double K2 = dMu2 * exp(s) + dMu1 * exp(-s);
  // double u2 = 2.0 * sinh(0.5 * s) * pow(K2, 0.5);
  // double w2 = pow(2 * (s * xm - K), 0.5);
  // if (s < 0) {
  //   w2 = -w2;
  // }
  // double xe = (xm + (dMu1 - dMu2))/pow(dMu1 + dMu2, 0.5);
  // double g1 = (dMu1 - dMu2)/pow(dMu1 + dMu2, 1.5);
  //
  // double dCDF = 0.0;
  //
  // if (abs3(xe) < 1e-4) {
  //   dCDF = Rf_pnorm5(-xe, 0.0, 1.0, 1, 0) + Rf_dnorm4(xe, 0.0, 1.0, 0) * g1/6.0 * (1.0 - pow(xe, 2.0));
  // } else {
  //   dCDF = Rf_pnorm5(-w2, 0.0, 1.0, 1, 0) - Rf_dnorm4(w2, 0.0, 1.0, 0) * (1.0/w2 - 1.0/u2);
  // }

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double dCDF = 0.0;

  if (dY < 0) {
    dCDF = Rf_pnchisq(2.0 * dMu2, -2 * dY, 2.0 * dMu1, 1, 0);
  } else {
    dCDF = Rf_pnchisq(2.0 * dMu1, 2 * (dY + 1.0), 2.0 * dMu2, 0, 0);
  }

  return dCDF;

}

// This code has been taken from the R code of the skellam package of
// Jerry W Lewis, Patrick E Brown and Michail Tsagris (see
// see http://r-forge.r-project.org/projects/healthqueues) and
// slightly modified to fit the GAS package.
double qSKELLAM(double dP, double dMu, double dSigma2) {

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double p = dP * (1 - 64 * 2.22e-16);

  double z = Rf_qnorm5(dP, 0.0, 1.0, 1, 0);
  double mu = dMu1 - dMu2;
  double vr = dMu1 + dMu2;

  double sg = pow(vr, 0.5);
  double c0 = mu + z * sg;
  double c1 = (pow(z, 2.0) - 1.0) * mu/vr/6.0;
  double c2 = -(c1 * mu - 2.0 * dMu1 * dMu2 * (pow(z, 2.0) - 3.0)/vr) *  z/12.0/vr/sg;

  double q0 = round(c0 + c1 + c2);
  double p0 = pSKELLAM(q0, dMu1, dMu2);

  if (p0 < p) {
  while (p0 < p) {
    q0 += 1;
    p0 = pSKELLAM(q0, dMu1, dMu2);
    }
  } else {

    p0 = pSKELLAM(q0 - 1, dMu1, dMu2);

    while (p0 > p) {
      q0 -= 1;
      p0 = pSKELLAM(q0 - 1, dMu1, dMu2);
    }
  }

  return q0;

}

double rSKELLAM(double dMu, double dSigma2) {

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double dY1 = Rf_rpois(dMu1);
  double dY2 = Rf_rpois(dMu2);

  double dY = dY1 - dY2;

  return dY;

}

arma::vec mSKELLAM(double dMu, double dSigma2) {

  arma::vec vMoments(4);

  vMoments(0) = dMu;
  vMoments(1) = dSigma2;
  vMoments(2) = dMu/pow(dSigma2, 1.5);
  vMoments(3) = 1/dSigma2 + 3.0;

  return vMoments;

}

arma::vec skellam_Score(double dY, arma::vec vTheta) {



  double dMu = vTheta(0);
  double dSigma2 = vTheta(1);

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double dDerivBessel = ModBesselFirst_Deriv(2.0 * pow(dMu1 * dMu2, 0.5), abs3(dY));

  double dBessel = ModBesselFirst(2.0 * pow(dMu1 * dMu2, 0.5), abs3(dY));

  double dMu1_s = 0.5 * dY / dMu1 - 1.0 + pow(dMu2/dMu1, 0.5) * dDerivBessel / dBessel;
  double dMu2_s = - 0.5 * dY / dMu2 - 1.0 + pow(dMu1/dMu2, 0.5) * dDerivBessel / dBessel;

  arma::vec vScore(2);

  vScore(0) = dMu1_s;
  vScore(1) = dMu2_s;

  arma::mat mJ(2, 2);

  mJ(0, 0) = 0.5;
  mJ(0, 1) = 0.5;
  mJ(1, 0) = -0.5;
  mJ(1, 1) = 0.5;

  return mJ.t() * vScore;

}


arma::mat skellam_IM(arma::vec vTheta) {
  arma::mat mI = eye(2, 2);
  return mI;
}

