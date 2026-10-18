#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

static bool validSkellam(double mu, double v) {
  return R_FINITE(mu) && R_FINITE(v) && v > fabs(mu);
}

double dSKELLAM(double dY, double dMu, double dSigma2, bool bLog=false) {
  if (!validSkellam(dMu,dSigma2) || ISNAN(dY)) return NA_REAL;
  if (!R_FINITE(dY) || dY != floor(dY)) return bLog ? R_NegInf : 0.0;
  const double l1 = 0.5*(dSigma2+dMu), l2 = 0.5*(dSigma2-dMu);
  const double a = sqrt(l1), b = sqrt(l2), z = 2*a*b;
  const double scaledI = Rf_bessel_i(z,fabs(dY),2);
  const double ld = -(a-b)*(a-b)+0.5*dY*(log(l1)-log(l2))+log(scaledI);
  return bLog ? ld : exp(ld);
}


// This code has been taken from the R code of the skellam package of
// Jerry W Lewis, Patrick E Brown and Michail Tsagris (see
// see https://r-forge.r-project.org/projects/healthqueues) and slightly modified to fit the GAS package.
double pSKELLAM(double dY, double dMu, double dSigma2) {
  if (!validSkellam(dMu,dSigma2) || ISNAN(dY)) return NA_REAL;
  if (dY == R_NegInf) return 0.0;
  if (dY == R_PosInf) return 1.0;
  const double y=floor(dY);
  const double l1=0.5*(dSigma2+dMu), l2=0.5*(dSigma2-dMu);
  return y < 0 ? Rf_pnchisq(2*l2,-2*y,2*l1,1,0)
               : Rf_pnchisq(2*l1,2*(y+1),2*l2,0,0);
}

// This code has been taken from the R code of the skellam package of
// Jerry W Lewis, Patrick E Brown and Michail Tsagris (see
// see https://r-forge.r-project.org/projects/healthqueues) and
// slightly modified to fit the GAS package.
double qSKELLAM(double dP, double dMu, double dSigma2) {
  if (!validSkellam(dMu,dSigma2) || !R_FINITE(dP) || dP < 0 || dP > 1)
    return NA_REAL;
  if (dP == 0) return R_NegInf;
  if (dP == 1) return R_PosInf;
  // Limit the search to exactly representable integer arithmetic.
  const double limit=4503599627370496.0;
  double q=floor(dMu+sqrt(dSigma2)*Rf_qnorm5(dP,0,1,1,0));
  if (!R_FINITE(q) || fabs(q)>=limit) return NA_REAL;
  double f=pSKELLAM(q,dMu,dSigma2), lo=q, hi=q, step=1.0;
  if (!R_FINITE(f)) return NA_REAL;
  bool bracketed=false;
  for (int j=0;j<64;++j) {
    if (f < dP) {
      hi=q+step;
      if (fabs(hi)>=limit) return NA_REAL;
      const double fh=pSKELLAM(hi,dMu,dSigma2);
      if (!R_FINITE(fh)) return NA_REAL;
      if (fh>=dP) { bracketed=true; break; }
      lo=hi;
    } else {
      lo=q-step;
      if (fabs(lo)>=limit) return NA_REAL;
      const double fl=pSKELLAM(lo,dMu,dSigma2);
      if (!R_FINITE(fl)) return NA_REAL;
      if (fl<dP) { bracketed=true; break; }
      hi=lo;
    }
    step*=2.0;
  }
  if (!bracketed) return NA_REAL;
  for (int j=0;j<64 && hi-lo>1.0;++j) {
    const double mid=lo+floor((hi-lo)/2.0);
    const double fm=pSKELLAM(mid,dMu,dSigma2);
    if (!R_FINITE(fm)) return NA_REAL;
    if (fm>=dP) hi=mid; else lo=mid;
  }
  return hi-lo<=1.0 ? hi : NA_REAL;
}

double rSKELLAM(double dMu, double dSigma2) {
  if (!validSkellam(dMu,dSigma2)) return NA_REAL;

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  double dY1 = Rf_rpois(dMu1);
  double dY2 = Rf_rpois(dMu2);

  double dY = dY1 - dY2;

  return dY;

}

arma::vec mSKELLAM(double dMu, double dSigma2) {

  arma::vec vMoments(4); vMoments.fill(NA_REAL);
  if (!validSkellam(dMu,dSigma2)) return vMoments;

  vMoments(0) = dMu;
  vMoments(1) = dSigma2;
  vMoments(2) = dMu/pow(dSigma2, 1.5);
  vMoments(3) = 1/dSigma2 + 3.0;

  return vMoments;

}

arma::vec skellam_Score(double dY, arma::vec vTheta) {

  double dMu = vTheta(0);
  double dSigma2 = vTheta(1);
  if (!validSkellam(dMu,dSigma2) || !R_FINITE(dY) || dY != floor(dY)) {
    arma::vec out(2); out.fill(NA_REAL); return out;
  }

  double dMu1 = 0.5 * (dSigma2 + dMu);
  double dMu2 = 0.5 * (dSigma2 - dMu);

  const double z=2.0*sqrt(dMu1)*sqrt(dMu2);
  const double dBessel=Rf_bessel_i(z,fabs(dY),2);
  const double dDerivBessel=fabs(dY)/z*dBessel+Rf_bessel_i(z,fabs(dY)+1.0,2);

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

