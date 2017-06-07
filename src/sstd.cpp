#include <RcppArmadillo.h>
#include "Utils.h"

using namespace arma;
using namespace Rcpp;


//// This Functions come principally from the rugarch package of Ghalanos (2016) and have been slightly
//// modified to fit the GAS package


double xdt(const double x, const double nu)
{
  double a, b, pdf;
  a = Rf_gammafn((nu+1.0)/2.0)/sqrt(PI*nu);
  b = Rf_gammafn(nu/2.0)*pow((1.0+(x*x)/nu),((nu+1.0)/2.0));
  pdf = a/b;
  return pdf;
}

double dstdstd(const double x, const double nu)
{
  double pdf, s;
  if(nu<=2){
    pdf = 999;
  } else{
    s = sqrt(nu/(nu-2.0));
    pdf = s*xdt(x*s,nu);
  }
  return pdf;
}

double qstd(const double p, const double mu, const double sigma, const double nu)
{
  double s = sqrt(nu/(nu-2.0));
  double q = Rf_qt(p, nu, 1, 0) * sigma/s + mu;
  return( q );
}

double rstd(const double nu)
{
  double ans = 0;
  if(nu > 2.0)
  {
    double s = sqrt(nu/(nu-2));
    ans = Rf_rt(nu) * 1.0 / s;
  }
  return(ans);
}

double pstd(const double q, const double mu, const double sigma, const double nu)
{
  double s = sqrt(nu/(nu-2.0));
  double z = (q-mu)/sigma;
  double p = Rf_pt(z*s, nu, 1, 0);
  return( p );
}

double rsstd(const double xi, const double nu)
{
  double weight, z, rr, m1, mu, sigma, xx, ans;
  ans = 0.0;
  weight = xi / (xi + 1.0/xi);
  z = Rf_runif(-1.0 * weight, 1.0 - weight);
  xx = (z < 0)? 1.0/xi : xi;
  rr = -1.0 * fabs(rstd(nu))/xx * sign_C(z);
  m1 = 2.0 * sqrt(nu - 2.0) / (nu - 1.0) /  Rf_beta(0.5, 0.5 * nu);
  mu = m1 * (xi - 1.0/xi);
  sigma =  sqrt((1.0 - (m1 * m1)) * ((xi * xi) + 1.0/(xi * xi)) + 2 * (m1 * m1) - 1.0);
  ans =  (rr - mu ) / sigma;
  return(ans);
}

double dsstdstd(const double x, const double xi, const double nu)
{
  double mu, m1,beta, sigma, z, g,pdf,a,b, xxi;
  xxi=xi;
  a = 1.0/2.0;
  b = nu/2.0;
  beta = (Rf_gammafn(a)/Rf_gammafn(a+b))*Rf_gammafn(b);
  m1 = 2.0*sqrt(nu-2.0)/(nu-1.0)/beta;
  mu = m1*(xi-1.0/xi);
  sigma = sqrt( (1.0-pow(m1,2))*( pow(xi,2)+1.0/(pow(xi,2)) )+2.0*pow(m1,2)-1.0 );
  z = x*sigma + mu;
  if(z==0){
    xxi=1;
  }
  if(z<0){
    xxi = 1/xi;
  }
  g = 2.0/(xi+1.0/xi);
  pdf = g*dstdstd(z/xxi,nu)*sigma;
  return pdf;
}

double psstd(const double q, const double mu, const double sigma, const double xi, const double nu)
{
  double qx = (q-mu)/sigma;
  double m1 = 2.0 * sqrt(nu-2.0) / (nu-1.0) / Rf_beta(0.5, nu/2.0);
  double mux = m1*(xi-1.0/xi);
  double sig =  sqrt((1-m1*m1)*(xi*xi+1/(xi*xi)) + 2*m1*m1 - 1);
  double z = qx*sig+mux;
  double Xi = (z<0)?1.0/xi:xi;
  double g = 2.0 / (xi + 1.0/xi);
  double p = Heaviside(z, 0) - signum(z) * g * Xi * pstd(-fabs(z)/Xi, 0, 1, nu);
  return( p );
}


double qsstd(const double p, const double xi, const double nu)
{
  double m1 = 2.0 * sqrt(nu-2.0) / (nu-1.0) / Rf_beta(0.5, nu/2.0);
  double mu = m1*(xi-1.0/xi);
  double sigma =  sqrt((1-m1*m1)*(xi*xi+1/(xi*xi)) + 2*m1*m1 - 1);
  double g = 2.0 / (xi + 1.0/xi);
  double z = p-0.5;
  double Xi = (z<0)?1.0/xi:xi;
  double tmp = (Heaviside(z, 0) - signum(z)*p)/(g*Xi);
  double q = (-signum(z)*qstd(tmp, 0, 1, nu)*Xi - mu)/sigma;
  return( q );
}


double sstdskew(double dXi, double dNu){
  // Theoretical moments based on bijection betweeen Fernandez and Steel verions
  // and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
  double m3 = 0.0;
  if (dNu > 2.0) {
    double eta = dNu;
    double k2  = pow(dXi, 2.0);
    double lda = (k2-1.0)/(k2+1.0);
    double ep1 = (eta+1.0)/2.0;
    double lnc = Rf_lgammafn(ep1) - Rf_lgammafn(eta/2.0) -0.5*log( M_PI*(eta-2.0));
    double c   = exp(lnc);
    double a   = 4.0*lda*c*(eta-2.0)/(eta-1.0);
    double b   = pow(1.0+3.0*pow(lda,2.0)-pow(a,2.0), 0.5);
    double my2 = 1.0+3.0*pow(lda,2.0);
    double my3 = 16.0*c*lda*(1.0+pow(lda,2.0))*(  pow(eta-2.0,2.0))/((eta-1.0)*(eta-3.0));
    // double my4 = 3.0*(eta-2.0)*(1.0+10.0*pow(lda,2.0)+5.0*pow(lda,4.0) )/(eta-4.0);
    m3  = (my3-3.0*a*my2+2.0*pow(a,3.0) )/pow(b,3.0);
  } else{
    m3 = NA_REAL;
  }
  return m3;
}


double sstdexkurt( double dXi, double dNu )
{
  // Theoretical moments based on bijection betweeen Fernandez and Steel verions
  // and Hansen's Generalized Skew-T (credit due to Michael Rockinger)
  double m4 = 0.0;
  if(dNu > 4.0){
    double eta  = dNu;
    double k2   = pow(dXi,2.0);
    double lda  = (k2-1.0)/(k2+1.0);
    double ep1 = (eta+1.0)/2.0;
    double lnc = Rf_lgammafn(ep1) - Rf_lgammafn(eta/2.0) -0.5*log( M_PI*(eta-2.0));
    double c   = exp(lnc);
    double a   = 4.0*lda*c*(eta-2.0)/(eta-1.0);
    double b   = pow(1.0+3.0*pow(lda,2.0)-pow(a,2.0),0.5 );
    double my2 = 1.0+3.0*pow(lda,2.0);
    double my3 = 16.0*c*lda*(1.0+pow(lda,2.0) )*(pow(eta-2.0,2.0))/((eta-1.0)*(eta-3.0));
    double my4 = 3.0*(eta-2.0)*(1.0+10.0*pow(lda,2.0)+5.0*pow(lda,4.0) )/(eta-4.0);
    // double m3  = (my3-3.0*a*my2+2.0*pow(a,3.0) )/pow(b,3.0);
    m4  = -3.0 + (my4-4.0*a*my3+6.0*pow(a,2.0)*my2-3.0*pow(a,4.0) )/pow(b,4.0);
  } else{
    m4 = NA_REAL;
  }
  return m4;
}

double dSSTD(double dY, double dMu, double dSigma, double dXi, double dNu, bool bLog = false ){
  double dZ = (dY - dMu)/dSigma;
  double dPDF = dsstdstd(dZ, dXi, dNu)/dSigma;

  if(bLog) dPDF = log(dPDF);

  return dPDF;
}

double pSSTD(double dY, double dMu, double dSigma, double dXi, double dNu){

  double dP = psstd(dY, dMu, dSigma, dXi, dNu);
  return dP;
}

double rSSTD(double dMu, double dSigma, double dXi, double dNu){

  double dY = dMu + rsstd(dXi, dNu)*dSigma;
  return dY;
}

double qSSTD(double dP, double dMu, double dSigma, double dXi, double dNu){

  double dQ = qsstd(dP, dXi, dNu)*dSigma + dMu;

  return dQ;
}

arma::vec mSSTD(double dMu, double dSigma, double dXi, double dNu){
  arma::vec vMoments(4);

  vMoments(0) = dMu;
  vMoments(1) = pow(dSigma, 2.0);
  vMoments(2) = sstdskew( dXi, dNu );
  vMoments(3) = sstdexkurt( dXi, dNu ) + 3.0;

  return vMoments;
}

arma::vec sstd_Score(double dY, arma::vec vTheta){

  double dMu    = vTheta(0);
  double dSigma = vTheta(1);
  double dXi    = vTheta(2);
  double dNu    = vTheta(3);

  // double dMu1 = 2.0*pow(dNu - 2.0, 0.5)/(dNu - 1.0)  * Rf_gammafn(0.5*(dNu + 1.0))/(Rf_gammafn(0.5*dNu)*Rf_gammafn(0.5));

  double dLogMu1 = log(2.0) + 0.5 * log(dNu - 2.0) - log(dNu - 1.0) + Rf_lgammafn(0.5 * (dNu + 1.0)) - Rf_lgammafn(0.5 * dNu) - Rf_lgammafn(0.5);
  double dMu1    = exp(dLogMu1);

  double dMu_tilde = dMu1*(dXi - 1.0/dXi);
  double dSigma_tilde = pow( (1.0 - pow(dMu1, 2.0))*(pow(dXi, 2.0) + pow(dXi, -2.0)) + 2.0*pow(dMu1, 2.0) -1.0, 0.5);
  double dLogSigma_tilde = 0.5 * log((1.0 - pow(dMu1, 2.0))*(pow(dXi, 2.0) + pow(dXi, -2.0)) + 2.0*pow(dMu1, 2.0) -1.0);
  double dZ = (dY - dMu)/dSigma * dSigma_tilde + dMu_tilde;

  double dXi_star  = dXi;
  double ddXi_star2 = 2.0 * dXi;

  if (dZ == 0.0){
    dXi_star = 1.0;
    ddXi_star2 = 0;
  }
  if (dZ < 0.0) {
    dXi_star = 1.0/dXi;
    ddXi_star2 = -2.0/pow(dXi, 3.0);
  }

  double dC = 1.0 + pow(dZ, 2.0)/(pow(dXi_star, 2.0) * (dNu - 2.0));

  double dA = 0.5*Rf_gammafn(0.5*(dNu + 1.0)) * ( pow(dNu - 2.0, -0.5) + pow(dNu - 2.0, 0.5)*Rf_digamma(0.5*(dNu + 1.0)));
  double dB = Rf_gammafn(0.5 * dNu) + (dNu - 1.0)*Rf_digamma(0.5*dNu)*Rf_gammafn(0.5*dNu) * 0.5;
  double ddMu1_nu = (dA*(dNu - 1.0)*Rf_gammafn(0.5*dNu) - dB*pow(dNu-2.0, 0.5)*Rf_gammafn(0.5*(dNu+1.0)))/pow( (dNu - 1.0)*Rf_gammafn(0.5*dNu), 2.0)*2.0/Rf_gammafn(0.5);


  double ddSigma_tilde_nu = (ddMu1_nu*( -(pow(dXi, 2.0) + pow(dXi, -2.0)) + 2.0 )*dMu1)/dSigma_tilde;

  double ddL_nu = ddSigma_tilde_nu*(dY - dMu)/dSigma + (dXi - 1.0/dXi)*ddMu1_nu;

  // score

  double ddMu    = (dNu + 1.0)/dC  * dSigma_tilde/dSigma * dZ/( pow(dXi_star, 2.0)*(dNu - 2.0) );

  // old one
  // double ddSigma = -1.0/dSigma + (dNu + 1.0)/dC * dSigma_tilde/pow(dSigma, 2.0) * dZ*(dY - dMu)/( pow(dXi_star, 2.0)*(dNu - 2.0) );

  double dQ = -1.0/dSigma;

  double dFoo1_W = 2.0 * log(dXi_star) + log(dNu - 2.0);
  double dFoo2_W = 2.0 * log(abs3(dZ));
  double dW = -3.0 * log(dSigma) + 2.0 * dLogSigma_tilde + 2.0 * log(abs3(dY - dMu)) + log(dNu + 1.0) - LogSum(dFoo1_W, dFoo2_W);

  double dFoo_E = dMu_tilde*(dY - dMu);

  double dSgn = 0.0;

  if (dFoo_E < 0) {

    dSgn = -1.0;

  } else {

    dSgn = +1.0;

  }

  double dE = log(dNu + 1.0) + dLogSigma_tilde + log(abs3(dMu_tilde * (dY - dMu))) - 2.0 * log(dSigma) - LogSum(dFoo1_W, dFoo2_W);

  double ddSigma = dQ + exp(dW) + dSgn * exp(dE);


  double ddNu     = ddSigma_tilde_nu/dSigma_tilde + 0.5*(1.0/dNu - 1.0/(dNu-2.0)) + Rf_digamma(0.5*(dNu+1.0))*0.5 -
                    0.5*Rf_digamma(0.5*dNu) - 1.0/(2.0 * dNu) - 0.5*(log(dC) + dZ*(dNu + 1.0)*(2.0*(dNu-2.0)*ddL_nu - dZ)/(dC*pow(dXi_star*(dNu-2.0),2.0)));


  double ddLogSigma_tilde_xi = (1.0 + pow(dMu1, 2.0))*(dXi - pow(dXi, -3.0))/((1.0 + pow(dMu1, 2.0))*(pow(dXi, 2.0) + pow(dXi, -2.0)) + 2.0*pow(dMu1, 2.0) - 1.0 );
  double ddSigma_tilde_xi = ddLogSigma_tilde_xi * dSigma_tilde;

  double ddLogG_xi = -(1.0 - 1.0/pow(dXi, 2.0))/(dXi + 1.0/dXi);

  double dL1 = pow(dZ, 2.0);
  double dL2 = (dNu - 2.0) * pow(dXi_star, 2.0);

  double ddL1 = 2.0 * dZ * ( (dY - dMu)/dSigma *ddSigma_tilde_xi + dMu1*(1.0 + 1.0/pow(dXi, 2.0)) );
  double ddL2 = (dNu - 2.0)*ddXi_star2;

  double ddXi    = ddLogG_xi + ddLogSigma_tilde_xi - 0.5 * (dNu + 1.0) *  1.0 / dC * ( ddL1*dL2 -  dL1* ddL2 )/pow(dL2, 2.0);

  arma::vec vScore(4);

  vScore(0) = ddMu;
  vScore(1) = ddSigma;
  vScore(2) = ddXi;
  vScore(3) = ddNu;

  return vScore;

}

arma::mat sstd_IM(arma::vec vTheta){

  arma::mat mIM = eye(4,4);

  return mIM;

}

