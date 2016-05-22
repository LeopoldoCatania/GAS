#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double dSTD(double dY, double dMu, double dPhi2 , double dNu, bool bLog=false) {

  double dLPDF = Rf_lgammafn((dNu+1.0)/2.0) -  Rf_lgammafn(dNu/2.0) - 0.5*log(dPhi2) - 0.5*log(M_PI*dNu) - (dNu+1.0)/2.0 * log(1.0 + pow(dY-dMu,2.0)/(dNu*dPhi2));

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;

}

arma::vec std_Score(double dY, arma::vec vTheta){

  double dMu   = vTheta(0);
  double dPhi2 = vTheta(1);
  double dNu   = vTheta(2);

  double dNu_s = R::digamma((dNu + 1.0)/2.0)*0.5 - R::digamma(dNu/2.0)*0.5 - 1.0/(2.0*dNu) - 0.5*log(1.0+ pow(dY-dMu,2.0)/(dNu*dPhi2) )+
    0.5*(dNu+1.0)*pow(dY-dMu,2.0)/(dPhi2*pow(dNu,2.0) )/(1.0 +  pow(dY-dMu,2.0)/(dNu*dPhi2) );

  double dPhi2_s = -1.0/(2.0*dPhi2) + (dNu+1.0)*pow(dY-dMu,2.0)/(2.0*pow(dPhi2,2.0)*dNu*(1.0 +  pow(dY-dMu,2.0)/(dNu*dPhi2)));

  double dMu_s = (dNu+1.0)*(dY-dMu)/(dNu*dPhi2 + pow(dY-dMu,2.0));

  arma::vec vScore(3);

  vScore(0)=dMu_s;
  vScore(1)=dPhi2_s;
  vScore(2)=dNu_s;

  return vScore;

}
arma::mat std_IM( arma::vec vTheta){

  double dPhi2 = vTheta(1);
  double dNu   = vTheta(2);

  arma::mat mIM=zeros(3,3);

  double uu = (dNu+1.0)/(dPhi2*(dNu+3.0));
  double dd = dNu/(2.0*pow(dPhi2,2.0)*(dNu+3.0));
  double tt = 0.5*( 0.5* R::trigamma(0.5*dNu) - 0.5* R::trigamma( 0.5*(dNu+1.0) ) - (dNu+5.0)/(dNu*(dNu+3.0)*(dNu+1.0)));
  double td = -1.0/(dPhi2*(dNu+3.0)*(dNu+1.0));

  mIM(0,0)=uu;
  mIM(1,1)=dd;
  mIM(2,2)=tt;
  mIM(2,1)=td;
  mIM(1,2)=td;

  return mIM;
}
