#include <RcppArmadillo.h>
#include "Utils.h"

using namespace Rcpp;
using namespace arma;

double Kast(double dX){

  double dOut = Rf_lgammafn((dX+1.0)/2.0) - 0.5*log(M_PI*dX) - Rf_lgammafn(dX*0.5);

  return exp(dOut);
}
double dAST(double dY, double dMu, double dSigma, double dAlpha, double dNu1,double dNu2 , bool bLog=false) {

  double dLPDF=0.0;

  if(dY<=dMu)   dLPDF= -log(dSigma) - (dNu1+1.0)/2.0 * log( 1.0 + pow(  (dY-dMu)/(2.0*dAlpha*dSigma*Kast(dNu1)) ,2.0) /dNu1) ;
  if(dY>dMu)    dLPDF= -log(dSigma) -(dNu2+1.0)/2.0 *log(1.0 + pow(  (dY-dMu)/(2.0*(1.0-dAlpha)*dSigma*Kast(dNu2)) ,2.0)/dNu2) ;

  if(!bLog) dLPDF=exp(dLPDF);

  return dLPDF;
}
double dKast(double dX){

  double out=(0.5*Rf_gammafn((dX+1)/2)*Rf_digamma((dX+1)/2) *sqrt(M_PI*dX)*Rf_gammafn(dX/2) - Rf_gammafn((dX+1)/2)* ( sqrt(M_PI)*Rf_gammafn(dX/2)*(dX*Rf_digamma(dX/2)+1)/(2*sqrt(dX))
  )  )/(M_PI*dX*pow(Rf_gammafn(0.5*dX),2));

  return out;
}
arma::vec ast_Score(double dY, arma::vec vTheta){

  double mu=vTheta(0);
  double sigma=vTheta(1);
  double alpha = vTheta(2);
  double df1=vTheta(3);
  double df2=vTheta(4);

  double K1=Kast(df1);
  double K2=Kast(df2);

  double dK1= dKast(df1);
  double dK2= dKast(df2);

  double A1= 1.0 + pow((dY-mu)/(2.0*alpha*sigma*K1),2.0)/df1;
  double A2= 1.0 + pow((dY-mu)/(2.0*(1.0-alpha)*sigma*K2),2.0)/df2;

  double dmu=0;
  double dsigma=0;
  double dalpha=0;
  double ddf1=0;
  double ddf2=0;

  if(dY<=mu){

    dmu =  (df1+1.0)/A1 * (dY-mu)/(df1*pow(2.0*alpha*sigma*K1,2.0));
    dsigma = -1.0/sigma + (df1+1.0)/A1 * (pow((dY-mu)/(2.0*alpha*K1),2.0))/(df1*pow(sigma,3.0));
    dalpha = (df1+1.0)/A1 * pow((dY-mu)/(2.0*sigma*K1),2.0)/(df1*pow(alpha,3.0));
    ddf1= -(  0.5*log(A1) + (df1+1.0)/(2.0*A1)*((-1.0/pow(df1,2.0) * pow((dY-mu)/(2.0*alpha*sigma*K1),2.0) + 1.0/df1 * pow((dY-mu)/(2.0*alpha*sigma),2.0)/pow(K1,3.0)
                                                   *(-2.0*dK1))  ));

  }else{

    dmu =  (df2+1.0)/A2 * (dY-mu)/(df2*pow(2.0*(1.0-alpha)*sigma*K2,2.0));
    dsigma = -1.0/sigma + (df2+1.0)/A2 * (pow((dY-mu)/(2.0*(1.0-alpha)*K2),2.0))/(df2*pow(sigma,3.0));
    dalpha = -(df2+1.0)/A2 * pow((dY-mu)/(2.0*sigma*K2),2.0)/(df2*pow(1.0-alpha,3.0));
    ddf2= -(  0.5*log(A2) + (df2+1.0)/(2.0*A2)*((-1.0/pow(df2,2.0) * pow((dY-mu)/(2.0*(1.0-alpha)*sigma*K2),2.0) +
      1.0/df2 * pow((dY-mu)/(2.0*(1.0-alpha)*sigma),2.0)/pow(K2,3.0)  *(-2.0*dK2))  ));
  }


  arma::vec vScore(5);

  vScore(0)=dmu;
  vScore(1)=dsigma;
  vScore(2)=dalpha;
  vScore(3)=ddf1;
  vScore(4)=ddf2;

  return vScore;


}
arma::vec ast1_Score(double dY, arma::vec vTheta){

  double mu=vTheta(0);
  double sigma=vTheta(1);
  double alpha = vTheta(2);
  double df1=vTheta(3);

  double K1=Kast(df1);

  double dK1= dKast(df1);

  double A1= 1.0 + pow((dY-mu)/(2.0*alpha*sigma*K1),2.0)/df1;
  double A2= 1.0 + pow((dY-mu)/(2.0*(1.0-alpha)*sigma*K1),2.0)/df1;

  double dmu=0;
  double dsigma=0;
  double dalpha=0;
  double ddf1=0;

  if(dY<=mu){

    dmu =  (df1+1.0)/A1 * (dY-mu)/(df1*pow(2.0*alpha*sigma*K1,2.0));
    dsigma = -1.0/sigma + (df1+1.0)/A1 * (pow((dY-mu)/(2.0*alpha*K1),2.0))/(df1*pow(sigma,3.0));
    dalpha = (df1+1.0)/A1 * pow((dY-mu)/(2.0*sigma*K1),2.0)/(df1*pow(alpha,3.0));
    ddf1= -(  0.5*log(A1) + (df1+1.0)/(2.0*A1)*((-1.0/pow(df1,2.0) * pow((dY-mu)/(2.0*alpha*sigma*K1),2.0) +
      1.0/df1 * pow((dY-mu)/(2.0*alpha*sigma),2.0)/pow(K1,3.0)  *(-2.0*dK1))  ));

  }else{

    dmu =  (df1+1.0)/A2 * (dY-mu)/(df1*pow(2.0*(1.0-alpha)*sigma*K1,2.0));
    dsigma = -1.0/sigma + (df1+1.0)/A2 * (pow((dY-mu)/(2.0*(1.0-alpha)*K1),2.0))/(df1*pow(sigma,3.0));
    dalpha = -(df1+1.0)/A2 * pow((dY-mu)/(2.0*sigma*K1),2.0)/(df1*pow(1.0-alpha,3.0));
    ddf1= -(  0.5*log(A2) + (df1+1.0)/(2.0*A2)*((-1.0/pow(df1,2.0) * pow((dY-mu)/(2.0*(1.0-alpha)*sigma*K1),2.0) +
      1.0/df1 * pow((dY-mu)/(2.0*(1.0-alpha)*sigma),2.0)/pow(K1,3.0)  *(-2.0*dK1))  ));
  }


  arma::vec vScore(4);

  vScore(0)=dmu;
  vScore(1)=dsigma;
  vScore(2)=dalpha;
  vScore(3)=ddf1;

  return vScore;


}
arma::mat ast_IM(arma::vec vTheta){

  double sigma = vTheta(1);
  double alpha = vTheta(2);
  double df1   = vTheta(3);
  double df2   = vTheta(4);

  arma::mat mIM=zeros(5,5);

  double DV1 = Rf_digamma((df1+1.0)/2.0) - Rf_digamma(df1/2.0);
  double DV2 = Rf_digamma((df2+1.0)/2.0) - Rf_digamma(df2/2.0);

  double dDV1 = Rf_trigamma((df1+1.0)/2.0)*0.5 - Rf_trigamma(df1/2.0)*0.5;
  double dDV2 = Rf_trigamma((df2+1.0)/2.0)*0.5 - Rf_trigamma(df2/2.0)*0.5;

  double KAST1 = Kast(df1);
  double KAST2 = Kast(df2);

  double uu = 3.0*( (df1+1.0)/(alpha*(df1+3.0)) + (df2+1.0)/((1.0-alpha)*(df2+3.0))  );
  double ud =- 1.0/(df1+1.0)  + (df1*DV1)/(df1+3.0);
  double ut = 1.0/(df2+1.0) - df2*DV2/(df2+3.0) ;
  double uq = -2.0*uu/(3.0*sigma);
  double uc = 2.0/sigma * (df1/(df1+3.0) - df2/(df2+3.0));
  double dd = alpha / 2.0 * (df1*pow(DV1,2.0)/(df1+3.0) - 2.0*DV1/(df1+1.0) -  dDV1 ) ;
  double dc = alpha*ud/sigma;
  double dq =  (  1.0/(df1+1.0) - (df1+1.0)/(df1+3.0) * DV1 ) /sigma;
  double tq = -(  1.0/(df2+1.0) - (df2+1.0)/(df2+3.0) * DV2 )/sigma;
  double tt = 0.5*(1.0-alpha)*( (df2*pow(DV2,2.0))/(df2+3.0) - 2.0*DV2/(df2+1.0) - dDV2  );
  double tc = -(1.0-alpha)*ut/sigma;
  double qq = ( (df1+1.0)/(alpha*(df1+3.0)*pow(KAST1,2.0)) + (df2+1.0)/((1.0-alpha)*(df2+3.0)*pow(KAST2,2.0)))/(4.0*pow(sigma,2.0));
  double qc = -2.0*uc/(3.0*sigma);
  double cc = 2.0*(alpha*df1/(df1+3.0) + (1.0-alpha)*df2/(df2+3.0)  )/pow(sigma,2.0);


  mIM(0,0)=qq;
  mIM(0,1)=qc;
  mIM(1,0)=qc;
  mIM(0,2)=uq;
  mIM(2,0)=uq;
  mIM(0,3)=dq;
  mIM(3,0)=dq;
  mIM(0,4)=tq;
  mIM(4,0)=tq;
  mIM(1,1)=cc;
  mIM(1,2)=uc;
  mIM(2,1)=uc;
  mIM(1,3)=dc;
  mIM(3,1)=dc;
  mIM(1,4)=tc;
  mIM(4,1)=tc;
  mIM(2,2)=uu;
  mIM(2,3)=ud;
  mIM(3,2)=ud;
  mIM(2,4)=ut;
  mIM(4,2)=ut;
  mIM(3,3)=dd;
  mIM(4,4)=tt;

  return mIM;

}
arma::mat ast1_IM(arma::vec vTheta){

  double sigma = vTheta(1);
  double alpha = vTheta(2);
  double df1   = vTheta(3);

  double DV1  = Rf_digamma((df1+1.0)/2.0) - Rf_digamma(df1/2.0);
  double dDV1 = Rf_trigamma((df1+1.0)/2.0)*0.5 - Rf_trigamma(df1/2.0)*0.5;

  double KAST1 = Kast(df1);

  double uu = 3.0*( (df1+1.0)/(alpha*(df1+3.0)) + (df1+1.0)/((1.0-alpha)*(df1+3.0))  );
  double ud =- 1.0/(df1+1.0)  + (df1*DV1)/(df1+3.0);
  double ut = 1.0/(df1+1.0) - df1*DV1/(df1+3.0) ;
  double uq = -2.0*uu/(3.0*sigma);
  double uc = 2.0/sigma * (df1/(df1+3.0) - df1/(df1+3.0));
  double dd = alpha / 2.0 * (df1*pow(DV1,2.0)/(df1+3.0) - 2.0*DV1/(df1+1.0) -  dDV1 ) ;
  double dc = alpha*ud/sigma;
  double dq =  (  1.0/(df1+1.0) - (df1+1.0)/(df1+3.0) * DV1 ) /sigma;
  double tq = -(  1.0/(df1+1.0) - (df1+1.0)/(df1+3.0) * DV1 )/sigma;
  double tt = 0.5*(1.0-alpha)*( (df1*pow(DV1,2.0))/(df1+3.0) - 2.0*DV1/(df1+1.0) - dDV1  );
  double tc = -(1.0-alpha)*ut/sigma;
  double qq = ( (df1+1.0)/(alpha*(df1+3.0)*pow(KAST1,2.0)) + (df1+1.0)/((1.0-alpha)*(df1+3.0)*pow(KAST1,2.0)))/(4.0*pow(sigma,2.0));
  double qc = -2.0*uc/(3.0*sigma);
  double cc = 2.0*(alpha*df1/(df1+3.0) + (1.0-alpha)*df1/(df1+3.0)  )/pow(sigma,2.0);

  arma::mat mIM=zeros(4,4);

  mIM(0,0)=qq;
  mIM(0,1)=qc;
  mIM(1,0)=qc;
  mIM(0,2)=uq;
  mIM(2,0)=uq;
  mIM(0,3)=dq+tq;
  mIM(3,0)=dq+tq;
  mIM(1,1)=cc;
  mIM(1,2)=uc;
  mIM(2,1)=uc;
  mIM(1,3)=dc+tc;
  mIM(3,1)=dc+tc;
  mIM(2,2)=uu;
  mIM(2,3)=ud+ut;
  mIM(3,2)=ud+ut;
  mIM(3,3)=dd+tt;

  return mIM;


}

double rAST(double dMu, double dSigma,double dAlpha, double dNu1, double dNu2){
  double dU  = Rf_runif(0.0,1.0);
  double dT1 = Rf_rt(dNu1);
  double dT2 = Rf_rt(dNu2);

  double K1=Kast(dNu1);
  double K2=Kast(dNu2);

  double dAlpha_star = dAlpha * K1 / (dAlpha*K1+(1-dAlpha)*K2);

  double dSign = sign_C(dU-dAlpha);

  double dY = dMu + dSigma * (dAlpha_star * abs3(dT1)* (dSign-1) +
                              (1-dAlpha_star)*abs3(dT2) * (dSign+1) );

  return dY;
}
double pAST(double dY, double dMu, double dSigma, double dAlpha, double dNu1, double dNu2) {

  double dLogK = 0.0;
  double dZ = 0.0;
  double dCdf = 0.0;
  double foo = 0.0;

  if (dY <= dMu){
    dLogK = log(tgamma((dNu1 + 1.0) / 2.0)) - log(tgamma(dNu1 / 2.0)) - 0.5 * log(M_PI * dNu1);
    dZ    = (dY - dMu) / (2.0 * dAlpha * (dSigma) * exp(dLogK));
    foo = R::pt(dZ, dNu1, 1, 0);
    dCdf  = 2.0 * dAlpha * foo;
  }else{
    dLogK = log(tgamma((dNu2 + 1.0) / 2.0)) - log(tgamma(dNu2 / 2.0)) - 0.5 * log(M_PI * dNu2);
    dZ    = (dY - dMu) / (2.0 * (1.0 - dAlpha) * (dSigma) * exp(dLogK));
    foo = R::pt(dZ, dNu2, 1, 0);
    dCdf  = dAlpha + 2.0 * (1.0 - dAlpha) * (foo - 0.5);
  }
  return dCdf;
}

arma::vec mAST(double dMu, double dSigma, double dAlpha, double dNu1, double dNu2){
  arma::vec vMoments(4);

  double dK1 = Kast(dNu1);
  double dK2 = Kast(dNu2);

  double dAlpha_star = dAlpha*dK1/(dAlpha*dK1+(1.0-dAlpha)*dK2);
  double dB          = dAlpha*dK1 + (1.0-dAlpha)*dK2;

  double dEY = 4.0 * dB* (-pow(dAlpha_star,2.0)*(dNu1/(dNu1-1.0))+ pow(1.0-dAlpha_star,2.0)*dNu2/(dNu2 - 1.0));

  vMoments(0) = dMu + dSigma * dEY;


  vMoments(1) = (4.0 * (dAlpha * pow(dAlpha_star,2.0)*dNu1/(dNu1-2.0) +
                (1.0-dAlpha)*pow(1.0-dAlpha_star,2.0)*dNu2/(dNu2-2.0)) -
                16.0*pow(dB,2.0)*pow(-pow(dAlpha_star,2.0)*dNu1/(dNu1-1.0) +
                pow(1.0 - dAlpha_star,2.0)*dNu2/(dNu2-1.0),2.0) ) * pow(dSigma, 2.0);

  // TODO
  vMoments(2) = 0.0;
  vMoments(3) = 0.0;
  return vMoments;
}

double qAST(double dP, double dMu, double dSigma, double dAlpha, double dNu1, double dNu2,
                        double lower=-150, double upper=150, int maxiter=1e4, double eps=1e-7) {
  double a=lower;
  double b=upper;

  double x=lower;
  double x1=upper;
  int iter = 1;
  double fa,fx;
  //check
  fa=pAST(a, dMu, dSigma, dAlpha, dNu1, dNu2) - dP;
  fx=pAST(x1, dMu, dSigma, dAlpha, dNu1, dNu2) - dP;

  if(fa*fx>0){
    Rprintf("Bisection Error: upper and lower function evaluations have same sign");
    return NA_LOGICAL;
  }

  do
  {
    fa=pAST(a, dMu, dSigma, dAlpha, dNu1, dNu2) - dP;
    fx=pAST(x, dMu, dSigma, dAlpha, dNu1, dNu2) - dP;

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
