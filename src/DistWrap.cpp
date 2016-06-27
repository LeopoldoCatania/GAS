#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"
#include "mvnorm.h"
#include "mvt.h"


using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double ddist_univ(double dY, arma::vec vTheta, std::string Dist, bool bLog){
  double dLPDF=0.0;
  if(Dist == "norm") dLPDF = dNORM(dY, vTheta(0), vTheta(1), bLog );
  if(Dist == "std")  dLPDF = dSTD(dY, vTheta(0), vTheta(1), vTheta(2), bLog);
  if(Dist == "ast")  dLPDF = dAST(dY, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4), bLog);
  if(Dist == "ast1") dLPDF = dAST(dY, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3), bLog);
  return dLPDF;
}

//[[Rcpp::export]]
double ddist_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist, bool bLog){
  double dLPDF=0.0;
  if(Dist == "mvnorm") dLPDF = dmvnorm_ThetaParam(vY, vTheta, iN, bLog);
  if(Dist == "mvt")    dLPDF = dmvt_ThetaParam(vY, vTheta, iN, bLog);
  return dLPDF;
}

//[[Rcpp::export]]
double rdist_univ(arma::vec vTheta, std::string Dist){
  double dY = 0.0;
  if(Dist == "norm") dY = vTheta(0) + pow(vTheta(1),0.5)*Rf_rnorm(0.0,1.0);
  if(Dist == "std")  dY = vTheta(0) + pow(vTheta(1),0.5)*Rf_rt(vTheta(2));
  if(Dist == "ast")  dY = rAST(vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4));
  if(Dist == "ast1") dY = rAST(vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3));
  return dY;
}

//[[Rcpp::export]]
arma::vec rdist_multi(arma::vec vTheta, int iN,std::string Dist){
  arma::vec vY(iN);
  if(Dist == "mvnorm") vY = arma::vectorise(rmvnorm_ThetaParam(vTheta,iN, 1));
  if(Dist == "mvt")    vY = arma::vectorise(rmvt_ThetaParam(vTheta,iN, 1));
  return vY;
}

//[[Rcpp::export]]
double pdist_univ(double dQ, arma::vec vTheta, std::string Dist){
  double dP=0.0;
  if(Dist == "norm") dP = Rf_pnorm5(dQ, vTheta(0), pow(vTheta(1),2.0), 1,0);
  if(Dist == "std")  dP = pSTD(dQ, vTheta(0), vTheta(1), vTheta(2));
  if(Dist == "ast")  dP = pAST(dQ, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4));
  if(Dist == "ast1") dP = pAST(dQ, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3));
  return dP;
}
//[[Rcpp::export]]
double qdist_univ(double dP, arma::vec vTheta, std::string Dist){
  double dQ=0.0;
  if(Dist == "norm") dQ = Rf_qnorm5(dP, vTheta(0), pow(vTheta(1),2.0), 1,0);
  if(Dist == "std")  dQ = qSTD(dP, vTheta(0), vTheta(1), vTheta(2));
  if(Dist == "ast")  dQ = qAST(dP, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4));
  if(Dist == "ast1") dQ = qAST(dP, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3));
  return dQ;
}
