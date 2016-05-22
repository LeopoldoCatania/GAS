#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"

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
double rdist_univ(arma::vec vTheta, std::string Dist){
  double dY = 0.0;
  if(Dist == "norm") dY = vTheta(0) + pow(vTheta(1),0.5)*Rf_rnorm(0.0,1.0);
  if(Dist == "std")  dY = vTheta(0) + pow(vTheta(1),0.5)*Rf_rt(vTheta(2));
  if(Dist == "ast")  dY = rAST(vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4));
  if(Dist == "ast1") dY = rAST(vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3));
  return dY;
}

