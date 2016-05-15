#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
double ddist_univ(double dY, arma::vec vTheta, std::string dist, bool bLog){
  double dLPDF=0.0;
  if(dist == "norm") dLPDF = dNORM(dY, vTheta(0), vTheta(1), bLog );
  if(dist == "std")  dLPDF = dSTD(dY, vTheta(0), vTheta(1), vTheta(2), bLog);
  if(dist == "ast")  dLPDF = dAST(dY, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(4), bLog);
  if(dist == "ast1") dLPDF = dAST(dY, vTheta(0), vTheta(1), vTheta(2),vTheta(3),vTheta(3), bLog);
  return dLPDF;
}
