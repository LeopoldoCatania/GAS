#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"
#include "mvnorm.h"
#include "mvt.h"


using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec Score_univ(double dY, arma::vec vTheta,std::string Dist){
  arma::vec vScore;
  if(Dist == "norm") vScore = norm_Score(dY,vTheta);
  if(Dist == "std")  vScore = std_Score(dY,vTheta);
  if(Dist == "ast")  vScore = ast_Score(dY,vTheta);
  if(Dist == "ast1") vScore = ast1_Score(dY,vTheta);
  return vScore;
}
//[[Rcpp::export]]
arma::vec Score_multi(arma::vec vY, arma::vec vTheta, int iN,std::string Dist){
  arma::vec vScore;
  if(Dist == "mvnorm") vScore = mvnorm_Score(vY, vTheta, iN);
  if(Dist == "mvt")    vScore = mvt_Score(vY, vTheta, iN);
  return vScore;
}
