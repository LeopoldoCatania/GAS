#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::vec Score_univ(double dY, arma::vec vTheta,std::string dist){
  arma::vec vScore;
  if(dist == "norm") vScore = norm_Score(dY,vTheta);
  if(dist == "std")  vScore = std_Score(dY,vTheta);
  if(dist == "ast")  vScore = ast_Score(dY,vTheta);
  if(dist == "ast1") vScore = ast1_Score(dY,vTheta);
  return vScore;
}
