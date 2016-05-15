#include <RcppArmadillo.h>
#include "norm.h"
#include "std.h"
#include "ast.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat IM_univ(arma::vec vTheta,std::string Dist){
  arma::mat mIM;
  if(Dist == "norm") mIM = norm_IM(vTheta);
  if(Dist == "std")  mIM = std_IM(vTheta);
  if(Dist == "ast")  mIM = ast_IM(vTheta);
  if(Dist == "ast1") mIM = ast1_IM(vTheta);
  return mIM;
}

